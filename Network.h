#include "Eigen/Eigen"
#include "induct.h"

#ifndef _CMATH_
#include <cmath>
#endif 

#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#include<future>
#include<thread>
#include<vector>
//ya,I know that looks stupid, but it SHUANG SI LE.

#define __GII__ induct_network(i,i).real()
#define __GIJ__ induct_network(i,j).real()
#define __BII__ induct_network(i,i).imag()
#define __BIJ__ induct_network(i,j).imag()
#define __EI__ node[i].e
#define __EJ__ node[j].e
#define __FI__ node[i].f
#define __FJ__ node[j].f
constexpr auto LOCAL_INF = 114514;


using namespace Eigen;
using std::complex;
using std::async;
using std::future;
using std::vector;
/*
*class Network.
*OOP
*/
//!this flow calculate library obey "first PQ then PV" arrange.
//!balance node must be the last one
class Network
{
private:
    unsigned int cpu_core_number=1u;
public:
    /*data*/
    enum Node_type
    {
        pv,
        pq,
        balance
    };
    struct node_arg
    {
        double p = 0.0;
        double q = 0.0;
        double e = 1.0;
        double f = 0.0;
        double r;
        double angle;
        complex<double> u;
        double e_balance()
        {
            return r * cos(angle);
        }
        double f_balance()
        {
            return r * sin(angle);
        }
        Node_type type;
    };

    int pv_node_number = 0;
    int pq_node_number = 0;
    int matrix_length = 0;
    int node_number = 0;
    int balance_no;
    complex<double> total_loss;

    induct induct_network;
    MatrixXd jacobi;
    MatrixXd delta_y;
    MatrixXd delta_x;
    double delta_e_max = LOCAL_INF; //!1 1 4 5 1 4!
    double delta_f_max = LOCAL_INF; //!1 1 4 5 1 4!
    node_arg* node;                //? nodes

    /*function*/
    void gen_delta_y();
    void gen_delta_x();
    int  gen_jacobi_part(const int start,const int end);
    void gen_jacobi();
    void gen_jacobi_mulit_thread(unsigned int thread_number);
    void gen_jacobi_mulit_thread();
    void renew_node_u();
    void set_node_u(int);
    void set_node_u();
    double get_e_delta_max();
    double get_f_delta_max();
    void gen_flow();
    complex<double>& flow(int i, int j);
    void gen_u();
    void init_network(int num_number);
    Network(int, int, int);
    ~Network();
};
void Network::renew_node_u()
{
    for (int i = 0; i < node_number - 1; ++i)
    {
        __EI__ += delta_x(2 * i + 1, 0);
        __FI__ += delta_x(2 * i, 0);
    }
    return;
}
void Network::set_node_u(int balance_node)
{
    for (int i = 0; i < node_number; ++i)
    {
        node[i].e = node[balance_node].e;
        node[i].f = node[balance_node].f;
    }

    return;
}
void Network::set_node_u()
{
    set_node_u(balance_no);

    return;
}
double Network::get_f_delta_max()
{

    double res = abs(delta_x(0, 0));
    for (int i = 0; i < matrix_length; i += 2)
    {
        if (abs(delta_x(i, 0)) > res)
            res = abs(delta_x(i, 0));
    }
    delta_f_max = res;
    return res;
}
double Network::get_e_delta_max()
{

    double res = abs(delta_x(1, 0));
    for (int i = 1; i < matrix_length; i += 2)
    {
        if (abs(delta_x(i, 0)) > res)
            res = abs(delta_x(i, 0));
    }
    delta_e_max = res;
    return res;
}
void Network::gen_delta_y()
{
    auto get_p_delta = [&](int i) -> double {
        double res = node[i].p;
        for (int j = 0; j < node_number; j++)
        {
            res -= __EI__ * (__GIJ__ * __EJ__ - __BIJ__ * __FJ__) + __FI__ * (__GIJ__ * __FJ__ + __BIJ__ * __EJ__);
        }
        return res;
    };
    auto get_q_delta = [&](int i) -> double {
        double res = node[i].q;
        for (int j = 0; j < node_number; j++)
        {
            res -= __FI__ * (__GIJ__ * __EJ__ - __BIJ__ * __FJ__) - __EI__ * (__GIJ__ * __FJ__ + __BIJ__ * __EJ__);
        }
        return res;
    };
    auto get_u_delta = [&](int i) -> double { return node[i].r * node[i].r - __EI__ * __EI__ - __FI__ * __FI__; };
    for (int i = 0; i < node_number - 1; ++i)
    {
        if (node[i].type==Node_type::pq) //?PQ nodes.
        {
            delta_y(2 * i) = get_p_delta(i);
            delta_y(2 * i + 1) = get_q_delta(i);                    
            assert(2 * i + 1 < matrix_length);//debug
        }
        if (node[i].type == Node_type::pv)//? PV nodes
        {
            delta_y(2 * i) = get_p_delta(i);
            delta_y(2 * i + 1) = get_u_delta(i);
            assert(2 * i + 1 < matrix_length);//debug
        }
    }

    return;
}
void Network::gen_delta_x()
{
     delta_x = jacobi.colPivHouseholderQr().solve(delta_y);
    //delta_x = jacobi.inverse() * delta_y;

    return;
}
void Network::gen_jacobi() {
    gen_jacobi_part(0, node_number - 1);
    return;
}
void Network::gen_jacobi_mulit_thread() {
    gen_jacobi_mulit_thread(cpu_core_number);
    return;
}
void Network::gen_jacobi_mulit_thread(unsigned int thread_number) {
    if (thread_number > cpu_core_number) {
        thread_number = cpu_core_number;
    }
    if (thread_number > node_number - 1) {
        thread_number = node_number - 1;
    }
    int per_thread = (node_number - 1) / thread_number;
    int last_thread = node_number-1-per_thread * thread_number+per_thread;
    int start = 0;
    int end = per_thread;
    vector<future<int>> fu;
    for (int i = 0;i < thread_number - 1;++i) {
        fu.push_back(std::async(&Network::gen_jacobi_part,this, start, end));
        start += per_thread;
        end += per_thread;
    }
    fu.push_back(std::async(&Network::gen_jacobi_part,this, start, node_number-1));
    for (auto& i : fu)
    {
        i.wait();
        i.get();
    }

    return;
}

int Network::gen_jacobi_part(const int start,const int end)//[start,end)
{
    ////tobedone
    auto get_a_ii = [&](int i) -> double {
        double sum = __GII__ * __EI__ - __BII__ * __FI__;
        for (int j = 0; j < node_number; ++j)
        {
            if (i != j)
            {
                sum += __GIJ__ * __EJ__ - __BIJ__ * __FJ__;
            }
        }
        return sum;
    };
    auto get_b_ii = [&](int i) -> double {
        double sum = __GII__ * __FI__ + __BII__ * __EI__;
        for (int j = 0; j < node_number; ++j)
        {
            if (i != j)
            {
                sum += __GIJ__ * __FJ__ + __BIJ__ * __EJ__;
            }
        }
        return sum;
    };


    for (int i = start; i < end; ++i)
    {
        double a = get_a_ii(i);
        double b = get_b_ii(i);
        for (int j = 0; j < node_number - 1; ++j)
        {
            auto get_h_ij = [&](int i, int j) -> double { return -__BIJ__ * __EI__ + __GIJ__ * __FI__; };
            auto get_n_ij = [&](int i, int j) -> double { return __GIJ__ * __EI__ + __BIJ__ * __FI__; };
            auto get_h_ii = [&](int i) -> double { return -__BII__ * __EI__ + __GII__ * __FI__ + b; };
            auto get_n_ii = [&](int i) -> double { return __GII__ * __EI__ + __BII__ * __FI__ + a; };
            auto get_j_ii = [&](int i) -> double { return -__GII__ * __EI__ - __BII__ * __FI__ + a; }; //!@
            auto get_l_ii = [&](int i) -> double { return -__BII__ * __EI__ + __GII__ * __FI__ - b; };
            auto get_r_ii = [&](int i) -> double { return 2 * __FI__; };
            auto get_s_ii = [&](int i) -> double { return 2 * __EI__; };

            if ( Node_type::pq==node[i].type) //?PQ nodes
            {

                if (i != j) //?non-diag
                {
                    /////!this can be rewriten as inline function to make it look batter.
                    jacobi(2 * i, 2 * j) = get_h_ij(i, j);                //? witch is H_ij
                    jacobi(2 * i, 2 * j + 1) = get_n_ij(i, j);            //? witch is N_ij
                    jacobi(2 * i + 1, 2 * j) = -jacobi(2 * i, 2 * j + 1); //? witch is J_ij
                    jacobi(2 * i + 1, 2 * j + 1) = jacobi(2 * i, 2 * j);  //? witch is L_ij
                    assert(2 * i +1 < matrix_length);//debug
                    assert(2 * j + 1 < matrix_length);//debug

                }
                else //?diag
                {
                    jacobi(2 * i, 2 * j) = get_h_ii(i);         //? witch is H_ii
                    jacobi(2 * i, 2 * j + 1) = get_n_ii(i);     //? witch is N_ii
                    jacobi(2 * i + 1, 2 * j) = get_j_ii(i);     //? witch is J_ii
                    jacobi(2 * i + 1, 2 * j + 1) = get_l_ii(i); //? witch is L_ii
                    assert(2 * i + 1 < matrix_length);//debug
                    assert(2 * j + 1 < matrix_length);//debug

                }
            }
            else if(Node_type::pv == node[i].type)//?PV nodes
            {
                //?non-diag
                if (i != j)
                {
                    jacobi(2 * i, 2 * j) = get_h_ij(i, j);     //? witch is H_ij
                    jacobi(2 * i, 2 * j + 1) = get_n_ij(i, j); //? witch is N_ij
                    jacobi(2 * i + 1, 2 * j) = 0;              //? witch is R_ij
                    jacobi(2 * i + 1, 2 * j + 1) = 0;          //? witch is S_ij
                    assert(2 * i + 1 < matrix_length);//debug
                    assert(2 * j + 1 < matrix_length);//debug
                }
                else //?diag
                {
                    jacobi(2 * i, 2 * j) = get_h_ii(i);         //? witch is H_ii
                    jacobi(2 * i, 2 * j + 1) = get_n_ii(i);     //? witch is N_ii
                    jacobi(2 * i + 1, 2 * j) = get_r_ii(i);     //? witch is R_ii
                    jacobi(2 * i + 1, 2 * j + 1) = get_s_ii(i); //? witch is S_ii
                    assert(2 * i + 1 < matrix_length);//debug
                    assert(2 * j + 1 < matrix_length);//debug
                }
            }
        }
    }
    return 0;
}
void Network::gen_flow() {
    gen_u();
    induct_network.gen_tree();

    for (int i = 0;i < node_number;++i) {
        if (node[i].type == Node_type::balance) {//balance node power 
            complex<double> sum;
            for (int j = 0;j < node_number;++j) {
                sum += induct_network[j][i] * node[j].u;
            }
            sum *= complex<double>{__EI__, __FI__};
            node[i].p = sum.real();
            node[i].q = sum.imag();

        }
        for (int j = 0;j < node_number;++j) {//line power loss
            flow(i, j) = node[i].u * ((conj(node[i].u) - conj(node[j].u) )* conj(induct_network[i][j]));
            total_loss += flow(i, j);
        }
    }
    return;
}

void Network::init_network(int num_number)
{
    return;
}
Network::Network(int pq, int pv, int total) : balance_no(total - 1), matrix_length(2 * total - 2), pv_node_number(pv), pq_node_number(pq), induct_network(total), node_number(total),
                                              jacobi(matrix_length, matrix_length), delta_x(matrix_length, 1), delta_y(matrix_length, 1)
{
    cpu_core_number = std::thread::hardware_concurrency();
    node = new node_arg[node_number];
    init_network(node_number);
}
Network::~Network()
{
    //try {
    //    delete[] node;
    //    induct_network.~induct();
    //}
    //catch (const std::exception & e) {}
}
complex<double>& Network::flow(int i,int j) {
    return (induct_network.flow)[i][j];
}
void Network::gen_u() {
    for (int i = 0;i < node_number;++i) {
    
        node[i].u = complex<double>{ __EI__,__FI__ };
 
    }

    return;
}