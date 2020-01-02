#include "Eigen/Eigen"
#include "induct.h"

#include <cmath>
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#include <complex>
//ya,I know that looks stupid, but it SHUANG SI LE.

#define GII induct_network[i][i].imag()
#define GIJ induct_network[i][j].imag()
#define BII induct_network[i][i].real()
#define BIJ induct_network[i][j].real()
#define EI node[i].e
#define EJ node[j].e
#define FI node[i].f
#define FJ node[j].f
#define EXIST(i, j) induct_network.book[i][j]
#define LOCAL_INF 114514

using namespace Eigen;
using std::complex;
/*
*class Network.
*OOP
*/
//!this flow calculate library obey "first PQ then PV" arrange.
//!balance node must be the last one
class Network
{
private:
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

    induct induct_network;
    MatrixXd jacobi;
    MatrixXd delta_y;
    MatrixXd delta_x;
    double delta_e_max = LOCAL_INF; //!1 1 4 5 1 4!
    double delta_f_max = LOCAL_INF; //!1 1 4 5 1 4!
    node_arg node[];                //? nodes

    /*function*/
    void gen_delta_y();
    void gen_delta_x();
    void gen_jacobi();
    void renew_node_u();
    void set_node_u(int);
    void set_node_u();
    double get_e_delta_max();
    double get_f_delta_max();

    void init_network(int num_number);
    Network(int, int, int);
    ~Network();
};
void Network::renew_node_u()
{
    for (int i = 0; i < node_number; ++i)
    {
        EI = delta_x(2 * i, 0);
        FI = delta_x(2 * i + 1, 0);
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
}
void Network::set_node_u()
{
    set_node_u(balance_no);

    return;
}
double Network::get_f_delta_max()
{
    /////todo work in prograss
    double res = delta_x(0, 0);
    for (int i = 0; i < matrix_length; i += 2)
    {
        if (delta_x(i, 0) > res)
            res = delta_x(i, 0);
    }
    delta_f_max = res;
    return res;
}
double Network::get_e_delta_max()
{
    /////todo work in prograss
    double res = delta_x(1, 0);
    for (int i = 1; i < matrix_length; i += 2)
    {
        if (delta_x(i, 0) > res)
            res = delta_x(i, 0);
    }
    delta_e_max = res;
    return res;
}
void Network::gen_delta_y()
{
    auto get_p_delta = [&](int i) {
        double res = node[i].p;
        for (int j = 0; j < node_number; j++)
        {
            res -= EI * (GIJ * EJ - BIJ * FJ) + FI * (GIJ * FJ + BIJ * EJ);
        }
        return res;
    };
    auto get_q_delta = [&](int i) {
        double res = node[i].q;
        for (int j = 0; j < node_number; j++)
        {
            res -= FI * (GIJ * EJ - BIJ * FJ) - EI * (GIJ * FJ + BIJ * EJ);
        }
        return res;
    };
    //!issue s in get u delta
    auto get_u_delta = [&](int i) -> double { return node[i].r * node[i].r - EI * EI - FI * FI; };
    for (int i = 0; i < matrix_length; i += 2)
    {
        if (i < pq_node_number) //?PQ nodes.
        {
            delta_y(i) = get_p_delta(i);
            delta_y(i + 1) = get_q_delta(i);
        }
        else //? PV nodes
        {
            delta_y(i) = get_p_delta(i);
            delta_y(i + 1) = get_u_delta(i);
        }
    }

    return;
}
void Network::gen_delta_x()
{
    delta_x = jacobi.inverse() * delta_y;

    return;
}
void Network::gen_jacobi()
{
    ////tobedone
    auto get_a_ii = [&](int i) {
        double sum = GII * EI - BII * FI;
        for (int j = 0; j < matrix_length; ++j)
        {
            if (EXIST(i, j))
            {
                sum += GIJ * EJ - BIJ * FJ;
            }
        }
        return sum;
    };
    auto get_b_ii = [&](int i) {
        double sum = GII * FI + BII * EI;
        for (int j = 0; j < matrix_length; ++j)
        {
            if (EXIST(i, j))
            {
                sum += GIJ * FJ + BIJ * EJ;
            }
        }
        return sum;
    };

    int k = pq_node_number;

    for (int i = 0; i < matrix_length; i += 2)
    {
        double a = get_a_ii(i);
        double b = get_b_ii(i);
        for (int j = 0; j < matrix_length; j += 2)
        {
            auto get_h_ij = [&](int i, int j) -> double { return EXIST(i, j) ? -BIJ * EI + GIJ * FI : 0; };
            auto get_j_ij = [&](int i, int j) -> double { return EXIST(i, j) ? GIJ * EI + BIJ * FI : 0; };
            auto get_h_ii = [&](int i) -> double { return -BII * EI + GII * FI + b; };
            auto get_n_ii = [&](int i) -> double { return GII * EI + BII * FI + a; };
            auto get_j_ii = [&](int i) -> double { return -GII * EI - BII * FI + a; };
            auto get_l_ii = [&](int i) -> double { return -BII * EI + GII * FI - b; };
            auto get_r_ii = [&](int i) -> double { return 2 * FI; };
            auto get_s_ii = [&](int i) -> double { return 2 * EI; };

            if (k--) //?PQ nodes
            {

                if (i != j) //?non-diag
                {
                    /////!this can be rewriten as inline function to make it look batter.
                    jacobi(i, j) = get_h_ij(i, j);        //? witch is H_ij
                    jacobi(i + 1, j) = get_j_ij(i, j);    //? witch is N_ij
                    jacobi(i, j + 1) = -jacobi(i + 1, j); //? witch is J_ij
                    jacobi(i + 1, j + 1) = jacobi(i, j);  //? witch is L_ij
                }
                else //?diag
                {
                    jacobi(i, j) = get_h_ii(i);         //? witch is H_ii
                    jacobi(i + 1, j) = get_n_ii(i);     //? witch is N_ii
                    jacobi(i, j + 1) = get_j_ii(i);     //? witch is J_ii
                    jacobi(i + 1, j + 1) = get_l_ii(i); //? witch is L_ii
                }
            }
            else //?PV nodes
            {
                //?non-diag
                if (i != j)
                {
                    jacobi(i, j) = get_h_ij(i, j);     //? witch is H_ij
                    jacobi(i + 1, j) = get_j_ij(i, j); //? witch is N_ij
                    jacobi(i, j + 1) = 0;              //? witch is R_ij
                    jacobi(i + 1, j + 1) = 0;          //? witch is S_ij
                }
                else //?diag
                {
                    jacobi(i, j) = get_h_ii(i);         //? witch is H_ii
                    jacobi(i + 1, j) = get_n_ii(i);     //? witch is N_ii
                    jacobi(i, j + 1) = get_r_ii(i);     //? witch is R_ii
                    jacobi(i + 1, j + 1) = get_s_ii(i); //? witch is S_ii
                }
            }
        }
    }
    return;
}

void Network::init_network(int num_number)
{
    return;
}
Network::Network(int pv, int pq, int total) : balance_no(total-1), matrix_length(2 * total - 2), pv_node_number(pv), pq_node_number(pq), induct_network(total), node_number(total),
                                              jacobi(matrix_length, matrix_length), delta_x(matrix_length, 1), delta_y(2 * total - 2, 1)
{
    init_network(matrix_length);
}
Network::~Network()
{
    induct_network.~induct();
}
