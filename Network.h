#include "Eigen/Eigen"
#include "induct.h"

#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif

#ifndef MAX_NODE_NUMBER
#define MAX_NODE_NUMBER 100
#endif
#define GII induct_network[i][i].imag
#define GIJ induct_network[i][j].imag
#define BII induct_network[i][i].real
#define BIJ induct_network[i][j].real
#define EI u[i].e
#define EJ u[j].e
#define FI u[i].f
#define FJ u[j].f
#define EXIST(i, j) induct_network.book[i][j]

using Eigen::MatrixXd;
using std::complex;
/*
*class Network.
*storage network conductance
*/
//!this flow calculate library obey "first PQ then PV" arrange.
class Network
{
private:
    /* data */
public:
    struct node_u
    {
        double e = 1.0;
        double f = 0.0;

        /* data */
    };
    int pv_node_number = 0;
    int pq_node_number = 0;
    int node_number = 0;

    induct induct_network;
    MatrixXd jacobi;
    MatrixXd delta_y;
    MatrixXd delta_x;
    node_u u[]; //the potiential of nodes

    void gen_delta_y();
    void gen_delta_x();
    void gen_jacobi();

    void init_network(int num_number);
    Network(int, int, int);
    ~Network();
};

void Network::gen_delta_y()
{
}
void Network::gen_delta_x()
{
    delta_x = delta_y.inverse() * delta_y;
}
void Network::gen_jacobi()
{
    //!tobedone
    auto get_a_ii = [&](int i) {
        double sum = GII * EI - BII * FI;
        for (int j = 0; j < node_number; ++j)
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
        for (int j = 0; j < node_number; ++j)
        {
            if (EXIST(i, j))
            {
                sum += GIJ * FJ + BIJ * EJ;
            }
        }
        return sum;
    };

    int k = pq_node_number;

    for (int i = 0; i < node_number; i += 2)
    {
        double a = get_a_ii(i);
        double b = get_b_ii(i);
        for (int j = 0; j < node_number; j += 2)
        {
            auto get_h_ij = [&](int i, int j) { return EXIST(i, j) ? -BIJ * EI + GIJ * FI : 0; };
            auto get_j_ij = [&](int i, int j) { return EXIST(i, j) ? GIJ * EI + BIJ * FI : 0; };
            auto get_h_ii = [&](int i) { return -BII * EI + GII * FI + b; };
            auto get_n_ii = [&](int i) { return GII * EI + BII * FI + a; };
            auto get_j_ii = [&](int i) { return -GII * EI - BII * FI + a; };
            auto get_l_ii = [&](int i) { return -BII * EI + GII * FI - b; };
            auto get_r_ii = [&](int i) { return 2 * FI; };
            auto get_s_ii = [&](int i) { return 2 * EI; };
            
            if (k--)//?PQ nodes
            {
                
                if (i != j)//?non-diag
                { 
                    /////!this can be rewriten as inline function to make it look batter.
                    jacobi(i, j) = get_h_ij(i, j);        //? witch is H_ij
                    jacobi(i + 1, j) = get_j_ij(i, j);    //? witch is N_ij
                    jacobi(i, j + 1) = -jacobi(i + 1, j); //? witch is J_ij
                    jacobi(i + 1, j + 1) = jacobi(i, j);  //? witch is L_ij
                }
                else//?diag
                {                                       
                    jacobi(i, j) = get_h_ii(i);         //? witch is H_ii
                    jacobi(i + 1, j) = get_n_ii(i);     //? witch is N_ii
                    jacobi(i, j + 1) = get_j_ii(i);     //? witch is J_ii
                    jacobi(i + 1, j + 1) = get_l_ii(i); //? witch is L_ii
                }
            }
            else//?PV nodes
            {
                //?non-diag
                if (i != j)
                {
                    jacobi(i, j) = get_h_ij(i, j);     //? witch is H_ij
                    jacobi(i + 1, j) = get_j_ij(i, j); //? witch is N_ij
                    jacobi(i, j + 1) = 0;              //? witch is R_ij
                    jacobi(i + 1, j + 1) = 0;          //? witch is S_ij
                }
                else//?diag
                {
                    jacobi(i, j) = get_h_ii(i);         //? witch is H_ii
                    jacobi(i + 1, j) = get_n_ii(i);     //? witch is N_ii
                    jacobi(i, j + 1) = get_r_ii(i);     //? witch is R_ii
                    jacobi(i + 1, j + 1) = get_s_ii(i); //? witch is S_ii
 
                }
            }
        }
    }
}

void Network::init_network(int num_number)
{
}
Network::Network(int pv, int pq, int total) : node_number(total), pv_node_number(pv), pq_node_number(pq), induct_network(total),
                                              jacobi((2 * total - 2), (2 * total - 2)), delta_x(1, 2 * total - 2), delta_y(1, 2 * total - 2)
{
    init_network(node_number);
}
Network::~Network()
{
}
