#include "Eigen/Dense"
#include "induct.h"

#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif

#ifndef MAX_NODE_NUMBER
#define MAX_NODE_NUMBER 100
#endif
//some compromise
using Eigen::MatrixXd;
using std::complex;
/*
*class Network.
*storage network conductance
*/
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
    for(int i=0;i<node_number;i+=2){
        for (int j=0;j<node_number;j+=2){
            
        }
    }
}

void Network::init_network(int num_number)
{
}
Network::Network(int pv, int pq, int total) : node_number(total), pv_node_number(pv), pq_node_number(pq),
                                              jacobi((2 * total - 2), (2 * total - 2)), delta_x(1, 2 * total - 2), delta_y(1, 2 * total - 2)
{
    init_network(node_number);
}
Network::~Network()
{
}
