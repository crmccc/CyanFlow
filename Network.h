#include<complex>
#include"Eigen/Dense"
#ifndef MAX_NODE_NUMBER
#define MAX_NODE_NUMBER 100
#endif
//some compromise
using Eigen::MatrixXd;
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
        double e;
        double f;
        /* data */
    };
    int pv_node_number=0;
    int pq_node_number=0;
    int num_number=0;
    
    node_u u [MAX_NODE_NUMBER];
    induct induct_network;
    Matrix<double,MAX_NODE_NUMBER,MAX_NODE_NUMBER> jacobi;

    int init_network(int node_number);
    
    Network(/* args */);
    ~Network();
};
int Network::init_network(int node_number){
    for(int i=0;i<node_number;++i){
        u[node_number].e=1.0;
        u[node_number].f=0.0;
        
    }
}
Network::Network(/* args */)
{
}

Network::~Network()
{
}
