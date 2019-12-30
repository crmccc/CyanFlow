#include<complex>
#include<map>
#include<list>
#include<vector>
//some compromise
using namespace std;
/*
*class Network.
*storage network conductance
*/
class Network
{
private:
    /* data */
public:
    const int MAX_NODE_NUMBER=100; 
    complex<double> data[MAX_NODE_NUMBER][MAX_NODE_NUMBER];
    int line_number=0;
    int node_number=0;

    int add_node(vector<map<int,double>>);
    
    Network(/* args */);
    ~Network();
};


Network::Network(/* args */)
{
}

Network::~Network()
{
}

int Network::add_node(vector<map<int,double>> conduct){
//return -1 for data out of range
//reutrn 0 for normal exit
    if (this.node_number>=this.MAX_NODE_NUMBER){
        return -1;
    }
    ++this.node_number;
    int i=0;
    for(auto it=conduct.begin();it!=conduct.end()-1;++it){
        data[i][node_number]=data[node_number][i]+=
        ++i;
    }
}