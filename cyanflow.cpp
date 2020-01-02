#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <cmath>
#include "Network.h"

using namespace Eigen;
using namespace std;

int node_number{0};
int pv_number{0};
int pq_number{0};
int line_number{0};

double precision;
int main()
{
    FILE *input_file = fopen("input.txt", "r");
    fscanf(input_file, "%d %d %d %d %f", &node_number, &pv_number, &pq_number, &line_number, &precision);
    int temp;

    Network network(pv_number, pq_number, node_number);

    for (int i = 0; i < line_number; ++i)
    {
        int from, to;
        double real, imag;
        fscanf(input_file, "%d %d %d %f %f", &temp, &from, &to, &real, &imag);
        network.induct_network.add_line((complex<double>{real, imag}), to, from);
    }
    for (int i = 0; i < node_number; ++i)
    {
        char type_c;
        double real, imag;
        fscanf(input_file, "%d %c %f %f", &temp, &type_c, &real, &imag);
        if ('Q' == type_c)
        {
            network.node[temp].p = real;
            network.node[temp].q = imag;
            network.node[temp].node_type = Network::Node_type::pq;
        }
        else if ('V' == type_c)
        {
            network.node[temp].p=real;
            network.node[temp].r= imag;
            network.node[temp].node_type = Network::Node_type::pv;
        }
        else
        {
            network.node[temp].r=real;
            network.node[temp].angle=imag;
            network.node[temp].node_type = Network::Node_type::balance;
        }
    }

    return 0;
}