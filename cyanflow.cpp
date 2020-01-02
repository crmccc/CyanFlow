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
    // FILE *input_file = fopen("input.txt", "r");
    // fscanf(input_file, "%d %d %d %d %f", &node_number, &pv_number, &pq_number, &line_number, &precision);
    fstream fin("input.txt");
    fin >> node_number >> pv_number >> pq_number >> line_number >> precision;
    int temp;

    extern Network network(pv_number, pq_number, node_number);

    for (int i = 0; i < line_number; ++i)
    {
        int from, to;
        double real, imag;
        // fscanf(input_file, "%d %d %d %f %f", &temp, &from, &to, &real, &imag);
        fin >> temp >> from >> to >> real >> imag;
        complex<double> tempcd{real, imag};
        network.induct_network.add_line(tempcd, to, from);
    }
    for (int i = 0; i < node_number; ++i)
    {
        char type_c;
        double real, imag;
        // fscanf(input_file, "%d %c %f %f", &temp, &type_c, &real, &imag);
        fin >> temp >> type_c >> real >> imag;
        if ('Q' == type_c)
        {
            network.node[temp].p = real;
            network.node[temp].q = imag;
            network.node[temp].type = Network::Node_type::pq;
        }
        else if ('V' == type_c)
        {
            network.node[temp].p = real;
            network.node[temp].r = imag;
            network.node[temp].type = Network::Node_type::pv;
        }
        else
        {
            network.node[temp].r = real;
            network.node[temp].angle = imag;
            network.node[temp].type = Network::Node_type::balance;
        }
    }
    //start..
    while(1){
        
    }

    return 0;
}

void log_show_inductance(Network &net)
{
    for (int i = 0; i < node_number; ++i)
    {
        for (int j = 0; j < node_number; ++j)
        {
            cout << "inductance:\n"
                 << net.induct_network[i][j] << ' ';
        }
        cout << '\n';
    }

    return;
}
void log_show_jacobi(Network &net)
{
    for (int i = 0; i < net.matrix_length; ++i)
    {
        for (int j = 0; j < net.matrix_length; ++j)
        {
            cout << "jacobi:\n"
                 << net.jacobi(i, j) << ' ';
        }
        cout << '\n';
    }
    return;
}
void log_show_delta_x(Network &net)
{
    for (int i = 0; i < net.matrix_length; ++i)
    {
        cout << "jacobi:\n"
             << net.delta_x(i, 0) << '\n';
    }
    return;
}
void log_show_delta_y(Network &net)
{
    for (int i = 0; i < net.matrix_length; ++i)
    {
        cout << "jacobi:\n"
             << net.delta_y(i, 0) << '\n';
    }
    return;
}