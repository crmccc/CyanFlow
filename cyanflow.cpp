#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include "Network.h"
#define SHOW_WIDTH 8
using namespace Eigen;
using namespace std;

void log_show_inductance(Network &);
void log_show_jacobi(Network &);
void log_show_jacobi_inverse(Network &);
void log_show_delta_x(Network &);
void log_show_delta_y(Network &);
void log_show_current_precision(Network &);
void log_show_voltage(Network &);
int node_number{0};
int pv_number{0};
int pq_number{0};
int line_number{0};
double precision;

int main()
{
    cout << setw(SHOW_WIDTH); //?debug
    // FILE *input_file = fopen("input.txt", "r");
    // fscanf(input_file, "%d %d %d %d %f", &node_number, &pv_number, &pq_number, &line_number, &precision);
    fstream fin("input.txt");
    fin >> node_number >> pv_number >> pq_number >> line_number >> precision;
    int temp;

    Network network(pv_number, pq_number, node_number);

    for (int i = 0; i < line_number; ++i)
    {
        int from, to;
        double real, imag;
        // fscanf(input_file, "%d %d %d %f %f", &temp, &from, &to, &real, &imag);
        fin >> temp >> from >> to >> real >> imag;
        complex<double> tempcd{real, imag};
        network.induct_network.add_line(tempcd, to - 1, from - 1);
    }
    for (int i = 0; i < node_number; ++i)
    {
        char type_c;
        double real, imag;
        // fscanf(input_file, "%d %c %f %f", &temp, &type_c, &real, &imag);
        fin >> temp >> type_c >> real >> imag;
        temp -= 1;
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
            network.node[temp].e = real * cos(imag);
            network.node[temp].f = real * sin(imag);
            network.node[temp].type = Network::Node_type::balance;
        }
    }
    //start..
    network.set_node_u();
    int iteration = 1;
    log_show_inductance(network); //?debug
    while (1)
    {
        network.gen_delta_y();
        network.gen_jacobi();
        network.gen_delta_x();
        network.get_f_delta_max();
        network.get_e_delta_max();
        cout << "interation:" << iteration << '\n'; //?debug
        log_show_delta_x(network);                  //?debug
        log_show_delta_y(network);                  //?debug
        log_show_jacobi(network);                   //?debug
        log_show_jacobi_inverse(network);           //?debug
        log_show_current_precision(network);        //?debug
        log_show_voltage(network);
        if (network.delta_e_max < precision && network.delta_f_max < precision)
        {
            break;
        }
        network.renew_node_u();
        ++iteration;
    }

    return 0;
}

void log_show_inductance(Network &net)
{
    cout << setw(SHOW_WIDTH); //?debug

    cout << "inductance:\n";
    for (int i = 0; i < node_number; ++i)
    {
        for (int j = 0; j < node_number; ++j)
        {
            cout << net.induct_network[i][j] << ' ';
        }
        cout << '\n';
    }

    return;
}
void log_show_jacobi(Network &net)
{
    cout << setw(SHOW_WIDTH); //?debug

    cout << "jacobi:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        for (int j = 0; j < net.matrix_length; ++j)
        {
            cout << net.jacobi(i, j) << ' ';
        }
        cout << '\n';
    }
    return;
}
void log_show_delta_x(Network &net)
{
    cout << setw(SHOW_WIDTH); //?debug

    cout << "delta_x:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        cout << net.delta_x(i, 0) << '\n';
    }
    return;
}
void log_show_delta_y(Network &net)
{
    cout << setw(SHOW_WIDTH); //?debug

    cout << "delta_y:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        cout << net.delta_y(i, 0) << '\n';
    }
    return;
}
void log_show_current_precision(Network &net)
{
    cout << setw(SHOW_WIDTH); //?debug

    cout << "delta e:" << net.delta_e_max << '\n';
    cout << "delta f:" << net.delta_f_max << '\n';
    return;
}
void log_show_voltage(Network &net)
{
    cout << setw(SHOW_WIDTH) << "voltage:\n";
    for (int i = 0; i < net.node_number; ++i)
    {
        cout << "e: " << net.node[i].e << "f: " << net.node[i].f << '\n';
    }

    return;
}
void log_show_jacobi_inverse(Network &net)
{
    auto& inv = net.jacobi.inverse();

    cout << setw(SHOW_WIDTH); //?debug

    cout << "jacobi inverse:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        for (int j = 0; j < net.matrix_length; ++j)
        {
            cout << inv(i, j) << ' ';
        }
        cout << '\n';
    }
}