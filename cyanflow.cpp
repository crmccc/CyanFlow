#include <fstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include "Network.h"
#include <string>
constexpr auto SHOW_WIDTH = 5;
constexpr auto MAX_ITERATION = 1000;
using namespace Eigen;
using namespace std;

void log_show_inductance(Network &);
void log_show_jacobi(Network &);
void log_show_jacobi_inverse(Network &);
void log_show_delta_x(Network &);
void log_show_delta_y(Network &);
void log_show_current_precision(Network &);
void log_show_voltage(Network &);
void log_show_node_arg(Network &);
void log_show_final(Network&);

int node_number{0};
int pv_number{0};
int pq_number{0};
int line_number{0};
double precision;

int main()
{
    
    //const char file_name[]="example.txt";
    for (int current_file_number = 1;current_file_number < 90;++current_file_number) {
        string file_name = to_string(current_file_number) + string(".txt");
        cout << setw(SHOW_WIDTH); //?debug
        // FILE *input_file = fopen("input.txt", "r");
        // fscanf(input_file, "%d %d %d %d %f", &node_number, &pv_number, &pq_number, &line_number, &precision);
        fstream fin(file_name);
        fin >> node_number >> pq_number >> pv_number >> line_number >> precision;

        Network network(pq_number, pv_number, node_number);

        for (int i = 0; i < line_number; ++i)
        {
            int from, to, temp;
            double real, imag;
            // fscanf(input_file, "%d %d %d %f %f", &temp, &from, &to, &real, &imag);
            fin >> temp >> from >> to >> real >> imag;
            complex<double> tempcd{ real, imag };
            network.induct_network.add_line(tempcd, to - 1, from - 1);
        }
        for (int i = 0; i < node_number; ++i)
        {
            char type_c;
            int temp;
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
        fin.close();
        //start..
        network.set_node_u();
        int iteration = 1;
        //log_show_inductance(network); //?debug
        //log_show_node_arg(network);   //?debug
        cout << "start iteration:\n";
        while (iteration < MAX_ITERATION)
        {
            network.gen_jacobi();
            network.gen_delta_y();
            network.gen_delta_x();
            network.get_f_delta_max();
            network.get_e_delta_max();
            //cout << "interation:" << iteration << '\n'; //?debug
            //log_show_jacobi(network);                   //?debug
            //log_show_jacobi_inverse(network);           //?debug
            //log_show_delta_y(network);                  //?debug
            //log_show_delta_x(network);                  //?debug
            //log_show_current_precision(network);        //?debug
            //log_show_voltage(network);                  //?debug
            if (network.delta_e_max < precision && network.delta_f_max < precision)
            {
                break;
            }
            network.renew_node_u();
            //cout << "*******************************\n";
            ++iteration;
        }
        network.gen_flow();
        log_show_final(network);

        if (iteration >= MAX_ITERATION) {
            cout<<"data number : "<<current_file_number<<" doesn't converge\n";
            //todo it doesnt converage
        }
        else
        {
            cout << "data number : " << current_file_number << "converaged at " << iteration << '\n';
        }
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
            cout << setw(10)<< setprecision(3)<<net.induct_network[i][j] << ' ';
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
            cout << setw(SHOW_WIDTH)<<setprecision(3)<<net.jacobi(i, j) << ' ';
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
    auto &inv = net.jacobi.inverse();

    cout << setw(SHOW_WIDTH); //?debug

    cout << "jacobi inverse:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        for (int j = 0; j < net.matrix_length; ++j)
        {
            cout<< setw(SHOW_WIDTH)<<setprecision(SHOW_WIDTH)<<inv(i, j) << ' ';
        }
        cout << '\n';
    }
}
void log_show_node_arg(Network &net)
{
    cout<<"node_number: "<<net.node_number;
    cout<<"\nPV_number: "<<net.pv_node_number;
    cout<<"\nPQ_number: "<<net.pq_node_number<<'\n';

    int i = 0;
    for (int i=0; i < node_number; ++i)
    {
        if (net.node[i].type == Network::Node_type::pq) {
            cout << '[' << i << "]: p = " << net.node[i].p << " q = " << net.node[i].q << '\n';
            cout << '[' << i << "]: e = " << net.node[i].e << " f = " << net.node[i].f << '\n';
        }
        if (net.node[i].type == Network::Node_type::pv) {
            cout << '[' << i << "]: p = " << net.node[i].p << " r = " << net.node[i].r << '\n';
            cout << '[' << i << "]: e = " << net.node[i].e << " f = " << net.node[i].f << '\n';

        }  
        if (net.node[i].type == Network::Node_type::balance) {
            cout << '[' << i << "]: r = " << net.node[i].r << " angle = " << net.node[i].r << '\n';
        }
    }
    return;
}

void log_show_final(Network& net) {
    cout << "node_voltage:\n";
    for (int i = 0;i < net.node_number;++i) {
        cout <<i+1<<"node: "<< net.node[i].u<<'\n';
    }
    cout << "node_power:\n";
    for (int i = 0;i < net.node_number;++i) {
        cout << i + 1 << "node: " << complex<double>{net.node[i].p,net.node[i].q} << '\n';
    }
    cout << "line loss:\n";
    int iter = net.induct_network.start;
    cout << iter+1;
    while (iter != net.induct_network.end) {

        for (int i = 0;i < node_number;++i) {
            if (net.induct_network.book[iter][i]) {
                cout <<"----" << net.flow(i, iter) << "---" << i+1;
                net.induct_network.book[iter][i] = net.induct_network.book[i][iter]= 0;
                iter = i;
                break;
            }
        }
    }

    return;

}