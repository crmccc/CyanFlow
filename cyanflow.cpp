/*
*Author:crmccc@126.com
*Create time:2020.1
*Github repo:https://github.com/crmccc/CyanFlow
*
*This program is licensed under the GNU General Public License v3.0
*
*I think it will handle 100+ nodes with ease.but I had never tryed.
*The parser.py (py3.6+) will convert any input file to my program input format.
*
*/
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "Network.h"

constexpr auto SHOW_WIDTH = 5;
constexpr auto MAX_ITERATION = 1000;

using namespace Eigen;
using namespace std;

const char file_name[] = "input.txt";
const char output_file[] = "output.txt";
ofstream fout("output.txt");

void log_show_inductance(Network &);        //print inductance network.
void log_show_jacobi(Network &);            //print jacobi matrix.
void log_show_jacobi_inverse(Network &);    //print jacobi matrix inverses.
void log_show_delta_x(Network &);           //print delta x.
void log_show_delta_y(Network &);           //print delta x.
void log_show_current_precision(Network &); //print current precision.
void log_show_voltage(Network &);           //print current voltage.
void log_show_node_arg(Network &);          //print node data.
void log_show_final(Network &);             //print final stage data.

class highground
{ //just for out put format
    complex<double> &temp;

public:
    friend ostream &operator<<(ostream &stream, highground &cd)
    { // overload << operator
        stream << cd.temp.real() << (cd.temp.imag() >= 0 ? "+j" : "-j") << fabs(cd.temp.imag());
        return stream;
    };
    highground(complex<double> &a) : temp(a){};
};

int main()
{

    int node_number{0}; //As its name suggested.
    int pv_number{0};   //As its name suggested.
    int pq_number{0};   //As its name suggested.
    int line_number{0}; //As its name suggested.
    double precision;   //As its name suggested.

    //string file_name = to_string(current_file_number) + string(".txt");
    fout << setw(SHOW_WIDTH); //?debug
    // FILE *input_file = fopen("input.txt", "r");
    // fscanf(input_file, "%d %d %d %d %f", &node_number, &pv_number, &pq_number, &line_number, &precision);
    ifstream fin(file_name);
    fin >> node_number >> pq_number >> pv_number >> line_number >> precision;

    Network network(pq_number, pv_number, node_number);

    for (int i = 0; i < line_number; ++i)
    {
        int from, to, temp;
        double real, imag;
        // fscanf(input_file, "%d %d %d %f %f", &temp, &from, &to, &real, &imag);
        fin >> temp >> from >> to >> real >> imag;
        complex<double> tempcd{real, imag};
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
        if ('Q' == type_c) //Deal with PQ nodes.
        {
            network.node[temp].p = real;
            network.node[temp].q = imag;
            network.node[temp].type = Network::Node_type::pq;
        }
        else if ('V' == type_c) //Deal with PV nodes.
        {
            network.node[temp].p = real;
            network.node[temp].r = imag;
            network.node[temp].type = Network::Node_type::pv;
        }
        else //Deal with balance node.
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
    network.set_node_u(); //set all node voltage to balance node's.
    int iteration = 1;
    log_show_inductance(network); //?debug
    log_show_node_arg(network);   //?debug
    while (iteration < MAX_ITERATION)
    {
        network.gen_jacobi();                                                   //As its name suggested.
        network.gen_delta_y();                                                  //As its name suggested.
        network.gen_delta_x();                                                  //As its name suggested.
        network.get_f_delta_max();                                              //As its name suggested.
        network.get_e_delta_max();                                              //As its name suggested.
        fout << "interation:" << iteration << '\n';                             //?debug
        log_show_jacobi(network);                                               //?debug
        log_show_jacobi_inverse(network);                                       //?debug
        log_show_delta_y(network);                                              //?debug
        log_show_delta_x(network);                                              //?debug
        log_show_current_precision(network);                                    //?debug
        log_show_voltage(network);                                              //?debug
        if (network.delta_e_max < precision && network.delta_f_max < precision) //Done?
        {
            break;
        }
        network.renew_node_u(); //As its name suggested.
        fout << "*******************************\n";
        ++iteration;
    }
    network.gen_flow();
    log_show_final(network);
    if (iteration >= MAX_ITERATION) //If not converage.
    {
        fout << "DOESN'T CONVERAGE!";
        //todo it doesn't converage
    }
    return 0;
}

void log_show_inductance(Network &net)
{
    fout << setw(SHOW_WIDTH); //?debug

    fout << "inductance:\n";
    for (int i = 0; i < net.node_number; ++i)
    {
        for (int j = 0; j < net.node_number; ++j)
        {
            fout << setw(15) << setprecision(5) << highground(net.induct_network[i][j]) << '\t'; //Its complex<double>,take the highground.
        }
        fout << '\n';
    }

    return;
}
void log_show_jacobi(Network &net)
{
    fout << setw(SHOW_WIDTH); //?debug

    fout << "jacobi:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        for (int j = 0; j < net.matrix_length; ++j)
        {
            fout << setw(15) << setprecision(5) << net.jacobi(i, j) << ' ';
        }
        fout << '\n';
    }
    return;
}
void log_show_delta_x(Network &net)
{
    fout << setw(SHOW_WIDTH); //?debug

    fout << "delta_x:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        fout << net.delta_x(i, 0) << '\n';
    }
    return;
}
void log_show_delta_y(Network &net)
{
    fout << setw(SHOW_WIDTH); //?debug

    fout << "delta_y:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        fout << net.delta_y(i, 0) << '\n';
    }
    return;
}
void log_show_current_precision(Network &net)
{
    fout << setw(SHOW_WIDTH); //?debug

    fout << "max delta e:" << net.delta_e_max << '\n';
    fout << "max delta f:" << net.delta_f_max << '\n';
    return;
}
void log_show_voltage(Network &net)
{
    fout << setw(SHOW_WIDTH) << "voltage:\n";
    for (int i = 0; i < net.node_number; ++i)
    {
        fout << "node " << i + 1 << " : e: " << net.node[i].e << " f: " << net.node[i].f << '\n';
    }

    return;
}
void log_show_jacobi_inverse(Network &net)
{
    auto &inv = net.jacobi.inverse();

    fout << setw(SHOW_WIDTH); //?debug

    fout << "jacobi inverse:\n";
    for (int i = 0; i < net.matrix_length; ++i)
    {
        for (int j = 0; j < net.matrix_length; ++j)
        {
            fout << setw(15) << setprecision(5) << inv(i, j) << ' ';
        }
        fout << '\n';
    }
}
void log_show_node_arg(Network &net)
{
    fout << "node_number: " << net.node_number;
    fout << "\nPV_number: " << net.pv_node_number;
    fout << "\nPQ_number: " << net.pq_node_number << '\n';

    int i = 0;
    for (int i = 0; i < net.node_number; ++i)
    {
        if (Network::Node_type::pq == net.node[i].type)
        {
            fout << '[' << i + 1 << "]: PQ node : p = " << net.node[i].p << " q = " << net.node[i].q << '\n';
        }
        if (Network::Node_type::pv == net.node[i].type)
        {
            fout << '[' << i + 1 << "]: PV node : p = " << net.node[i].p << " r = " << net.node[i].r << '\n';
        }
        if (Network::Node_type::balance == net.node[i].type)
        {
            fout << '[' << i + 1 << "]: balance node : r = " << net.node[i].r << " angle = " << net.node[i].angle << '\n';
        }
        fout << '[' << i + 1 << "]: e = " << net.node[i].e << " f = " << net.node[i].f << '\n';
    }
    return;
}

void log_show_final(Network &net)
{
    fout << "\nFinal result:\n";
    fout << "node voltage:\n";
    for (int i = 0; i < net.node_number; ++i)
    {
        fout << i + 1 << " node: " << highground(net.node[i].u) << '\n'; //Its complex<double>,take the highground.
    }
    fout << "node power:\n";
    for (int i = 0; i < net.node_number; ++i)
    {
        fout << i + 1 << " node: " << highground(complex<double>{net.node[i].p, net.node[i].q}) << '\n'; //Its complex<double>,take the highground.
    }
    fout << "line power flow:\n";
    int iter = net.induct_network.start;
    fout << iter + 1;
    for (int j = 0; j < net.node_number; ++j)
    {

        for (int i = 0; i < net.node_number; ++i)
        {
            if (net.induct_network.book[iter][i])
            {
                fout << " ---" << highground(net.flow(iter, i)) << "---> " << i + 1;
                net.induct_network.book[iter][i] = net.induct_network.book[i][iter] = 0;
                iter = i;
                break;
            }
        }
    }
    fout << '\n';
    fout << "total power loss: " << highground(net.total_loss); //Its complex<double>,take the highground.
    return;
}
//just for output complex<double> in ideal format.
