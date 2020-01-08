/*
*Author:crmccc@126.com
*Create time:2020.1
*Github repo:https://github.com/crmccc/CyanFlow
*
*This program is licensed under the GNU General Public License v3.0
*
*This head file defined and implemented an inductance network
*/

#ifndef _GLIBCXX_MAP
#include <map>
#endif
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#ifndef _INDUCT_CYAN_
#define _INDUCT_CYAN_
#endif

using std::map;
//#include<list>
using std::complex;
//using std::list;

/*
*!THE FIRST NODE IS MARKED AS 0,NOT 1!
*This class genurate inductance network for power network power flow calculate.
*This class is specifically designed for the specific lession.
*Example: 
*   induct ind(node_number);
*   ind.add_line({1,4},1,1);
*/
class induct
{
private:
    complex<double> zero_complex{0, 0};    //To avoid the cost of copy complex,make a zero.
    map<int, complex<double>> *inductance; //The data of the inductance network.
public:
    bool **book; //whether i,j has a line or not?

    // bool book[100][100]{0};
    int node_number; //as its name suggested
    //vector<map<int, complex<double>>> inductance;
    //complex<double>** inductance;
    complex<double> **flow; //the flow data
    int start;
    //The start point of the network ,mainly used in generating hte tree of the network.
    //complex<double> inductance[100][100];

    int add_line(complex<double> &, int, int);
    void add_line_with_transformer(complex<double> &ind, int f, int t, double k);
    map<int, complex<double>> &operator[](int i);
    complex<double> &operator()(int i, int j);
    void gen_tree();

    induct(int);
    ~induct();
};
//Overload the () operator for a batter user experience.
//Now you can get the inductance of two point like that:
//  induct inductance(node_number);
//  cout<<inductance(1,5);
complex<double> &induct::operator()(int i, int j)
{
    if (book[i][j] == false && i != j)
    {
        return zero_complex;
    }
    else
    {
        return inductance[i][j];
    }
}
//Overload [] operator for a batter user experience.
//Legency.
map<int, complex<double>> &induct::operator[](int i)
{
    return inductance[i];
}
//A tree generator (or plan to be).
//I soon realized that it's not necessary.
//Now it will find the start point of the network.
//Use a typical tree generate algorithm.
inline void induct::gen_tree()
{
    int *tree_book = new int[node_number]{0};
    for (int i = 0; i < node_number; ++i)
    {
        for (int j = 0; j < node_number; ++j)
        {
            if (book[i][j])
            {
                ++tree_book[i];
            }
        }
    }
    for (start = 0; start < node_number; ++start)
    {
        if (tree_book[start] == 1)
        {
            break;
        }
    }
    if (start >= node_number)
        start = node_number - 1;
    //    delete[] tree_book;

    return;
}
//The main generate function of the class.
//ind: the resistance of the line
//a,b: 2 ends of the line
//example :
//  ind.add_line({1,0},1,1);
void induct::add_line_with_transformer(complex<double> &ind, int f, int t, double k)
{
    ind = ind / k;
    add_line(ind, f, t);

    return;
}

int induct::add_line(complex<double> &ind, int a, int b)
{
    inductance[a][b] -= 1.0 / ind;
    inductance[b][a] -= 1.0 / ind;
    inductance[a][a] += 1.0 / ind;
    inductance[b][b] += 1.0 / ind;
    book[a][b] = book[b][a] = 1; //! self-self will not be counted.
    return 0;
}
induct::induct(int node_number) : node_number(node_number)
{
    book = new bool *[node_number];
    inductance = new map<int, complex<double>>[node_number];
    flow = new complex<double> *[node_number];

    for (int i = 0; i < node_number; ++i)
    {
        //inductance[i] = new complex<double> [node_number];
        flow[i] = new complex<double>[node_number];
        book[i] = new bool[node_number]{0};
        //memset(book[i], 0, node_number * sizeof(bool));
    }
}

induct::~induct()
{
    //try {
    //    for (int i = 0; i < node_number; ++i)
    //    {
    //        delete[] inductance[i];
    //    }
    //    delete[] inductance;
    //    for (int i = 0; i < node_number; ++i)
    //    {
    //        delete[] flow[i];
    //    }
    //    delete[] flow;
    //    for (int i = 0; i < node_number; ++i)
    //    {
    //        delete[] book[i];
    //    }
    //    delete[] book;
    //}
    //catch (const std::exception &  e) {}
}
