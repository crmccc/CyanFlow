#ifndef _GLIBCXX_MAP
#include <map>
#endif
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#include<cstring>
#include<vector>
using std::vector;
using std::complex;
using std::map;

//!THE FIRST NODE IS MARKED AS 0,NOT 1!
/*
*
* inductance ---> node
*
* --data: a list of maps of int and complex double.
*
*/
class induct
{
private:
    /* data */
public:
    
    // bool **book;
    bool book[100][100]{0};
    int node_number;
    // vector<map<int, complex<double>>> inductance;
    // complex<double>** inductance;
    complex<double> inductance[100][100];
    
    int add_line(complex<double>& , int , int );
    complex<double>* operator[](int i);

    induct(int);
    ~induct();
};
complex<double>* induct::operator[](int i)
{
    return inductance[i];
}
int induct::add_line(complex<double>& ind, int a, int b)
{
    inductance[a][b] -= ind;
    inductance[b][a] = inductance[a][b];
    inductance[a][a] += ind;
    inductance[b][b] += ind;
    book[a][b] = book[b][a] = 1; //! self-self will not be counted.
    return 0;
}
induct::induct(int node_number) : node_number(node_number)
{
    // book = new bool* [node_number];
    // inductance = new complex<double> *[node_number];
    // for (int i = 0; i < node_number; ++i)
    {
        // inductance[i]=new complex<double> [node_number];
        // book[i] = new bool[node_number]{0};
        // memset(book[i], 0, node_number * sizeof(bool));
    }
}

induct::~induct()
{
    // for (int i = 0; i < node_number; ++i)
    // {
        // delete[] book[i];
    // }
    // delete[] book;

}
