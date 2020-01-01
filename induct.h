#ifndef _GLIBCXX_MAP
#include <map>
#endif
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#include<cstring>
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
    
    bool **book;
    int node_number;
    map<int, complex<double>>* inductance;

    int add_line(complex<double> i, int a, int b);
    map<int, complex<double>> &operator[](int i);

    induct(int);
    ~induct();
};
map<int, complex<double>> &induct::operator[](int i)
{
    return inductance[i];
}
int induct::add_line(complex<double> ind, int a, int b)
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
    inductance = new map<int,complex<double>> [node_number];
    book = new bool* [node_number];
    for (int i = 0; i < node_number; ++i)
    {
        book[i] = new bool[node_number]{0};
        // memset(book[i], 0, node_number * sizeof(bool));
    }
}

induct::~induct()
{
    for (int i = 0; i < node_number; ++i)
    {
        delete[] book[i];
    }
    delete[] book;
    delete[] inductance;

}
