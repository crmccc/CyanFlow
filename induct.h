#ifndef _GLIBCXX_MAP
#include <map>
#endif
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#include <map>
#ifndef MAX_NODE_NUMBER
#define MAX_NODE_NUMBER 100
#endif

using std::complex;
using std::map;
//!IMPORTANT!
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
    map<int, complex<double>> inductance[];
    bool **book;
    int node_number;

    int add_line(complex<double> i, int a, int b);

    map<int, complex<double>> &operator[](int i);

    induct(int);
    ~induct();
};
map<int, complex<double>> &induct::operator[](int i)
{
    return inductance[i];
}
int induct::add_line(complex<double> i, int a, int b)
{
    inductance[a][b] -= i;
    inductance[b][a] = inductance[a][b];
    inductance[a][a] += i;
    inductance[b][b] += i;
    book[a][b] = book[b][a] = 1; //! self-self will not be counted.
    return 0;
}
induct::induct(int node_number) : node_number(node_number)
{
    book = new bool *[node_number];
    for (int i = 0; i < node_number; ++i)
    {
        book[i] = new bool[node_number];
        memset(book[i], 0, node_number * sizeof(bool));
    }
}

induct::~induct()
{
    for (int i = 0; i < node_number; ++i)
    {
        delete[] book[i];
    }
    delete[] book;
}
