#ifndef _GLIBCXX_MAP
#include <map>
#endif
#ifndef _GLIBCXX_COMPLEX
#include <complex>
#endif
#include<cstring>
#include<vector>
//#include<list>
using std::vector;
using std::complex;
//using std::list;

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
    // bool book[100][100]{0};
    int node_number;
    // vector<map<int, complex<double>>> inductance;
     complex<double>** inductance;
     complex<double>** flow;
     int start,end;
    //complex<double> inductance[100][100];
    
    int add_line(complex<double>& , int , int );
    complex<double>* operator[](int i);
    void gen_tree();

    induct(int);
    ~induct();
};
complex<double>* induct::operator[](int i)
{
    return inductance[i];
}
inline void induct::gen_tree()
{
    int* tree_book = new int[node_number] {0};
    for (int i = 0;i < node_number;++i) {
        for (int j = 0;j < node_number;++j) {
            if (book[i][j]) {
                ++tree_book[i];
            }
        }
    }
    for (start=0;start < node_number;++start) {
        if (tree_book[start] == 1) {
            break;
        }
    }
    for (end = 0;end < node_number;++end) {
        if (tree_book[end] == 1&&end!=start) {
            break;
        }
    }
    delete[] tree_book;

    return;
}
int induct::add_line(complex<double>& ind, int a, int b)
{
    inductance[a][b] -= 1.0/ind;
    inductance[b][a] -= 1.0/ind;
    inductance[a][a] += 1.0/ind;
    inductance[b][b] += 1.0/ind;
    book[a][b] = book[b][a] = 1; //! self-self will not be counted.
    return 0;
}
induct::induct(int node_number) : node_number(node_number)
{
     book = new bool* [node_number];
     inductance = new complex<double> *[node_number];
     flow = new complex<double> * [node_number];

     for (int i = 0; i < node_number; ++i)
    {
         inductance[i] = new complex<double> [node_number];
         flow[i] = new complex<double>[node_number];
         book[i] = new bool[node_number]{0};
         //memset(book[i], 0, node_number * sizeof(bool));
    }
}

induct::~induct()
{
     for (int i = 0; i < node_number; ++i)
     {
         delete[] inductance[i];
     }
     delete[] inductance;
     for (int i = 0; i < node_number; ++i)
     {
         delete[] flow[i];
     }
     delete[] flow;
     for (int i = 0; i < node_number; ++i)
     {
         delete[] book[i];
     }
     delete[] book;

}
