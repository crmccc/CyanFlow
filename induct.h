#include <complex>
#ifndef MAX_NODE_NUMBER
#define MAX_NODE_NUMBER 100
#endif
using std::complex;
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
    complex<double> inductance[MAX_NODE_NUMBER][MAX_NODE_NUMBER];
    int add_line(complex<double> i,int a,int b);
    complex<double>* operator[](int i);
    induct(/* args */);
    ~induct();
};
complex<double>* induct::operator[](int i){
    return inductance[i];
}
int induct::add_line(complex<double> i,int a,int b){
    inductance[a][b]-=i;
    inductance[b][a]=inductance[a][b];
    inductance[a][a]+=i;
    inductance[b][b]+=i;

    return 0;
}
induct::induct(/* args */)
{
}

induct::~induct()
{
}
