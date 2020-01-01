#ifndef _GLIBCXX_MAP
#include<map>
#endif
#ifndef _GLIBCXX_COMPLEX
#include<complex>
#endif
#include<map>
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
    map<int,complex<double>> inductance;

    int add_line(complex<double> i,int a,int b);
    
    complex<double>& operator[](int i);
    
    induct(/* args */);
    ~induct();
};
complex<double>& induct::operator[](int i){
    return inductance[i];
}
int induct::add_line(complex<double> i,int a,int b){
    inductance[a]-=i;
    inductance[b]=inductance[a];
    inductance[a]+=i;
    inductance[b]+=i;
    return 0;
}
induct::induct(/* args */)
{
}

induct::~induct()
{
}
