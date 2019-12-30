#include <complex>
using namespace std;

class Matrix
{
private:
    /* data */
public:
    complex<double> **data;
    int x_length=0;
    int y_length=0;
    Matrix(int, int);
    ~Matrix();
};
//x----
//y
//|
//|
Matrix::Matrix(int x, int y)
{
    data = new complex<double> *[x];
    for (int i = 0; i < y; ++i)
    {
        data[i] = new complex<double>[y];
    }
    x_length=x;
    y_length=y;
}

Matrix::~Matrix()
{
    for (int x=0;x<x_length;++x){
        delete[] data[x];
    }
    delete[] data;
}
