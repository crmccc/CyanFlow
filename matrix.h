#include <complex>
using std::complex;
int numbers=0;
class Matrix
//a fix length ugly complex matrix
{
private:
    /* data */
    complex<double> **make_a_matrix(int x, int y);

public:
    complex<double> **data;
    int x_length = 0;
    int y_length = 0;
    int flag;
    Matrix operator*(const Matrix &);
    void operator=(const Matrix &);
    complex<double>*& operator[](int i);
    int mutiplit(const Matrix &,const Matrix &,const Matrix &);
    
    Matrix(int, int);
    ~Matrix();
};

//x----
//y
//|
//|
//also fixed length
void Matrix::operator=(const Matrix &b){
    //this* must bigger than b
    // for (int x = 0; x < x_length; ++x)
    // {
    //     delete[] data[x];
    // }
    // delete[] data;
    for (int x=0;x<x_length;++x){
        for (int y=0;y<y_length;++y){
            data[x][y]=b.data[x][y];
        }
    }
    x_length=b.x_length;
    y_length=b.y_length;
}

Matrix Matrix::operator*(const Matrix &b) 
{    
    Matrix res(x_length,b.y_length);
    for (int i = 0; i < res.x_length; i++)
    {
        for (int j = 0; j < res.y_length; j++)
        {
            for (int k = 0; k < y_length; k++)
            {
                res.data[i][j] = data[i][k] * b.data[k][j]+res.data[i][j];
            }
        }
    }
    return (res);
}
complex<double>*& Matrix::operator[](int i){
    return data[i];
}

complex<double> **Matrix::make_a_matrix(int x, int y)
{
    complex<double> **res = new complex<double> *[x];
    for (int i = 0; i < y; ++i)
    {
        res[i] = new complex<double>[y]{0};
    }

    return res;
}

Matrix::Matrix(int x, int y)
{
    numbers++;
    flag=numbers;
    x_length = x;
    y_length = y;
    data = make_a_matrix(x, y);
}

Matrix::~Matrix()
{
    numbers--;
    for (int x = 0; x < x_length; ++x)
    {
       delete[] data[x];
    }
    delete[] data;
    data=NULL;
}

int Matrix::mutiplit(const Matrix &a,const Matrix &b,const Matrix &res){
    for (int i = 0; i < res.x_length; i++)
    {
        for (int j = 0; j < res.y_length; j++)
        {
            for (int k = 0; k < y_length; k++)
            {
                res.data[i][j] = data[i][k] * b.data[k][j]+res.data[i][j];
            }
        }
    }

}
