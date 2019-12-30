#include <complex>
using std::complex;

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
    Matrix(int, int);
    Matrix operator*(const Matrix &);
    void operator=(const Matrix &);
    complex<double>*& operator[](int i);
    
    ~Matrix();
};

//x----
//y
//|
//|
//also fixed length
void Matrix::operator=(const Matrix &b){
    for (int x = 0; x < x_length; ++x)
    {
        delete[] data[x];
    }
    delete[] data;
    data=b.data;
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
    return res;
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
    x_length = x;
    y_length = y;
    data = make_a_matrix(x, y);
    y_length = y;

}

Matrix::~Matrix()
{
    for (int x = 0; x < x_length; ++x)
    {
        delete[] data[x];
    }
    delete[] data;
    data=NULL;
}