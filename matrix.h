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
    Matrix operator*(const Matrix &b);
    ~Matrix();
};

//x----
//y
//|
//|
//also fixed length
Matrix Matrix::operator*(const Matrix &bM) 
{    
    Matrix res(bM.y_length, x_length);
    complex<double> **b = bM.data;
    complex<double> **a = data;
    complex<double> **c = res.data;
    for (int i = 0; i < x_length; i++)
    {
        for (int j = 0; j < bM.y_length; j++)
        {
            for (int k = 1; k <= y_length; a++)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return res;
}

complex<double> **Matrix::make_a_matrix(int x, int y)
{
    complex<double> **res = new complex<double> *[x];
    for (int i = 0; i < y; ++i)
    {
        res[i] = new complex<double>[y];
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

Matrix::Matrix(int x, int y)
{
    data = make_a_matrix(x, y);
}

Matrix::~Matrix()
{
    for (int x = 0; x < x_length; ++x)
    {
        delete[] data[x];
    }
    delete[] data;
}
