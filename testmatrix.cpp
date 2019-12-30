#include"matrix.h"
#include<stdio.h>
int main (){
    Matrix a(100,100);
    Matrix b(100,100);
    Matrix c(100,100);
    c=a*b;
    printf("%f %f\n%f %f\n\n",a.data[0][0].real(),a.data[0][1].real(),a.data[1][0].real(),a.data[1][1].real());
    printf("%f %f\n%f %f\n\n",b.data[0][0].real(),b.data[0][1].real(),b.data[1][0].real(),b.data[1][1].real());
    printf("%f %f\n%f %f\n\n",c.data[0][0].real(),c.data[0][1].real(),c.data[1][0].real(),c.data[1][1].real());
    
    a.~Matrix();
    b.~Matrix();
    c.~Matrix();
    return 0;
}