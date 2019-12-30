#include"matrix.h"
#include<stdio.h>
void show_matrix(Matrix &a);
int main (){
    Matrix a(100,100);
    Matrix b(100,100);
    Matrix c(100,100);
    c=a*b;
    for (int i=0;i<2;++i){
        for(int j=0;j<2;++j){
            a.data[i][j]=i*j+1;
            b.data[i][j]=i*j+1;
        }
    }
    c=a*b;
    show_matrix(c);

    a.~Matrix();
    b.~Matrix();
    c.~Matrix();
    return 0;
}

void show_matrix(Matrix &a){
    for(int i=0;i<2;++i){
        for(int j=0;j<2;++j){
            printf("%f ",a.data[i][j].real());
        }
        printf("\n");
    }
}