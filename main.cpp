#include <iostream>
#include "GuassDePyramid.h"
using namespace std;
int n=4;
int main() {
    int** p=new int* [n];
    for (int i = 0; i < n; ++i) {
        p[i]=new int[n];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            p[i][j]=i+j;
        }
    }
    GaussPyramid g(p,n,2);
    g.GenerateDoG();
    g.output();
    return 0;
}
