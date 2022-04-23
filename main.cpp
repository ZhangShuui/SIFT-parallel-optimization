#include <iostream>
#include "GuassDePyramid.h"
#include <sys/time.h>
using namespace std;
const int MAX=4096;
int n=128;
timeval start, final;
int main() {
    int** p=new int* [MAX];
    for (int i = 0; i < MAX; ++i) {
        p[i]=new int[MAX];
    }
    for (int i = 0; i < MAX; ++i) {
        for (int j = 0; j < MAX; ++j) {
            p[i][j]=i+j;
        }
    }

    for (int n = 32; n <= MAX; n*=2) {
        int count=0;
        gettimeofday(&start,NULL);
        gettimeofday(&final,NULL);
        while (final.tv_sec-start.tv_sec<5.0){
            count++;
            GaussPyramid g(p, n, 2);
            g.GenerateDoG();
            gettimeofday(&final,NULL);
        }
        cout<<n<<","<<(final.tv_sec-start.tv_sec)/(1000.0*count)<<endl;
    }
    return 0;
}
