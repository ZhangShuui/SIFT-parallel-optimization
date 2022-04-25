#include <iostream>
#include "GuassDePyramid.h"
#include "GaussDePyramid-pThread.h"
#include "GaussDePyramid-NEON.h"
#include "GaussDePyramid-NEONxPTHREAD.h"
#include <sys/time.h>
using namespace std;
const int MAX=4096;
int n=8;
timeval start, final;
int main() {
    int** p=new int* [MAX];
    for (int i = 0; i < MAX; ++i) {
        p[i]=new int[MAX];
    }
    for (int i = 0; i < MAX; ++i) {
        for (int j = 0; j < MAX; ++j) {
            p[i][j]=1;
        }
    }
    for (int n = 32; n <= MAX; n*=2) {
        int count=0;
        gettimeofday(&start,NULL);
        gettimeofday(&final,NULL);
        while (final.tv_sec-start.tv_sec<10.0){
            count++;
            GaussPyramid_nxp g(p, n, 2);
            g.GenerateDoG();
            gettimeofday(&final,NULL);
        }
        cout<<n<<","<<(final.tv_sec-start.tv_sec)/double (count)<<endl;
    }
    GaussPyramid g(p ,n ,2);
    g.GenerateDoG();
    g.output();
    return 0;
}
