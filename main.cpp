#include <iostream>
#include "GuassDePyramid.h"
#include "GaussDePyramid-OpenMP.h"
#include <vector>
//#include "GaussDePyramid-SSExPTHREAD.h"
//#include "GaussDePyramid-AVX512xPTHREAD.h"
//#include "GaussDePyramid-NEON.h"
//#include "GaussDePyramid-NEONxPTHREAD.h"
#include <sys/time.h>

using namespace std;
timeval start,final;
const int MAX = 4096;
int n = 8;
vector<int> aaaaa;
int main() {
    int **p = new int *[MAX];
    for (int i = 0; i < MAX; ++i) {
        p[i] = new int[MAX];
    }
    for (int i = 0; i < MAX; ++i) {
        for (int j = 0; j < MAX; ++j) {
            p[i][j] = 1;
        }
    }
    float total = 0.0;
    int count = 0;
    n = 8;
//    while (n<=MAX){
//        GaussPyramid_omp g(p, n, 2);
//        total =0.0;
//        count =0;
//        while (total< 10.0){
//            g.GaussPyInit();
//            gettimeofday(&start, NULL);
//            g.GenerateDoG();
//            gettimeofday(&final, NULL);
//            total += final.tv_sec - start.tv_sec;
//            count++;
//        }
//        float total2 =0.0;
//        int count2 =0;
//        while (total< 10.0){
//            g.GaussPyInit();
//            gettimeofday(&start, NULL);
//            g.GenerateDoG();
//            gettimeofday(&final, NULL);
//            total2 += final.tv_sec - start.tv_sec;
//            count2++;
//        }
//        std::cout  << total / float(count) <<","<<total2 / float(count2);
//        n*=2;
//    }
//
    GaussPyramid_omp g_omp(p,n,2);
    GaussPyramid g(p,n,2);
    g.GenerateDoG();
    g_omp.GenerateDoG();
    g.output();
    g_omp.output();
    cout<<"over!"<<endl;
    return 0;
}
