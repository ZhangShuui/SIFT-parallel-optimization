#include <iostream>
#include "GuassDePyramid.h"
//#include "GaussDePyramid-NEON.h"
#include "GaussDePyramid-OpenMP.h"
#include "GaussDePyramid-NEONxOpenMP.h"
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
    n = 128;
    int k = 1;
    while (n<=16384){
        GaussPyramid_nomp g_omp(p, n, 2);
        GaussPyramid_omp g(p,n,2);
        total =0.0;
        count =0;
        while (total< 10.0){
            g_omp.GaussPyInit();
            gettimeofday(&start, NULL);
            g_omp.GenerateDoG_nomp_dynamic();
            gettimeofday(&final, NULL);
            total += final.tv_sec - start.tv_sec;
            count++;
        }

        float total2 =0.0;
        int count2 =0;
        while (total2< 10.0){
            g.GaussPyInit();
            gettimeofday(&start, NULL);
            g.GenerateDoG_omp();
            gettimeofday(&final, NULL);
            total2 += final.tv_sec - start.tv_sec;
            count2++;
        }
        std::cout <<n<<","  << total / float(count) <<","<<total2 / float(count2)<<endl;
        n*=2;
//        k++;
    }
//    GaussPyramid_nomp g(p,8,2);
//    g.GaussPyInit();
//    g.GenerateDoG_nomp_dynamic();
//    g.output();
//    GaussPyramid gg(p,8,2);
//    gg.GaussPyInit();
//    gg.GenerateDoG();
//    gg.output();
//    GaussPyramid_nomp g(p,8,2);
//    g.GaussPyInit();
//    g.GenerateDoG_nomp_dynamic();
//    g.output();
    return 0;
}
