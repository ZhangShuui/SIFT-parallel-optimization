#include <iostream>
#include "GuassDePyramid.h"
//#include "GaussDePyramid-NEON.h"
#include "GaussDePyramid-AVXxOpenMP.h"
#include "GaussDePyramid-OpenMP.h"
#include "GaussDePyramid-AVX512xOpenMP.h"
//#include "GaussDePyramid-NEONxOpenMP.h"
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
    float total2 =0.0;
    int count2 =0;
    float total3 =0.0;
    int count3 =0;
    n = 8;
    int k = 1;
    while (n<=4096){
        GaussPyramid_aomp g_omp(p, n, 2);
        GaussPyramid_omp g(p,n,2);
        GaussPyramid_a512omp ga(p,n,2);
        total =0.0;
        count =0;
        while (total< 6.0){
            g_omp.GaussPyInit();
            gettimeofday(&start, NULL);
            g_omp.GenerateDoG_nomp_dynamic();
            gettimeofday(&final, NULL);
            total += final.tv_sec - start.tv_sec;
            count++;
        }

        total2 =0.0;
        count2 =0;
        while (total2< 6.0){
            g.GaussPyInit();
            gettimeofday(&start, NULL);
            g.GenerateDoG_omp();
            gettimeofday(&final, NULL);
            total2 += final.tv_sec - start.tv_sec;
            count2++;
        }

        total3 =0.0;
        count3 =0;
        while (total3< 6.0){
            ga.GaussPyInit();
            gettimeofday(&start, NULL);
            ga.GenerateDoG_nomp_dynamic();
            gettimeofday(&final, NULL);
            total3 += final.tv_sec - start.tv_sec;
            count3++;
        }
        std::cout <<n<<","<<total2 / float(count2)<<","  << total / float(count) <<","<<total3 / float(count3)<<endl;
        n*=2;
//        k++;
    }
//    GaussPyramid_omp g(p,16,2);
//    g.GaussPyInit();
//    g.GenerateDoG_omp();
//    g.output();
//    GaussPyramid gg(p,8,2);
//    gg.GaussPyInit();
//    gg.GenerateDoG();
//    gg.output();
//    GaussPyramid_aomp ggg(p,8,2);
//    ggg.GaussPyInit();
//    ggg.GenerateDoG_nomp_dynamic();
//    ggg.output();
//    GaussPyramid_a512omp ggg(p,16,2);
//    ggg.GaussPyInit();
//    ggg.GenerateDoG_nomp_dynamic();
//    ggg.output();
    return 0;
}
