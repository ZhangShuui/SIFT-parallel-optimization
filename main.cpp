#include <iostream>
#include "GuassDePyramid.h"
#include "GaussDePyramid-MPI.h"

//#include "GaussDePyramid-NEON.h"
//#include "GaussDePyramid-AVXxOpenMP.h"
//#include "GaussDePyramid-OpenMP.h"
//#include "GaussDePyramid-AVX512xOpenMP.h"
//#include "GaussDePyramid-NEONxOpenMP.h"
//#include "GaussDePyramid-SSExPTHREAD.h"
//#include "GaussDePyramid-AVX512xPTHREAD.h"
//#include "GaussDePyramid-NEON.h"
//#include "GaussDePyramid-NEONxPTHREAD.h"


using namespace std;
//timeval start,final;
const int MAX = 1024;
int n = 512;






int main(int argc,char* argv[]) {
    int **p = new int *[MAX];
    for (int i = 0; i < MAX; ++i) {
        p[i] = new int[MAX];
    }
    for (int i = 0; i < MAX; ++i) {
        for (int j = 0; j < MAX; ++j) {
            p[i][j] = 1;
        }
    }
//    float total = 0.0;
//    int count = 0;
//    float total2 =0.0;
//    int count2 =0;
//    float total3 =0.0;
//    int count3 =0;
//    n = 8;
//    int k = 15;
//    while (k <= 25){
//        g.chunk_size = k;
//        total2 =0.0;
//        count2 =0;
//        while (total2< 6.0){
//            g.GaussPyInit();
//            gettimeofday(&start, NULL);
//            g.GenerateDoG_omp_dynamic();
//            gettimeofday(&final, NULL);
//            total2 += final.tv_sec - start.tv_sec;
//            count2++;
//        }
//        cout<<k<<","<<total2/float (count2)<<endl;
////        n*=2;
//        k++;
//    }
    int times=0;
    GaussPyramid_mpi g(p,n,2);
    std::chrono::duration<double, std::milli> elapsed{};
    auto start= std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    elapsed +=end-start;
    while(elapsed.count()<100){
        start= std::chrono::high_resolution_clock::now();
        g.GenerateDoG_mpi(argc,argv);
        end = std::chrono::high_resolution_clock::now();
        elapsed += end-start;

        times+=1;
    }
    cout<<float (elapsed.count())/ float (times)<<endl;
    return 0;
}
