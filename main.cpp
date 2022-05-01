#include <iostream>
#include "GuassDePyramid.h"
#include "GaussDePyramid-pThread.h"
//#include "GaussDePyramid-SSExPTHREAD.h"
//#include "GaussDePyramid-AVX512xPTHREAD.h"
//#include "GaussDePyramid-NEON.h"
//#include "GaussDePyramid-NEONxPTHREAD.h"
#include <sys/time.h>

using namespace std;

const int MAX = 4096;
int n = 8;

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
    GaussPyramid_p g(p, n, 2);
    GaussPyramid_p::THREAD_COUNT = 5;

    g.GaussPyInit();
    gettimeofday(&start, NULL);
    g.GenerateDoG();
    gettimeofday(&final, NULL);
    total += final.tv_sec - start.tv_sec;
    count++;

    std::cout << GaussPyramid_p::THREAD_COUNT << "," << total / float(count) << std::endl;

//
    return 0;
}
