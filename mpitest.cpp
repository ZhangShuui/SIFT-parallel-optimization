//
// Created by zsr on 2022/6/15.
//



#include <iostream>
#include "GuassDePyramid.h"

#include <mpi.h>
#include <omp.h>
//#include <immintrin.h>
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
int thread_count = 4;
int chunk_size = 5;
const int MAX = 4096;
int n = 128;
int S=2;
float **** GaussPy;
int length=n;
int layer=0;
bool is_initialized= false;
void GenerateDoG_mpi_omp(int argc,char** argv) {
    int i,len=length,j=0,jj,m,n,k,MyLayer,i1;
    float l,sig;
    double begin,final;
    MPI_Init(&argc, &argv);
    begin = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    if (i < S+3) {
        while (len) {
            l = float(len - 1) / 2.0;
            sig = sigma / float (i + 1);
            {
#pragma omp for
                for (m = 0; m < len; ++m) {
                    for (n = 0; n < len; ++n) {
                        GaussPy[j][i][m][n] *= exp(-(n - l) * (n - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
                    }
                }
#pragma omp for
                for (m = 0; m < len; ++m) {
                    for (n = 0; n < len; ++n) {
                        GaussPy[j][i][m][n] *= exp(-(m - l) * (m - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
                    }
                    MPI_Send(GaussPy[j][i][m],len,MPI_FLOAT,S+3, j*length* 10+i*length+m,MPI_COMM_WORLD);
                }
            }
            j++;
            len /= 2;
        }
    }
    else if (i == S + 3) {
        for (j = 0; j < S+3; j++) {
            len=length;
            m=0;
            while (len) {
                for (i1 = 0; i1 < len; ++i1) {
                    MPI_Recv(GaussPy[m][j][i1], len, MPI_FLOAT, j, m * length * 10 + j * length + i1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                }
                len /= 2;
                m++;
            }
        }
        for (MyLayer = 0; MyLayer < layer; MyLayer++) {
            len = length;
            k = MyLayer;
            while (k) {
                k--;
                len /= 2;
            }
            for (jj = 0; jj < S + 2; ++jj) {
#pragma omp for
                for (j = 0; j < len; ++j) {
                    for (m = 0; m < len; ++m) {
                        GaussPy[MyLayer][jj][j][m] -= GaussPy[MyLayer][jj + 1][j][m];
                    }
                }
            }
        }
        len = length;
        final = MPI_Wtime();
        cout<<final- begin<<endl;
//        for (jj = 0; jj < layer; ++jj) {
//            for (j = 0; j < len; ++j) {
//                for (k = 0; k < len; ++k) {
//                    std::cout << GaussPy[jj][0][j][k] << " ";
//                }
//                std::cout << std::endl;
//            }
//            for (k = 0; k < len; ++k) {
//                std::cout << "==";
//            }
//            std::cout << std::endl;
//            len /= 2;
//        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
void GenerateDoG_mpi(int argc,char** argv) {
    int i,len=length,j=0,jj,m,n,k,MyLayer,i1;
    float l,sig;
    double begin,final;
    MPI_Init(&argc, &argv);
    begin = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    if (i < S+3) {
        while (len) {
            l = float(len - 1) / 2.0;
            sig = sigma / float (i + 1);
            {
                for (m = 0; m < len; ++m) {
                    for (n = 0; n < len; ++n) {
                        GaussPy[j][i][m][n] *= exp(-(n - l) * (n - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
                    }
                }
                for (m = 0; m < len; ++m) {
                    for (n = 0; n < len; ++n) {
                        GaussPy[j][i][m][n] *= exp(-(m - l) * (m - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
                    }
                    MPI_Send(GaussPy[j][i][m],len,MPI_FLOAT,S+3, j*length* 10+i*length+m,MPI_COMM_WORLD);
                }
            }
            j++;
            len /= 2;
        }
    }
    else if (i == S + 3) {
        for (j = 0; j < S+3; j++) {
            len=length;
            m=0;
            while (len) {
                for (i1 = 0; i1 < len; ++i1) {
                    MPI_Recv(GaussPy[m][j][i1], len, MPI_FLOAT, j, m * length * 10 + j * length + i1, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                }
                len /= 2;
                m++;
            }
        }
        for (MyLayer = 0; MyLayer < layer; MyLayer++) {
            len = length;
            k = MyLayer;
            while (k) {
                k--;
                len /= 2;
            }
            for (jj = 0; jj < S + 2; ++jj) {
                for (j = 0; j < len; ++j) {
                    for (m = 0; m < len; ++m) {
                        GaussPy[MyLayer][jj][j][m] -= GaussPy[MyLayer][jj + 1][j][m];
                    }
                }
            }
        }
        len = length;
        final = MPI_Wtime();
        cout<<final- begin<<endl;
//        for (jj = 0; jj < layer; ++jj) {
//            for (j = 0; j < len; ++j) {
//                for (k = 0; k < len; ++k) {
//                    std::cout << GaussPy[jj][0][j][k] << " ";
//                }
//                std::cout << std::endl;
//            }
//            for (k = 0; k < len; ++k) {
//                std::cout << "==";
//            }
//            std::cout << std::endl;
//            len /= 2;
//        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
//void GenerateDoG_mpi_SSE(int argc,char** argv) {
//    int i,len=length,j=0,jj,m,n,k,MyLayer,i1,kk,ll;
//    __m128 va,vf;
//    double begin,final;
//    float* filter = new float [length];
//    float l,sig;
//    MPI_Init(&argc, &argv);
//    begin = MPI_Wtime();
//    MPI_Comm_rank(MPI_COMM_WORLD, &i);
//    if (i < S+3) {
//        while (len) {
//            if (len<=2){
//                    l = float(len - 1) / 2.0;
//                    sig = sigma / float (i + 1);
//                    {
//                        for (n = 0; n < len; ++n) {
//                            filter[n] = exp(-(n - l) * (n - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
//                        }
//                        for (m = 0; m < len; ++m) {
//                            for (n = 0; n < len; ++n) {
//                                GaussPy[j][i][m][n] *= filter[n];
//                            }
//                        }
//                        for (m = 0; m < len; ++m) {
//                            for (n = 0; n < len; ++n) {
//                                GaussPy[j][i][m][n] *= filter[m];
//                            }
//                            MPI_Send(GaussPy[j][i][m],len,MPI_FLOAT,S+3, j*length* 10+i*length+m,MPI_COMM_WORLD);
//                        }
//                    }
//                    j++;
//                    len /= 2;
//                } else {
//                l = float(len - 1) / 2.0;
//                sig = sigma / float (i + 1);
//                for (n = 0; n < len; ++n) {
//                    filter[n] = exp(-(n - l) * (n - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
//                }
//                for (jj = 0; jj < len; jj+=4) {
//                    vf = _mm_loadu_ps(filter+jj);
//                    for (kk = 0; kk < len; ++kk) {
//                        va = _mm_loadu_ps(GaussPy[j][i][kk] + jj);
//                        va = _mm_mul_ps(va,vf);
//                        _mm_store_ps(GaussPy[j][i][kk] + jj, va);
//                    }
//                }
//                for (jj = 0; jj < len; jj++) {
//                    vf = _mm_set1_ps(filter[jj]);
//                    for (kk = 0; kk < len; kk+=4) {
//                        va = _mm_loadu_ps(GaussPy[j][i][jj] + kk);
//                        va = _mm_mul_ps(va,vf);
//                        _mm_store_ps(GaussPy[j][i][jj] + kk, va);
//                    }
//                    MPI_Send(GaussPy[j][i][jj],len,MPI_FLOAT,S+3, j*length* 10+i*length+jj,MPI_COMM_WORLD);
//                }
//                j++;
//                len /= 2;
//            }
//        }
//    }
//    else if (i == S + 3) {
//        for (j = 0; j < S+3; j++) {
//            len=length;
//            m=0;
//            while (len) {
//                for (i1 = 0; i1 < len; ++i1) {
//                    MPI_Recv(GaussPy[m][j][i1], len, MPI_FLOAT, j, m * length * 10 + j * length + i1, MPI_COMM_WORLD,
//                             MPI_STATUS_IGNORE);
//                }
//                len /= 2;
//                m++;
//            }
//        }
//        for (MyLayer = 0; MyLayer < layer; MyLayer++) {
//            len = length;
//            k = MyLayer;
//            while (k) {
//                k--;
//                len /= 2;
//            }
//            if (len<=2)
//                for (jj = 0; jj < S + 2; ++jj) {
//                    for (j = 0; j < len; ++j) {
//                        for (m = 0; m < len; ++m) {
//                            GaussPy[MyLayer][jj][j][m] -= GaussPy[MyLayer][jj + 1][j][m];
//                        }
//                    }
//                }
//            else {
//                for (j = 0; j < S + 2; ++j) {
//                    for (kk = 0; kk < len; ++kk) {
//                        for (ll = 0; ll < len; ll+=4) {
//                            va = _mm_loadu_ps(GaussPy[MyLayer][j][kk] + ll);
//                            vf = _mm_loadu_ps(GaussPy[MyLayer][j+1][kk] + ll);
//                            va = _mm_sub_ps(va,vf);
//                            _mm_store_ps(GaussPy[MyLayer][j][kk] + ll, va);
//                        }
//                    }
//                }
//            }
//        }
//        len = length;
//        final = MPI_Wtime();
//        cout<<final- begin<<endl;
////        for (jj = 0; jj < layer; ++jj) {
////            for (j = 0; j < len; ++j) {
////                for (k = 0; k < len; ++k) {
////                    std::cout << GaussPy[jj][0][j][k] << " ";
////                }
////                std::cout << std::endl;
////            }
////            for (k = 0; k < len; ++k) {
////                std::cout << "==";
////            }
////            std::cout << std::endl;
////            len /= 2;
////        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    MPI_Finalize();
//}
//void GenerateDoG_mpi_AVX(int argc,char** argv) {
//    int i,len=length,j=0,jj,m,n,k,MyLayer,i1,kk,ll;
//    __m256 va,vf;
//    double begin,final;
//    float* filter = new float [length];
//    float l,sig;
//    MPI_Init(&argc, &argv);
//    begin = MPI_Wtime();
//    MPI_Comm_rank(MPI_COMM_WORLD, &i);
//    if (i < S+3) {
//        while (len) {
//            if (len<=4){
//                l = float(len - 1) / 2.0;
//                sig = sigma / float (i + 1);
//                {
//                    for (n = 0; n < len; ++n) {
//                        filter[n] = exp(-(n - l) * (n - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
//                    }
//                    for (m = 0; m < len; ++m) {
//                        for (n = 0; n < len; ++n) {
//                            GaussPy[j][i][m][n] *= filter[n];
//                        }
//                    }
//                    for (m = 0; m < len; ++m) {
//                        for (n = 0; n < len; ++n) {
//                            GaussPy[j][i][m][n] *= filter[m];
//                        }
//                        MPI_Send(GaussPy[j][i][m],len,MPI_FLOAT,S+3, j*length* 10+i*length+m,MPI_COMM_WORLD);
//                    }
//                }
//                j++;
//                len /= 2;
//            } else {
//                l = float(len - 1) / 2.0;
//                sig = sigma / float (i + 1);
//                for (n = 0; n < len; ++n) {
//                    filter[n] = exp(-(n - l) * (n - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
//                }
//                for (jj = 0; jj < len; jj+=8) {
//                    vf = _mm256_load_ps(filter+jj);
//                    for (kk = 0; kk < len; ++kk) {
//                        va = _mm256_load_ps(GaussPy[j][i][kk] + jj);
//                        cout<<"i="<<i<<",here2"<<endl;
//                        va = _mm256_mul_ps(va,vf);
//                        cout<<"i="<<i<<",here"<<endl;
//                        _mm256_store_ps(GaussPy[j][i][kk] + jj, va);
//                    }
//                }
//
//                for (jj = 0; jj < len; jj++) {
//                    vf = _mm256_set1_ps(filter[jj]);
//                    for (kk = 0; kk < len; kk+=8) {
//                        va = _mm256_load_ps(GaussPy[j][i][jj] + kk);
//                        va = _mm256_mul_ps(va,vf);
//                        _mm256_store_ps(GaussPy[j][i][jj] + kk, va);
//                    }
//                    MPI_Send(GaussPy[j][i][jj],len,MPI_FLOAT,S+3, j*length* 10+i*length+jj,MPI_COMM_WORLD);
//                }
//                j++;
//                len /= 2;
//            }
//        }
//
//    }
//    else if (i == S + 3) {
//        for (j = 0; j < S+3; j++) {
//            len=length;
//            m=0;
//            while (len) {
//                for (i1 = 0; i1 < len; ++i1) {
//                    MPI_Recv(GaussPy[m][j][i1], len, MPI_FLOAT, j, m * length * 10 + j * length + i1, MPI_COMM_WORLD,
//                             MPI_STATUS_IGNORE);
//                }
//                len /= 2;
//                m++;
//            }
//        }
//        for (MyLayer = 0; MyLayer < layer; MyLayer++) {
//            len = length;
//            k = MyLayer;
//            while (k) {
//                k--;
//                len /= 2;
//            }
//            if (len<=4)
//                for (jj = 0; jj < S + 2; ++jj) {
//                    for (j = 0; j < len; ++j) {
//                        for (m = 0; m < len; ++m) {
//                            GaussPy[MyLayer][jj][j][m] -= GaussPy[MyLayer][jj + 1][j][m];
//                        }
//                    }
//                }
//            else {
//                for (j = 0; j < S + 2; ++j) {
//                    for (kk = 0; kk < len; ++kk) {
//                        for (ll = 0; ll < len; ll+=8) {
//                            va = _mm256_load_ps(GaussPy[MyLayer][j][kk] + ll);
//                            vf = _mm256_load_ps(GaussPy[MyLayer][j+1][kk] + ll);
//                            va = _mm256_sub_ps(va,vf);
//                            _mm256_store_ps(GaussPy[MyLayer][j][kk] + ll, va);
//                        }
//                    }
//                }
//            }
//        }
//        len = length;
//        final = MPI_Wtime();
//        cout<<final- begin<<endl;
////        for (jj = 0; jj < layer; ++jj) {
////            for (j = 0; j < len; ++j) {
////                for (k = 0; k < len; ++k) {
////                    std::cout << GaussPy[jj][0][j][k] << " ";
////                }
////                std::cout << std::endl;
////            }
////            for (k = 0; k < len; ++k) {
////                std::cout << "==";
////            }
////            std::cout << std::endl;
////            len /= 2;
////        }
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    MPI_Finalize();
//}
void GaussPyInit(int* data[MAX]) {
    int step = 1;
    length=n;
    int ll=length;
    if (!is_initialized){
        while (ll){
            layer++;
            ll/=2;
        }
        GaussPy=new float***[layer];
        for (int i = 0; i < layer; ++i) {
            GaussPy[i] = new float** [S + 3];
            for (int j = 0; j < S + 3; ++j) {
                GaussPy[i][j] = new float* [length / step];
                for (int k = 0; k < length / step; ++k) {
                    GaussPy[i][j][k] = new float[length / step];
                }
            }
            step *= 2;
        }
    }
    is_initialized = true;
    int len = length;
    step = 1;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            for (int k = 0; k < len; ++k) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[i][j][k][l] = float (data[k * step][l * step]);
                }
            }
        }
        len /= 2;
        step *= 2;
    }
}
void delete_mpi() {
    int step = 1;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            for (int k = 0; k < length / step; ++k) {
                delete[]GaussPy[i][j][k];
            }
        }
        step *= 2;
    }
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            delete[]GaussPy[i][j];
        }
    }
    for (int i = 0; i < layer; ++i) {
        delete[]GaussPy[i];
    }
    delete[]GaussPy;
}


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
//    int times=0;
//    std::chrono::duration<double, std::milli> elapsed{};
//    std::chrono::duration<double, std::milli> elapsed1{};
//    auto start= std::chrono::high_resolution_clock::now();
//    auto end = std::chrono::high_resolution_clock::now();
//    GaussPyInit(p);
//    start= std::chrono::high_resolution_clock::now();
//    GenerateDoG_mpi_omp(argc,argv);
//    end = std::chrono::high_resolution_clock::now();
//    elapsed1 += end-start;
//    cout<<float (elapsed1.count())<<endl;
//    GaussPyInit(p);
//    start= std::chrono::high_resolution_clock::now();
//    GenerateDoG_mpi(argc,argv);
//    end = std::chrono::high_resolution_clock::now();
//    elapsed += end-start;
//    delete_mpi();
//    cout<<float (elapsed.count())<<endl;
    n=256;
//    GaussPyInit(p);
//    GenerateDoG_mpi_AVX(argc,argv);
//    GaussPyInit(p);
//    GenerateDoG_mpi(argc,argv);
    GaussPyInit(p);
    GenerateDoG_mpi_omp(argc,argv);


    return 0;
}

