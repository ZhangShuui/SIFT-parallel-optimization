//
// Created by zsr on 2022/5/31.
//

#ifndef SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_MPI_H
#define SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_MPI_H

#include "GuassDePyramid.h"
#include <iostream>
#include <mpi.h>
#include <chrono>
#include <math.h>
#include <sys/time.h>


class GaussPyramid_mpi {
public:
    int** data; //记录图像灰度数据
    GaussPyramid_mpi();

    GaussPyramid_mpi(int** img, int len/*,int wid*/, int S);

    float**** GaussPy;

    void GaussPyInit();

    void output();

    void GaussFilter(int theLayer);

    void GenerateDoG();

    void GenerateDoG_mpi_normal();

    void GenerateDoG_mpi(int argc,char** argv);

    ~GaussPyramid_mpi();

    int thread_count;

    int chunk_size;

    int all_time;

protected:
    int length;
    //int width;  先尝试宽度相同的版本
    int S;  //提取图像特征后,需要进行对比的图片数
    int layer;
    float* filter;
    bool is_initialized;

};

GaussPyramid_mpi::GaussPyramid_mpi() {
    data = nullptr;
    is_initialized = false;
}

GaussPyramid_mpi::GaussPyramid_mpi(int** img, int len, int S) {
    length = len;
    is_initialized = false;
    thread_count = 8;
    chunk_size = 5;
    data = new int* [len];
    for (int i = 0; i < len; ++i) {
        data[i] = new int[len];
    }
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            data[i][j] = img[i][j];
        }
    }
    this->S = S;
    int x = 0;
    while (len) {
        x++;
        len /= 2;
    }
    layer = x;
    GaussPy = new float*** [layer];
    filter = new float[length];
    GaussPyInit();
}

//初始化高斯金字塔，尚未进行高斯滤波操作
void GaussPyramid_mpi::GaussPyInit() {
    int step = 1;
    if (!is_initialized)
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
    is_initialized = true;
    int len = length;
    step = 1;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            for (int k = 0; k < len; ++k) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[i][j][k][l] = data[k * step][l * step];
                }
            }
        }
        len /= 2;
        step *= 2;
    }
}

void GaussPyramid_mpi::output() {
    int len = length;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < len; ++j) {
            for (int k = 0; k < len; ++k) {
                std::cout << GaussPy[i][0][j][k] << " ";
            }
            std::cout << std::endl;
        }
        for (int k = 0; k < len; ++k) {
            std::cout << "==";
        }
        std::cout << std::endl;
        len /= 2;
    }
}

void GaussPyramid_mpi::GaussFilter(int theLayer) {//采用双边滤波
    float len = length;
    int t = theLayer;

    while (theLayer != 0) {
        theLayer--;
        len /= 2;
    }
    theLayer = t;
    int MyLen = len;
    len = (len - 1) / 2;
    int i;
    float sig = sigma;
    float* fillter = new float[length];

    for (i = 0; i < S + 3; ++i) {
        {
            sig = sigma / (i + 1);
            for (int j = 0; j < MyLen; ++j) {
                fillter[j] = exp(-(j - len) * (j - len) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
            }
        }
        for (int j = 0; j < MyLen; ++j) {
            for (int k = 0; k < MyLen; ++k) {
                GaussPy[theLayer][i][j][k] *= fillter[k];
            }
        }
        for (int j = 0; j < MyLen; ++j) {
            for (int k = 0; k < MyLen; ++k) {
                GaussPy[theLayer][i][k][j] *= fillter[k];
            }
        }
    }
    delete[]fillter;
}

void GaussPyramid_mpi::GenerateDoG() {
    int len = length;
    for (int i = 0; i < layer; ++i) {
        GaussFilter(i);
        int j = 0;
        for (j = 0; j < S + 2; j++) {
            for (int k = 0; k < len; ++k) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[i][j][k][l] -= GaussPy[i][j + 1][k][l];
                }
            }
        }
        len /= 2;
    }
}

GaussPyramid_mpi::~GaussPyramid_mpi() {
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



//void GaussPyramid_mpi::GenerateDoG_omp_dynamic() {
//    for (int i = 0; i < S; i++) {
//        int len = length;
//        int k = i;
//        int j = 0;
//        float* fil = new float[len];
//        while (len) {
//
//            float l = float(len - 1) / 2.0;
//            float sig = sigma / (i + 1);
//            {
//                for (int jj = 0; jj < len; ++jj) {
//                    fil[jj] = exp(-(jj - l) * (jj - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
//                }
//                for (int m = 0; m < len; ++m) {
//                    for (int n = 0; n < len; ++n) {
//                        GaussPy[j][i][m][n] *= fil[n];
//                    }
//                }
//                for (int m = 0; m < len; ++m) {
//                    for (int n = 0; n < len; ++n) {
//                        GaussPy[j][i][m][n] *= fil[m];
//                    }
//                }
//
//                {
//                    len /= 2;
//                    j++;
//
//                }
//            }
//
//        }
//        delete[]fil;
//    }
//    int len;
//    for (int MyLayer = 0; MyLayer < layer; MyLayer++) {
//        {
//            len = length;
//            int k = MyLayer;
//            while (k) {
//                k--;
//                len /= 2;
//            }
//        }
//        for (int i = 0; i < S - 1; ++i) {
//            for (int j = 0; j < len; ++j) {
//                for (int l = 0; l < len; ++l) {
//                    GaussPy[MyLayer][i][j][l] -= GaussPy[MyLayer][i + 1][j][l];
//                }
//            }
//        }
//    }
//}



void GaussPyramid_mpi::GenerateDoG_mpi(int argc,char** argv) {
    int i,len=length,j=0,jj,m,n,k,MyLayer,i1;
    float l,sig;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &i);
    if (i < S+3) {
        while (len) {
            l = float(len - 1) / 2.0;
            sig = sigma / (i + 1);
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
                    MPI_Send(GaussPy[j][i][m],len,MPI_INT,S+3, i*len+m,MPI_COMM_WORLD);
                }
            }
            j++;
            len /= 2;
        }
    }
    else if (i == S + 3) {
        len=length;
        m=0;
        while(len){
            for (j = 0; j < S+3; j++){
                for (i1 = 0; i1 < len; ++i1) {
                    MPI_Recv(GaussPy[m][j][i1],len, MPI_INT,j,j*len+i1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                }
            }
            len/=2;
            m++;
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
    MPI_Finalize();
}

void GaussPyramid_mpi::GenerateDoG_mpi_normal() {

}

#endif //SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_MPI_H
