//
// Created by zsr on 2022/4/23.
//

#ifndef SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_PTHREAD_H
#define SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_PTHREAD_H
#define MAX_COUNT 7

//
// Created by zsr on 2022/4/17.
//
//const float sigma = 2.0;
//const float PI = 3.1414926;
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

timeval start, final;

class GaussPyramid_p;

struct ThreadParam_t {
    GaussPyramid_p *self;
    int theLayer;//记录层数
};//用在高斯金字塔同层间相减的参数


class GaussPyramid_p {
public:
    static int THREAD_COUNT;
    pthread_barrier_t barrier_Division;
    pthread_barrier_t barrier_Sub;
    int **data; //记录图像灰度数据
    GaussPyramid_p();

    GaussPyramid_p(int **img, int len/*,int wid*/, int S);

    float ****GaussPy;

    void GaussPyInit();

    void output();

    void GaussFilter(int theLayer);

    void GenerateDoG();

    void GenerateDoG_i();

    void Timing_ThreadCount();

    ~GaussPyramid_p();

    static void *thread_Filter_Sub(void *param);

    static void *thread_Filter_Sub_inner_layer(void *param);

protected:
    int length;
    //int width;  先尝试宽度相同的版本
    int S;  //提取图像特征后,需要进行对比的图片数
    int layer;
    float *filter;
};

int GaussPyramid_p::THREAD_COUNT = 7;

GaussPyramid_p::GaussPyramid_p() {
    data = nullptr;
}

GaussPyramid_p::GaussPyramid_p(int **img, int len, int S) {
    length = len;
    data = new int *[len];
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
    GaussPy = new float ***[layer];
    filter = new float[length];
    GaussPyInit();
}

//初始化高斯金字塔，尚未进行高斯滤波操作
void GaussPyramid_p::GaussPyInit() {
    int step = 1;
    for (int i = 0; i < layer; ++i) {
        GaussPy[i] = new float **[S + 3];
        for (int j = 0; j < S + 3; ++j) {
            GaussPy[i][j] = new float *[length / step];
            for (int k = 0; k < length / step; ++k) {
                GaussPy[i][j][k] = new float[length / step];
            }
        }
        step *= 2;
    }
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

void GaussPyramid_p::output() {
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

void GaussPyramid_p::GaussFilter(int theLayer) {//采用双边滤波
    float len = length;
    int t = theLayer;
    while (theLayer != 0) {
        theLayer--;
        len /= 2;
    }
    theLayer = t;
    int MyLen = len;
    len = (len - 1) / 2;

    for (int i = 0; i < S + 3; ++i) {
        float sig = sigma / (i + 1);
        for (int i = 0; i < MyLen; ++i) {
            filter[i] = exp(-(i - len) * (i - len) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
        }
        for (int j = 0; j < MyLen; ++j) {
            for (int k = 0; k < MyLen; ++k) {
                GaussPy[theLayer][i][j][k] *= filter[k];
            }
        }
        for (int j = 0; j < MyLen; ++j) {
            for (int k = 0; k < MyLen; ++k) {
                GaussPy[theLayer][i][k][j] *= filter[k];
            }
        }
    }

}

void GaussPyramid_p::GenerateDoG() {
    pthread_barrier_init(&barrier_Division, NULL, THREAD_COUNT);
    pthread_t handles[THREAD_COUNT];
    ThreadParam_t paras[THREAD_COUNT];
    for (int i = 0; i < THREAD_COUNT; ++i) {
        paras[i] = {this, i};
        pthread_create(handles + i, NULL, thread_Filter_Sub, paras + i);
    }
    for (int i = 0; i < THREAD_COUNT; ++i) {
        pthread_join(handles[i], nullptr);
    }
    pthread_barrier_destroy(&barrier_Division);
}

GaussPyramid_p::~GaussPyramid_p() {
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

void *GaussPyramid_p::thread_Filter_Sub(void *param) {
    ThreadParam_t *p = (ThreadParam_t *) param;
    GaussPyramid_p *self = p->self;
    int layer = self->layer;
    int length = self->length;
    float ****GaussPy = self->GaussPy;
    int S = self->S;
    int MyLayer = p->theLayer;
    for (int i = MyLayer; i < layer; i += THREAD_COUNT) {
        int len = length;
        int k = i;
        while (k--) {
            len /= 2;
        }
        float *fil = new float[len];
        float l = float(len - 1) / 2.0;
        for (int i = 0; i < S + 3; ++i) {
            float sig = sigma / (i + 1);
            for (int i = 0; i < len; ++i) {
                fil[i] = exp(-(i - l) * (i - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
            }
            for (int j = 0; j < len; ++j) {
                for (int k = 0; k < len; ++k) {
                    GaussPy[MyLayer][i][j][k] *= fil[k];
                }
            }
            for (int j = 0; j < len; ++j) {
                for (int k = 0; k < len; ++k) {
                    GaussPy[MyLayer][i][k][j] *= fil[k];
                }
            }
        }
        for (int i = 0; i < S + 2; ++i) {
            for (int j = 0; j < len; ++j) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[MyLayer][i][j][l] -= GaussPy[MyLayer][i + 1][j][l];
                }
            }
        }
        delete[]fil;
    }
    pthread_barrier_wait(&self->barrier_Division);
    pthread_exit(NULL);
}

void *GaussPyramid_p::thread_Filter_Sub_inner_layer(void *param) {
    ThreadParam_t *p = (ThreadParam_t *) param;
    GaussPyramid_p *self = p->self;
    int layer = self->layer;
    int length = self->length;
    float ****GaussPy = self->GaussPy;
    int S = self->S + 3;
    int MyS = p->theLayer;
    for (int i = MyS; i < S; i += THREAD_COUNT) {
        int len = length;
        int k = i;
        int j = 0;
        while (len) {
            float *fil = new float[len];
            float l = float(len - 1) / 2.0;
            float sig = sigma / (i + 1);
            for (int jj = 0; jj < len; ++jj) {
                fil[jj] = exp(-(jj - l) * (jj - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
            }
            for (int m = 0; m < len; ++m) {
                for (int n = 0; n < len; ++n) {
                    GaussPy[j][i][m][n] *= fil[n];
                }
            }
            for (int m = 0; m < len; ++m) {
                for (int n = 0; n < len; ++n) {
                    GaussPy[j][i][m][n] *= fil[m];
                }
            }
            len /= 2;
            j++;
            delete[]fil;
        }
    }
    pthread_barrier_wait(&self->barrier_Division);
    for (int MyLayer = MyS; MyLayer < layer; MyLayer += THREAD_COUNT) {
        int len = length;
        int k = MyLayer;
        while (k) {
            k--;
            len /= 2;
        }
        for (int i = 0; i < S-1; ++i) {
            for (int j = 0; j < len; ++j) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[MyLayer][i][j][l] -= GaussPy[MyLayer][i + 1][j][l];
                }
            }
        }
    }
    pthread_barrier_wait(&self->barrier_Sub);
    pthread_exit(NULL);
}

void GaussPyramid_p::Timing_ThreadCount() {
    float total = 0.0;
    int count = 0;
    for (int i = 1; i <= MAX_COUNT; ++i) {
        THREAD_COUNT = i;
        GaussPyInit();
        while (total < 10.0) {
            GaussPyInit();
            gettimeofday(&start, NULL);
            GenerateDoG();
            gettimeofday(&final, NULL);
            total += final.tv_sec - start.tv_sec;
            count++;
        }
        std::cout << i << "," << total / float(count) << std::endl;
    }
}

void GaussPyramid_p::GenerateDoG_i() {
    pthread_barrier_init(&barrier_Division, NULL, THREAD_COUNT);
    pthread_barrier_init(&barrier_Sub, NULL, THREAD_COUNT);
    pthread_t handles[THREAD_COUNT];
    ThreadParam_t paras[THREAD_COUNT];
    for (int i = 0; i < THREAD_COUNT; ++i) {
        paras[i] = {this, i};
        pthread_create(handles + i, NULL, thread_Filter_Sub_inner_layer, paras + i);
    }
    for (int i = 0; i < THREAD_COUNT; ++i) {
        pthread_join(handles[i], nullptr);
    }
    pthread_barrier_destroy(&barrier_Division);
    pthread_barrier_destroy(&barrier_Sub);
}


#endif //SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_PTHREAD_H
