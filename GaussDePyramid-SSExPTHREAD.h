//
// Created by zsr on 2022/4/26.
//
#ifndef SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_SSEXPTHREAD_H
#define SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_SSEXPTHREAD_H
#define THREAD_COUNT_st 7
#include <iostream>
#include <math.h>
#include <immintrin.h>
#include <pthread.h>
class GaussPyramid_sxp;
struct ThreadParam_st{
    GaussPyramid_sxp* self;
    int theLayer;//记录层数
};//用在高斯金字塔同层间相减的参数


class GaussPyramid_sxp{
public:
    pthread_barrier_t barrier_Division;
    int** data; //记录图像灰度数据
    GaussPyramid_sxp();
    GaussPyramid_sxp(int** img, int len/*,int wid*/, int S);
    float**** GaussPy;
    void GaussPyInit();
    void output();
    void GaussFilter(int theLayer);
    void GenerateDoG();
    ~GaussPyramid_sxp();
    static void* thread_Filter_Sub(void *param);
protected:
    int length;
    //int width;  先尝试宽度相同的版本
    int S;  //提取图像特征后,需要进行对比的图片数
    int layer;
    float* filter;
};

GaussPyramid_sxp::GaussPyramid_sxp() {
    data= nullptr;
}

GaussPyramid_sxp::GaussPyramid_sxp(int **img, int len, int S) {
    length=len;
    data= new int*[len];
    for (int i = 0; i < len; ++i) {
        data[i]=new int [len];
    }
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            data[i][j]=img[i][j];
        }
    }
    this->S=S;
    int x=0;
    while (len){
        x++;
        len/=2;
    }
    layer=x;
    GaussPy=new float***[layer];
    filter=new float[length];
    GaussPyInit();
}
//初始化高斯金字塔，尚未进行高斯滤波操作
void GaussPyramid_sxp::GaussPyInit() {
    int step=1;
    for (int i = 0; i < layer; ++i) {
        GaussPy[i]=new float**[S+3];
        for (int j = 0; j < S + 3; ++j) {
            GaussPy[i][j]=new float*[length/step];
            for (int k = 0; k < length/step; ++k) {
                GaussPy[i][j][k]=new float [length/step];
            }
        }
        step*=2;
    }
    int len=length;
    step=1;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            for (int k = 0; k < len; ++k) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[i][j][k][l]=data[k*step][l*step];
                }
            }
        }
        len/=2;
        step*=2;
    }
}

void GaussPyramid_sxp::output() {
    int len=length;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < len; ++j) {
            for (int k = 0; k < len; ++k) {
                std::cout<<GaussPy[i][0][j][k]<<" ";
            }
            std::cout<<std::endl;
        }
        for (int k = 0; k < len; ++k) {
            std::cout<<"==";
        }
        std::cout<<std::endl;
        len/=2;
    }
}

void GaussPyramid_sxp::GaussFilter(int theLayer) {//采用双边滤波
    float len=length;
    int t=theLayer;
    while (theLayer!=0){
        theLayer--;
        len/=2;
    }
    theLayer=t;
    int MyLen=len;
    len=(len-1)/2;

    for (int i = 0; i < S + 3; ++i) {
        float sig=sigma/(i+1);
        for (int i = 0; i < MyLen; ++i) {
            filter[i] = exp(-(i-len)*(i-len)/(2*sig*sig))/(sig*sqrt(2*PI));
        }
        for (int j = 0; j < MyLen; ++j) {
            for (int k = 0; k < MyLen; ++k) {
                GaussPy[theLayer][i][j][k]*=filter[k];
            }
        }
        for (int j = 0; j < MyLen; ++j) {
            for (int k = 0; k < MyLen; ++k) {
                GaussPy[theLayer][i][k][j]*=filter[k];
            }
        }
    }

}

void GaussPyramid_sxp::GenerateDoG() {
    pthread_barrier_init(&barrier_Division,NULL,THREAD_COUNT_st);
    pthread_t handles[THREAD_COUNT_st];
    ThreadParam_st paras[THREAD_COUNT_st];
    for (int i = 0; i < THREAD_COUNT_st; ++i) {
        paras[i]={this,i};
        pthread_create(handles+i,NULL, thread_Filter_Sub,paras+i);
    }
    for (int i = 0; i < THREAD_COUNT_st; ++i) {
        pthread_join(handles[i], nullptr);
    }
    pthread_barrier_destroy(&barrier_Division);
}

GaussPyramid_sxp::~GaussPyramid_sxp() {
    int step=1;
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            for (int k = 0; k < length/step; ++k) {
                delete []GaussPy[i][j][k];
            }
        }
        step*=2;
    }
    for (int i = 0; i < layer; ++i) {
        for (int j = 0; j < S + 3; ++j) {
            delete []GaussPy[i][j];
        }
    }
    for (int i = 0; i < layer; ++i) {
        delete []GaussPy[i];
    }
    delete []GaussPy;
}

void *GaussPyramid_sxp::thread_Filter_Sub(void *param) {
    ThreadParam_st * p=(ThreadParam_st *) param;
    GaussPyramid_sxp* self=p->self;
    int layer=self->layer;
    int length=self->length;
    float**** GaussPy=self->GaussPy;
    int S=self->S;
    int MyLayer=p->theLayer;
    float* fil =new float[length];
    for (int i = MyLayer; i < layer; i+=THREAD_COUNT_st) {
        int len=length;
        int k=i;
        while (k--){
            len/=2;
        }
        if (len<=2){
            float l=float (len-1)/2.0;
            for (int mm = 0; mm < S + 3; ++mm) {
                float sig=sigma/float (mm + 1);
                for (int kk = 0; kk < len; ++kk) {
                    fil[kk] = exp(-(float (kk) - l) * (float (kk) - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
                }
                for (int j = 0; j < len; ++j) {
                    for (int kk = 0; kk < len; ++kk) {
                        GaussPy[i][mm][j][kk]*=fil[kk];
                    }
                }
                for (int j = 0; j < len; ++j) {
                    for (int kk = 0; kk < len; ++kk) {
                        GaussPy[i][mm][kk][j]*=fil[kk];
                    }
                }
            }
            for (int j = 0; j < S + 2; ++j) {
                for (int kk = 0; kk < len; ++kk) {
                    for (int ll = 0; ll < len; ++ll) {
                        GaussPy[i][j][kk][ll]-=GaussPy[i][j + 1][kk][ll];
                    }
                }
            }
        } else {
            __m128 vf,va;
            float l=float (len-1)/2.0;
            for (int m = 0; m < S + 3; ++m) {
                float sig=sigma/(m+1);
                for (int kk= 0; kk < len; ++kk) {
                    fil[kk] = exp(-(kk - l) * (kk - l) / (2 * sig * sig)) / (sig * sqrt(2 * PI));
                }
                for (int j = 0; j < len; j+=4) {
                    vf = _mm_loadu_ps(fil+j);
                    for (int kk = 0; kk < len; ++kk) {
                        va = _mm_loadu_ps(GaussPy[i][m][kk] + j);
                        va = _mm_mul_ps(va,vf);
                         _mm_store_ps(GaussPy[i][m][kk] + j, va);
                    }
                }
                for (int j = 0; j < len; j++) {
                    vf = _mm_set1_ps(fil[j]);
                    for (int kk = 0; kk < len; kk+=4) {
                        va = _mm_loadu_ps(GaussPy[i][m][j] + kk);
                        va = _mm_mul_ps(va,vf);
                         _mm_store_ps(GaussPy[i][m][j] + kk, va);
                    }
                }
            }
            for (int j = 0; j < S + 2; ++j) {
                for (int kk = 0; kk < len; ++kk) {
                    for (int ll = 0; ll < len; ll+=4) {
                        va = _mm_loadu_ps(GaussPy[i][j][kk] + ll);
                        vf = _mm_loadu_ps(GaussPy[i][j+1][kk] + ll);
                        va = _mm_sub_ps(va,vf);
                         _mm_store_ps(GaussPy[i][j][kk] + ll, va);
                    }
                }
            }
        }
    }
    delete []fil;
    pthread_barrier_wait(&self->barrier_Division);
    self= nullptr;
    GaussPy= nullptr;
    p= nullptr;
    pthread_exit(NULL);
}









#endif //SIFT_GUASS_NORMAL_GAUSSDEPYRAMID_SSEXPTHREAD_H
