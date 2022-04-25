//
// Created by zsr on 2022/4/17.
//

#ifndef SIFT_GUASS_NORMAL_GUASSDEPYRAMID_H
#define SIFT_GUASS_NORMAL_GUASSDEPYRAMID_H
const float sigma = 2.0;
const float PI = 3.1414926;
#include <iostream>
#include <math.h>
class GaussPyramid{
public:
    int** data; //记录图像灰度数据
    GaussPyramid();
    GaussPyramid(int** img, int len/*,int wid*/, int S);
    float**** GaussPy;
    void GaussPyInit();
    void output();
    void GaussFilter(int theLayer);
    void GenerateDoG();
    ~GaussPyramid();
protected:
    int length;
    //int width;  先尝试宽度相同的版本
    int S;  //提取图像特征后,需要进行对比的图片数
    int layer;
    float* filter;
};

GaussPyramid::GaussPyramid() {
    data= nullptr;
}

GaussPyramid::GaussPyramid(int **img, int len, int S) {
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
void GaussPyramid::GaussPyInit() {
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

void GaussPyramid::output() {
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

void GaussPyramid::GaussFilter(int theLayer) {//采用双边滤波
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

void GaussPyramid::GenerateDoG() {
    int len=length;
    for (int i = 0; i < layer; ++i) {
        GaussFilter(i);
        for (int j = 0; j < S + 2; ++j) {
            for (int k = 0; k < len; ++k) {
                for (int l = 0; l < len; ++l) {
                    GaussPy[i][j][k][l]-=GaussPy[i][j+1][k][l];
                }
            }
        }
        len/=2;
    }
}

GaussPyramid::~GaussPyramid() {
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


#endif //SIFT_GUASS_NORMAL_GUASSDEPYRAMID_H
