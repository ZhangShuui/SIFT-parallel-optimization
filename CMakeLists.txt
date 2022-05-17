cmake_minimum_required(VERSION 3.20)
project(SIFT_Guass_normal)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_STANDARD 14)

add_executable(SIFT_Guass_normal main.cpp GuassDePyramid.h GaussDePyramid-pThread.h GaussDePyramid-NEON.h GaussDePyramid-NEONxPTHREAD.h GaussDePyramid-SSExPTHREAD.h GaussDePyramid-AVXxPTHREAD.h GaussDePyramid-AVX512xPTHREAD.h GaussDePyramid-OpenMP.h)
