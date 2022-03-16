cmake_minimum_required(VERSION 3.21)
project(RolloffModelLib)

set(CMAKE_CXX_STANDARD 14)

add_library(SingleModel SHARED SingleModel.cpp jni.h jni_md.h SingleModel.h)
add_library(TimeModel SHARED TimeModel.cpp jni.h jni_md.h TimeModel.h)


