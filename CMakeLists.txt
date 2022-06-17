cmake_minimum_required(VERSION 3.20)
project(MCCPS_luban)
set(CMAKE_CXX_STANDARD 14)

# openMP 配置
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${OpenMP_C_FLAGS}")
endif()


add_executable(MCCPS_luban luban_mc_nvt.cpp  luban_mc_nvt.h)