#ifndef METHODS_H
#define METHODS_H

#include<vector>
#include<cmath>
#include<limits>
#include"ioData.h"
#include<iostream>

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution,
const std::string &OUT_FILE_PATH = "empty", Type disturbance = 0.0);

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
const std::string &OUT_FILE_PATH = "empty", Type disturbance = 0.0);

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &A);

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &A);


#endif