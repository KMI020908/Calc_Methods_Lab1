#ifndef METHODS_H
#define METHODS_H

#include<vector>
#include<cmath>
#include<limits>
#include"ioData.h"
#include<iostream>

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = std::numeric_limits<Type>::epsilon());

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = std::numeric_limits<Type>::epsilon());

template<typename Type>
void findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = std::numeric_limits<Type>::epsilon());

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution); // Найти невязку

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &A); // Число обусловленности с октаэдоической метрикой

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &A); // Число обусловленности с кубической метрикой

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix); // Обратная матрица

#endif