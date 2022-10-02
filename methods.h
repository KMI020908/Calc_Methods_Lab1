#ifndef METHODS_H
#define METHODS_H

#include<vector>
#include<cmath>
#include"ioData.h"
#include<iostream>

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-6);

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
Type accuracy = 1e-6);

template<typename Type>
void findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy = 1e-6);

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution); // Найти невязку

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix); // Октаэдрическая норма матрицы

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &A); // Число обусловленности с октаэдоической метрикой

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix); // Кубическая норма матрицы

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &A); // Число обусловленности с кубической метрикой

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector); // Октаэдрическая нормам вектора

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector); // Кубическая норма вектора

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix); // Обратная матрица

template<typename Type>
void trasposeMatrix(std::vector<std::vector<Type>> &matrix);
#endif