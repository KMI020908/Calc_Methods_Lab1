#ifndef IO_DATA_H
#define IO_DATA_H

#include<string>
#include<fstream> 
#include<iomanip> 
#include<vector>
#include"PRNG.h"
#include"Flags.h"

template<typename Type>
FILE_FLAG generateRandomTest(size_t equationDim, Type minValue, Type maxValue, const std::string& FILE_PATH);

template<typename Type>
FILE_FLAG readData(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs,  const std::string& IN_FILE_PATH);

template<typename Type>
FILE_FLAG writeData(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, SOLUTION_FLAG FLAG = HAS_SOLUTION);

template<typename Type>
FILE_FLAG writeQRMatrix(const std::vector<std::vector<Type>> &Q, const std::vector<std::vector<Type>> &R, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeMatrixMultiply(const std::vector<std::vector<Type>> &B, const std::vector<std::vector<Type>> &A, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG writeResidual(Type residual, const std::string& OUT_FILE_PATH);

template<typename Type>
FILE_FLAG addPerturbation(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type perturbation, SOLUTION_FLAG FLAG = HAS_SOLUTION);

template<typename Type>
FILE_FLAG writeConds(Type cond_1, Type cond_inf, const std::string& OUT_FILE_PATH, SOLUTION_FLAG FLAG = HAS_SOLUTION);

#endif