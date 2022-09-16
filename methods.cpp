#include<vector>
#include<cmath>
#include<limits>
#include"ioData.h"
#include<iostream>


template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, const std::string &OUT_FILE_PATH, Type disturbance = 0.0){
    size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    if (disturbance >= std::numeric_limits<Type>::epsilon()){
        for (size_t i = 0; i < dimMatrix; i++){
            rCoefs[i] += disturbance; 
        }
    }
    std::vector<std::vector<Type>> A = lCoefs; // Левая матрица до преобразований
    std::vector<Type> b = rCoefs; // Правая матрица после преобразований
    for (size_t k = 0; k < dimMatrix; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        size_t mainRow = k; // Строка главного элемента
        for (size_t i = k + 1; i < dimMatrix; i++){   // Частичный выбор главного элемента
            if (abs(lCoefs[i][k]) > abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainValue != lCoefs[k][k]){ //Замена строк
            Type temp;
            for (size_t i = 0; i < dimMatrix; i++){
                temp = lCoefs[k][i];
                lCoefs[k][i] = lCoefs[mainRow][i];
                lCoefs[mainRow][i] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        for (size_t i = k + 1; i < dimMatrix; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (size_t j = k;  j < dimMatrix; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (abs(mainValue) < std::numeric_limits<Type>::epsilon()){ // detA = 0
            std::vector<Type> solution(0, 0);
            writeData(solution, OUT_FILE_PATH, NO_SOLUTION);
            return NO_SOLUTION;
        } 
    }
    std::vector<Type> solution(dimMatrix); // Обратный ход Гаусса
    for (int i = dimMatrix - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (size_t j = i + 1; j < dimMatrix; j++){
            sum += lCoefs[i][j] * solution[j]; 
        }
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    if (disturbance < std::numeric_limits<Type>::epsilon()){ // Ищем невязку для системы без возмущения
         std::vector<Type> b1(dimMatrix);  // Правая часть после подстановки полученного решения
        for (size_t i = 0; i < dimMatrix; i++){
                Type sum = 0.0;
                for (size_t k = 0; k < dimMatrix; k++){
                    sum += A[i][k] * solution[k];
                }
                b1[i] = sum;
        }
        std::vector<Type> discrepancyVector(dimMatrix); // Вектор невязки
        for (size_t i = 0; i < dimMatrix; i++){
            discrepancyVector[i] = b[i] - b1[i];
        }
        Type discrepancy = 0.0; // Невязка
        for (size_t i = 0; i < dimMatrix; i++){
            discrepancy += discrepancyVector[i] * discrepancyVector[i];     
        }
        discrepancy = sqrt(discrepancy);
        // Вывод данных
        writeData(solution, OUT_FILE_PATH);
        writeDiscrepancy(discrepancyVector, discrepancy, OUT_FILE_PATH); 
    }
    else{
        // Добавление возмущения
        addDisturbance(solution, OUT_FILE_PATH, disturbance);
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, const std::string &OUT_FILE_PATH, Type disturbance = 0.0){
    size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    if (disturbance >= std::numeric_limits<Type>::epsilon()){
            for (size_t i = 0; i < dimMatrix; i++){
                rCoefs[i] += disturbance; 
            }
        }
    std::vector<std::vector<Type>> A = lCoefs;
    std::vector<Type> b = rCoefs;
    for (size_t k = 0; k < dimMatrix; k++){
        Type mainValue = lCoefs[k][k];    // Главный элемент
        size_t mainRow = k; // Строка главного элемента
        for (size_t i = k + 1; i < dimMatrix; i++){   // Частичный выбор главного элемента
            if (abs(lCoefs[i][k]) > abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainValue != lCoefs[k][k]){ //Замена строк
            Type temp;
            for (size_t i = 0; i < dimMatrix; i++){
                temp = lCoefs[k][i];
                lCoefs[k][i] = lCoefs[mainRow][i];
                lCoefs[mainRow][i] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        if (abs(mainValue) < std::numeric_limits<Type>::epsilon()){ // detA = 0
                std::vector<Type> solution(0, 0);
                writeData(solution, OUT_FILE_PATH, NO_SOLUTION);
                return NO_SOLUTION;
            } 
        for (size_t i = k + 1; i < dimMatrix; i++){
            Type c = lCoefs[k][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
            Type s = lCoefs[i][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
            for (size_t j = k; j < dimMatrix; j++){
                Type temp = lCoefs[k][j];
                lCoefs[k][j] = c*lCoefs[k][j] + s*lCoefs[i][j];
                lCoefs[i][j] = -s*temp + c*lCoefs[i][j];
                if (abs(lCoefs[i][j]) < std::numeric_limits<Type>::epsilon())
                    lCoefs[i][j] = 0;
            }
            Type temp = rCoefs[k];
            rCoefs[k] = c*rCoefs[k] + s*rCoefs[i];
            rCoefs[i] = -s*temp + c*rCoefs[i];
        }
    }
    // Нахождение Q матрицы
    std::vector<std::vector<Type>> R_rev(dimMatrix); // обратная к R матрица
    for (size_t i = 0; i < dimMatrix; i++){
        R_rev[i].resize(dimMatrix);
    }
    for (size_t i = 0; i < dimMatrix; i++){
        R_rev[i][i] = 1/lCoefs[i][i];
    }
    for (size_t i = 0; i < dimMatrix; i++){
        for (size_t j = i + 1; j < dimMatrix; j++){
            Type sum = 0.0;
            for (size_t k = 0; k < j; k++){
                sum += R_rev[i][k] * lCoefs[k][j];
            }
            R_rev[i][j] = -R_rev[j][j]*sum;
        } 
    }
    std::vector<std::vector<Type>> Q(dimMatrix); // Искомая матрица Q
    for (size_t i = 0; i < dimMatrix; i++){
        Q[i].resize(dimMatrix);
    }
    for (size_t i = 0; i < dimMatrix; i++){
        for (size_t j = 0; j < dimMatrix; j++){
            Type sum = 0.0;
            for (size_t k = 0; k < dimMatrix; k++){
                sum += A[i][k] * R_rev[k][j];
            }
            Q[i][j] = sum;
        } 
    }
    std::vector<Type> solution(dimMatrix); // Обратный ход Гаусса
    for (int i = dimMatrix - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (size_t j = i + 1; j < dimMatrix; j++){
            sum += lCoefs[i][j] * solution[j]; 
            }
            solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
        }

    if (disturbance < std::numeric_limits<Type>::epsilon()){ // Ищем невязку для системы без возмущения
        std::vector<Type> b1(dimMatrix); // Правая часть после подстановки полученного решения
        for (size_t i = 0; i < dimMatrix; i++){
                Type sum = 0.0;
                for (size_t k = 0; k < dimMatrix; k++){
                    sum += A[i][k] * solution[k];
                }
                b1[i] = sum;
        }
        std::vector<Type> discrepancyVector(dimMatrix); // Вектор невязки
        for (size_t i = 0; i < dimMatrix; i++){
            discrepancyVector[i] = b[i] - b1[i];
        }
        Type discrepancy = 0.0;
        for (size_t i = 0; i < dimMatrix; i++){
            discrepancy += discrepancyVector[i] * discrepancyVector[i];     
        }
        discrepancy = sqrt(discrepancy); // Невязка
        // Вывод данных
        writeData(solution, OUT_FILE_PATH);
        writeDiscrepancy(discrepancyVector, discrepancy, OUT_FILE_PATH);
        writeQRMatrix(Q, lCoefs, OUT_FILE_PATH); 
    }
    else{
        addDisturbance(solution, OUT_FILE_PATH, disturbance);
        writeQRMatrix(Q, lCoefs, OUT_FILE_PATH);
    }
    return HAS_SOLUTION;
}