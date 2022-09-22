#include "methods.h"

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    solution.resize(dimMatrix);
    for (std::size_t k = 0; k < dimMatrix; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        for (std::size_t i = k + 1; i < dimMatrix; i++){   // Частичный выбор главного элемента
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainValue != lCoefs[k][k]){ //Замена строк
            Type temp;
            for (std::size_t i = 0; i < dimMatrix; i++){
                temp = lCoefs[k][i];
                lCoefs[k][i] = lCoefs[mainRow][i];
                lCoefs[mainRow][i] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        for (std::size_t i = k + 1; i < dimMatrix; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < dimMatrix; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = dimMatrix - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < dimMatrix; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    solution.resize(dimMatrix);
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
        if (abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
            
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
     // Обратный ход Гаусса
    for (int i = dimMatrix - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < dimMatrix; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }   
    return HAS_SOLUTION;
}

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution){
    std::size_t dimMatrix = lCoefs.size();
    std::vector<Type> b1(dimMatrix);  // Правая часть после подстановки полученного решения
        for (std::size_t i = 0; i < dimMatrix; i++){
                Type sum = 0.0;
                for (std::size_t k = 0; k < dimMatrix; k++){
                    sum += lCoefs[i][k] * solution[k];
                }
                b1[i] = sum;
        }
        std::vector<Type> discrepancyVector(dimMatrix); // Вектор невязки
        for (std::size_t i = 0; i < dimMatrix; i++){
            discrepancyVector[i] = rCoefs[i] - b1[i];
        }
        Type discrepancy = 0.0; // Невязка
        for (std::size_t i = 0; i < dimMatrix; i++){
            discrepancy += discrepancyVector[i] * discrepancyVector[i];     
        }
    return std::sqrt(discrepancy);
}

// Нахождение числа обусловности для октаэдрической нормы
template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &A){
    std::size_t dimMatrix = A.size();
    std::vector<std::vector<Type>> B; //Обратная к A матрица
    invertMatrix(A, B);
    Type infNormOfA = 0.0;
    Type infNormOfB = 0.0;
    // Норма A
    for (std::size_t j = 0; j < dimMatrix; j++){
        Type sum = 0.0;
        for (std::size_t i = 0; i < dimMatrix; i++){
            sum += std::abs(A[i][j]);
        }
        if (sum > infNormOfA)
            infNormOfA = sum;
    }
    // Норма обратной к A матрицы
    for (std::size_t j = 0; j < dimMatrix; j++){
        Type sum = 0.0;
        for (std::size_t i = 0; i < dimMatrix; i++){
            sum += std::abs(B[i][j]);
        }
        if (sum > infNormOfB)
            infNormOfB = sum;
    }
    Type cond = infNormOfA * infNormOfB;
    return cond;
}

// Нахождение числа обусловности для кубической нормы
template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &A){
    std::size_t dimMatrix = A.size();
    std::vector<std::vector<Type>> B; //Обратная к A матрица
    invertMatrix(A, B);
    Type infNormOfA = 0.0;
    Type infNormOfB = 0.0;
    // Норма A
    for (std::size_t i = 0; i < dimMatrix; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < dimMatrix; j++){
            sum += std::abs(A[i][j]);
        }
        if (sum > infNormOfA)
            infNormOfA = sum;
    }
    // Норма обратной к A
    for (std::size_t i = 0; i < dimMatrix; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < dimMatrix; j++){
            sum += std::abs(B[i][j]);
        }
        if (sum > infNormOfB)
            infNormOfB = sum;
    }
    Type cond = infNormOfA * infNormOfB;
    return cond;
}

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix){
    std::size_t dimMatrix = inputMatrix.size();
    resMatrix.resize(dimMatrix);
    std::vector<std::vector<Type>> E(dimMatrix);
    for (std::size_t i = 0; i < dimMatrix; i++){
        E[i].resize(dimMatrix);
        resMatrix[i].resize(dimMatrix);
    }
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            if (i == j)
                E[i][j] = 1.0;
            else
                E[i][j] = 0.0;   
        }
    }
    std::vector<Type> solution(dimMatrix);
    std::vector<std::vector<Type>> tempMatrix(dimMatrix);
    SOLUTION_FLAG flag;
    for (std::size_t i = 0; i < dimMatrix; i++){
        tempMatrix = inputMatrix;
        flag = gaussMethod<Type>(tempMatrix, E[i], solution);
        if (flag == NO_SOLUTION){
            for (std::size_t i = 0; i < dimMatrix; i++)
                resMatrix[i].clear();
            resMatrix.clear();
            return NOT_INVERTIBLE;
        }
        for (std::size_t j = 0; j < dimMatrix; j++)
            resMatrix[j][i] = solution[j];
    }
    return IS_INVERTIBLE;
}