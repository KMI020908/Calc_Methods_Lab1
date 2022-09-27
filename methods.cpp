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
        if (std::abs(mainValue) < accuracy) // detA = 0
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
    std::size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    solution.resize(dimMatrix);
    for (std::size_t k = 0; k < dimMatrix; k++){
        for (std::size_t i = k + 1; i < dimMatrix; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                Type s = lCoefs[i][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                for (std::size_t j = k; j < dimMatrix; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c*lCoefs[k][j] + s*lCoefs[i][j];
                    lCoefs[i][j] = -s*temp + c*lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < std::numeric_limits<Type>::epsilon())
                        lCoefs[i][j] = 0;
                }
                Type temp = rCoefs[k];
                rCoefs[k] = c*rCoefs[k] + s*rCoefs[i];
                rCoefs[i] = -s*temp + c*rCoefs[i];
            }
        }
    }
    if (std::abs(lCoefs[dimMatrix - 1][dimMatrix - 1]) < accuracy){  // detA = 0
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
void findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t dimMatrix = lCoefs.size();
    Q.resize(dimMatrix);
    for (std::size_t i = 0; i < dimMatrix; i++){
        Q[i].resize(dimMatrix, 0);
    }
    for (std::size_t i = 0; i < dimMatrix; i++){
        Q[i][i] = 1;
    }
    for (std::size_t k = 0; k < dimMatrix; k++){
        for (std::size_t i = k + 1; i < dimMatrix; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                Type s = lCoefs[i][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                for (std::size_t j = 0; j < dimMatrix; j++){
                    Type temp = Q[k][j];
                    Q[k][j] = c*Q[k][j] + s*Q[i][j];
                    Q[i][j] = -s*temp + c*Q[i][j];
                    if (std::abs(Q[i][j]) < std::numeric_limits<Type>::epsilon())
                        Q[i][j] = 0;
                }
                for (std::size_t j = k; j < dimMatrix; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c*lCoefs[k][j] + s*lCoefs[i][j];
                    lCoefs[i][j] = -s*temp + c*lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < std::numeric_limits<Type>::epsilon())
                        lCoefs[i][j] = 0;
                }
            }
        }
    }
    trasposeMatrix(Q);
}

template<typename Type>
void trasposeMatrix(std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    if (rows != 0){
        for (std::size_t i = 0; i < rows; i++){
            for (std::size_t j = i + 1; j < cols; j++){
                Type temp = matrix[i][j];
                matrix[i][j] = matrix[j][i];
                matrix[j][i] = temp;
            }
        }
    }
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

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix){
    std::size_t dimMatrix = matrix.size();
    Type norm1OfMatrix = 0;
    for (std::size_t j = 0; j < dimMatrix; j++){
        Type sum = 0.0;
        for (std::size_t i = 0; i < dimMatrix; i++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > norm1OfMatrix)
            norm1OfMatrix = sum;
    }
    return norm1OfMatrix;
}

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &A){
    std::size_t dimMatrix = A.size();
    std::vector<std::vector<Type>> B; //Обратная к A матрица
    invertMatrix(A, B);
    Type norm1OfA = findMatrixNorm1(A);
    Type norm1OfB = findMatrixNorm1(B);
    Type cond = norm1OfA * norm1OfB;
    return cond;
}

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix){
    std::size_t dimMatrix = matrix.size();
    Type normInfOfMatrix = 0;
    for (std::size_t i = 0; i < dimMatrix; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < dimMatrix; j++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > normInfOfMatrix)
            normInfOfMatrix = sum;
    }
    return normInfOfMatrix;
}

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &A){
    std::size_t dimMatrix = A.size();
    std::vector<std::vector<Type>> B; //Обратная к A матрица
    invertMatrix(A, B);
    Type normInfOfA = findMatrixNormInf(A);
    Type normInfOfB = findMatrixNormInf(B);
    Type cond = normInfOfA * normInfOfB;
    return cond;
}

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix){
    std::size_t dimMatrix = inputMatrix.size();
    resMatrix.resize(dimMatrix);
    std::vector<std::vector<Type>> E(dimMatrix);
    for (std::size_t i = 0; i < dimMatrix; i++){
        E[i].resize(dimMatrix, 0);
        resMatrix[i].resize(dimMatrix);
    }
    for (std::size_t i = 0; i < dimMatrix; i++){
        E[i][i] = 1;
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

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector){
    Type sum = 0;
    for (std::size_t i = 0; i < vector.size(); i++)
        sum += std::abs(vector[i]);
    return sum;
}

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector){
    Type max = std::abs(vector[0]);
    for (std::size_t i = 1; i < vector.size(); i++)
        if (std::abs(vector[i]) > max)
            max = std::abs(vector[i]);
    return max;
}