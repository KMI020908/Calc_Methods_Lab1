#include<vector>
#include<cmath>
#include<limits>
#include"ioData.h"


template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, const std::string OUT_FILE_PATH){
    size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
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
            sum = sum + lCoefs[i][j] * solution[j]; 
        }
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    writeData(solution, OUT_FILE_PATH);
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, const std::string OUT_FILE_PATH){
size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
for (size_t k = 0; k < dimMatrix - 1; k++){
    for (size_t i = k + 1; i < dimMatrix; i++){
        Type c = lCoefs[k][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
        Type s = lCoefs[i][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
        for (size_t j = k; j < dimMatrix; j++){
            Type temp = lCoefs[k][j];
            lCoefs[k][j] = c*lCoefs[k][j] + s*lCoefs[i][j];
            lCoefs[i][j] = -s*temp + c*lCoefs[i][j];
            if (abs(lCoefs[i][j]) - 1e-15 < 0)
                lCoefs[i][j] = 0;
        }
        rCoefs[k] = c*rCoefs[k] + s*rCoefs[i];
        rCoefs[i] = -s*rCoefs[k] + c*rCoefs[i];
    }
}
std::vector<Type> solution(dimMatrix); // Обратный ход Гаусса
for (int i = dimMatrix - 1; i >= 0 ; i--){
    Type sum = 0.0;
    for (size_t j = i + 1; j < dimMatrix; j++){
         sum = sum + lCoefs[i][j] * solution[j]; 
        }
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    writeData(solution, OUT_FILE_PATH);
    return HAS_SOLUTION;

return HAS_SOLUTION;
}