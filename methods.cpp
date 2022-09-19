#include "methods.h"

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution,
const std::string &OUT_FILE_PATH, Type disturbance){
    size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    solution.resize(dimMatrix);
    if (disturbance >= std::numeric_limits<Type>::epsilon()){
        for (size_t i = 0; i < dimMatrix; i++){
            rCoefs[i] += disturbance; 
        }
    }
    const std::vector<std::vector<Type>> A = lCoefs; // Левая матрица до преобразований
    const std::vector<Type> b = rCoefs; // Правая матрица после преобразований
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
            if (OUT_FILE_PATH != "empty"){
                writeData(solution, OUT_FILE_PATH, NO_SOLUTION);
                Type cond_1 = 0;
                Type cond_inf = 0;
                writeConds(cond_1, cond_inf, OUT_FILE_PATH ,NO_SOLUTION);
            }
            return NO_SOLUTION;
        } 
    }
    // Обратный ход Гаусса
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
        if(OUT_FILE_PATH != "empty"){
            writeData(solution, OUT_FILE_PATH);
            Type cond_1 = findCond_1(A);
            Type cond_inf = findCond_inf(A);
            writeConds(cond_1, cond_inf, OUT_FILE_PATH);
            writeDiscrepancy(discrepancyVector, discrepancy, OUT_FILE_PATH);
        }
    }
    else{
        // Добавление возмущения
        if(OUT_FILE_PATH != "empty")
            addDisturbance(solution, OUT_FILE_PATH, disturbance);
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, 
const std::string &OUT_FILE_PATH, Type disturbance){
    size_t dimMatrix = lCoefs.size(); // Размерность СЛАУ
    solution.resize(dimMatrix);
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
                if (OUT_FILE_PATH != "empty"){
                    writeData(solution, OUT_FILE_PATH, NO_SOLUTION);
                    Type cond_1 = 0;
                    Type cond_inf = 0;
                    writeConds(cond_1, cond_inf, OUT_FILE_PATH ,NO_SOLUTION);
                }
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
     // Обратный ход Гаусса
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
        if (OUT_FILE_PATH != "empty"){
            writeData(solution, OUT_FILE_PATH);
            Type cond_1 = findCond_1(A);
            Type cond_inf = findCond_inf(A);
            writeConds(cond_1, cond_inf, OUT_FILE_PATH);
            writeDiscrepancy(discrepancyVector, discrepancy, OUT_FILE_PATH);
            writeQRMatrix(Q, lCoefs, OUT_FILE_PATH); 
        }
    }
    else{
        if (OUT_FILE_PATH != "empty"){
            addDisturbance(solution, OUT_FILE_PATH, disturbance);
            writeQRMatrix(Q, lCoefs, OUT_FILE_PATH);
        }
    }
    return HAS_SOLUTION;
}

// Нахождение числа обусловности для октаэдрической нормы
template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &A){
    size_t dimMatrix = A.size();
    std::vector<std::vector<Type>> B(dimMatrix); //Обратная к A матрица
    std::vector<std::vector<Type>> E(dimMatrix);
    for (size_t i = 0; i < dimMatrix; i++){
        E[i].resize(dimMatrix);
        B[i].resize(dimMatrix);
    }
    for (size_t i = 0; i < dimMatrix; i++){
        for (size_t j = 0; j < dimMatrix; j++){
            if (i == j)
                E[i][j] = 1.0;
            else
                E[i][j] = 0.0;   
        }
    }
    std::vector<Type> solution(dimMatrix);
    std::vector<std::vector<Type>> tempMatrix(dimMatrix);
    for (size_t i = 0; i < dimMatrix; i++){
        tempMatrix = A;
        gaussMethod<Type>(tempMatrix, E[i], solution);
        for (size_t j = 0; j < dimMatrix; j++)
            B[j][i] = solution[j];
    }
    Type infNormOfA = 0.0;
    Type infNormOfB = 0.0;
    // Норма A
    for (size_t j = 0; j < dimMatrix; j++){
        Type sum = 0.0;
        for (size_t i = 0; i < dimMatrix; i++){
            sum += abs(A[i][j]);
        }
        if (sum > infNormOfA)
            infNormOfA = sum;
    }
    // Норма обратной к A матрицы
    for (size_t j = 0; j < dimMatrix; j++){
        Type sum = 0.0;
        for (size_t i = 0; i < dimMatrix; i++){
            sum += abs(B[i][j]);
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
    size_t dimMatrix = A.size();
    std::vector<std::vector<Type>> B(dimMatrix); //Обратная к A матрица
    std::vector<std::vector<Type>> E(dimMatrix);
    for (size_t i = 0; i < dimMatrix; i++){
        E[i].resize(dimMatrix);
        B[i].resize(dimMatrix);
    }
    for (size_t i = 0; i < dimMatrix; i++){
        for (size_t j = 0; j < dimMatrix; j++){
            if (i == j)
                E[i][j] = 1.0;
            else
                E[i][j] = 0.0;   
        }
    }
    std::vector<Type> solution(dimMatrix);
    std::vector<std::vector<Type>> tempMatrix(dimMatrix);
    for (size_t i = 0; i < dimMatrix; i++){
         tempMatrix = A;
        gaussMethod<Type>(tempMatrix, E[i], solution);
        for (size_t j = 0; j < dimMatrix; j++)
            B[j][i] = solution[j];
    }
    
    Type infNormOfA = 0.0;
    Type infNormOfB = 0.0;
    // Норма A
    for (size_t i = 0; i < dimMatrix; i++){
        Type sum = 0.0;
        for (size_t j = 0; j < dimMatrix; j++){
            sum += abs(A[i][j]);
        }
        if (sum > infNormOfA)
            infNormOfA = sum;
    }
    // Норма обратной к A
    for (size_t i = 0; i < dimMatrix; i++){
        Type sum = 0.0;
        for (size_t j = 0; j < dimMatrix; j++){
            sum += abs(B[i][j]);
        }
        if (sum > infNormOfB)
            infNormOfB = sum;
    }
    Type cond = infNormOfA * infNormOfB;
    return cond;
}
