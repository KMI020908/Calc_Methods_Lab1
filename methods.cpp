#include "methods.h"

template<typename Type>
SOLUTION_FLAG gaussMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    for (std::size_t k = 0; k < rows; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        for (std::size_t i = k + 1; i < rows; i++){   // Частичный выбор главного элемента
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainRow != k){ //Замена строк
            Type temp;
            for (std::size_t j = 0; j < cols; j++){
                temp = lCoefs[k][j];
                lCoefs[k][j] = lCoefs[mainRow][j];
                lCoefs[mainRow][j] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        for (std::size_t i = k + 1; i < rows; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < cols; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (std::abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG gaussMethodFull(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    std::vector<std::size_t> switchCols; // Вектор, хранящий индексы перемещенных столбцов исходной матрицы
    for (std::size_t k = 0; k < rows; k++){ // Прямой ход Гаусса
        Type mainValue = lCoefs[k][k];    // Главный элемент
        std::size_t mainRow = k; // Строка главного элемента
        std::size_t mainСol = k; // Столбец главного элемента
        // Полный выбор
        // Поиск главного элмента в k-ом столбце  
        for (std::size_t i = k + 1; i < rows; i++){ 
            if (std::abs(lCoefs[i][k]) > std::abs(mainValue)){
                mainValue = lCoefs[i][k]; 
                mainRow = i;
            }
        }
        if (mainRow != k){ //Замена строк
            Type temp;
            for (std::size_t j = 0; j < cols; j++){
                temp = lCoefs[k][j];
                lCoefs[k][j] = lCoefs[mainRow][j];
                lCoefs[mainRow][j] = temp;
            }
            temp =rCoefs[k];
            rCoefs[k] = rCoefs[mainRow];
            rCoefs[mainRow] = temp;
        }
        // Поиск главного элмента в k-ой строке 
        for (std::size_t j = k + 1; j < cols; j++){ 
            if (std::abs(lCoefs[k][j]) > std::abs(mainValue)){
                mainValue = lCoefs[k][j]; 
                mainСol = j;
            }
        }
        //Замена столбцов
        if (mainСol != k){ 
            Type temp;
            for (std::size_t i = 0; i < rows; i++){
                temp = lCoefs[i][k];
                lCoefs[i][k] = lCoefs[i][mainСol];
                lCoefs[i][mainСol] = temp;
            }
            switchCols.push_back(k);
            switchCols.push_back(mainСol);
        }

        // Прямой ход Гаусса 
        for (std::size_t i = k + 1; i < rows; i++){ 
            Type C = lCoefs[i][k]/lCoefs[k][k];
            rCoefs[i] = rCoefs[i] - C*rCoefs[k];
            for (std::size_t j = k;  j < cols; j++){
                lCoefs[i][j] = lCoefs[i][j] - C*lCoefs[k][j];
            }  
        }
        if (std::abs(mainValue) < accuracy) // detA = 0
            return NO_SOLUTION;
    }
    // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }
    // Обратная перестановка
    for (int i = switchCols.size() - 2; i >= 0; i -= 2){
        Type temp = solution[switchCols[i]];
        solution[switchCols[i]] = solution[switchCols[i + 1]];
        solution[switchCols[i + 1]] = temp;

    }
    return HAS_SOLUTION;
}

template<typename Type>
SOLUTION_FLAG qrMethod(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, std::vector<Type> &solution, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    solution.resize(rows); // Искомое решение
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NO_SOLUTION;
    if (rows != cols){
        solution.resize(0);
        return NO_SOLUTION;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                Type s = lCoefs[i][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                for (std::size_t j = k; j < cols; j++){
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
    if (std::abs(lCoefs[rows - 1][rows - 1]) < accuracy){  // detA = 0
        return NO_SOLUTION;
    }
    
     // Обратный ход Гаусса
    for (int i = rows - 1; i >= 0 ; i--){
        Type sum = 0.0;
        for (std::size_t j = i + 1; j < cols; j++)
            sum += lCoefs[i][j] * solution[j]; 
        solution[i] = (rCoefs[i] - sum)/lCoefs[i][i];
    }   
    return HAS_SOLUTION;
}

template<typename Type>
QUADRATIC_FLAG findQMatrix(std::vector<std::vector<Type>> &lCoefs, std::vector<std::vector<Type>> &Q, Type accuracy){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return IS_QUADRATIC;
    if (rows != cols)
        return NOT_QUADRATIC;
    Q.resize(rows);
    for (std::size_t i = 0; i < rows; i++){
        Q[i].resize(cols, 0);
    }
    for (std::size_t i = 0; i < rows; i++){
        Q[i][i] = 1;
    }
    for (std::size_t k = 0; k < rows; k++){
        for (std::size_t i = k + 1; i < rows; i++){
            if (std::abs(lCoefs[i][k]) >= accuracy){
                Type c = lCoefs[k][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                Type s = lCoefs[i][k]/sqrt(lCoefs[k][k]*lCoefs[k][k] + lCoefs[i][k]*lCoefs[i][k]);
                for (std::size_t j = 0; j < cols; j++){
                    Type temp = Q[k][j];
                    Q[k][j] = c*Q[k][j] + s*Q[i][j];
                    Q[i][j] = -s*temp + c*Q[i][j];
                    if (std::abs(Q[i][j]) < accuracy)
                        Q[i][j] = 0;
                }
                for (std::size_t j = k; j < cols; j++){
                    Type temp = lCoefs[k][j];
                    lCoefs[k][j] = c*lCoefs[k][j] + s*lCoefs[i][j];
                    lCoefs[i][j] = -s*temp + c*lCoefs[i][j];
                    if (std::abs(lCoefs[i][j]) < accuracy)
                        lCoefs[i][j] = 0;
                }
            }
        }
    }
    transposeMatrix(Q);
    return IS_QUADRATIC;
}

template<typename Type>
std::size_t transposeMatrix(std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return 0;
    for (std::size_t i = 0; i < rows; i++){
        for (std::size_t j = i + 1; j < cols; j++){
            Type temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
    return rows;
}

template<typename Type>
Type findResidual(const std::vector<std::vector<Type>> &lCoefs, const std::vector<Type> &rCoefs, const std::vector<Type> &solution){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    std::vector<Type> b1(rows);  // Правая часть после подстановки полученного решения
        for (std::size_t i = 0; i < rows; i++){
                Type sum = 0.0;
                for (std::size_t k = 0; k < cols; k++){
                    sum += lCoefs[i][k] * solution[k];
                }
                b1[i] = sum;
        }
        std::vector<Type> discrepancyVector(rows); // Вектор невязки
        for (std::size_t i = 0; i < rows; i++){
            discrepancyVector[i] = rCoefs[i] - b1[i];
        }
        Type discrepancy = 0.0; // Невязка
        for (std::size_t i = 0; i < rows; i++){
            discrepancy += discrepancyVector[i] * discrepancyVector[i];     
        }
    return std::sqrt(discrepancy);
}

template<typename Type>
Type findMatrixNorm1(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    Type norm1OfMatrix = 0;
    for (std::size_t j = 0; j < cols; j++){
        Type sum = 0.0;
        for (std::size_t i = 0; i < rows; i++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > norm1OfMatrix)
            norm1OfMatrix = sum;
    }
    return norm1OfMatrix;
}

template<typename Type>
Type findCond_1(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<std::vector<Type>> invMatrix; //Обратная к A матрица
    INVERTIBLE_FLAG flag = invertMatrix(matrix, invMatrix);
    Type norm1OfMatrix = findMatrixNorm1(matrix);
    Type norm1OfInvMatrix = 0;
    if (flag == IS_INVERTIBLE)
        norm1OfInvMatrix = findMatrixNorm1(invMatrix);
    else
        norm1OfInvMatrix = INFINITY;
    Type cond = norm1OfMatrix * norm1OfInvMatrix;
    return cond;
}

template<typename Type>
Type findMatrixNormInf(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    Type normInfOfMatrix = 0;
    for (std::size_t i = 0; i < rows; i++){
        Type sum = 0.0;
        for (std::size_t j = 0; j < cols; j++){
            sum += std::abs(matrix[i][j]);
        }
        if (sum > normInfOfMatrix)
            normInfOfMatrix = sum;
    }
    return normInfOfMatrix;
}

template<typename Type>
Type findCond_inf(const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<std::vector<Type>> invMatrix; //Обратная к A матрица
    INVERTIBLE_FLAG flag = invertMatrix(matrix, invMatrix);
    Type normInfOfMatrix = findMatrixNormInf(matrix);
    Type normInfOfInvMatrix = 0;
    if (flag == IS_INVERTIBLE)
        normInfOfInvMatrix = findMatrixNormInf(invMatrix);
    else
        normInfOfInvMatrix = INFINITY;
    Type cond = normInfOfMatrix * normInfOfInvMatrix;
    return cond;
}

template<typename Type>
INVERTIBLE_FLAG invertMatrix(const std::vector<std::vector<Type>> &inputMatrix, std::vector<std::vector<Type>> &resMatrix,
    SOLUTION_FLAG (*method)(std::vector<std::vector<Type>> &, std::vector<Type>&, std::vector<Type>&, Type accuracy)){
    std::size_t rows = inputMatrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = inputMatrix[0].size();
    else
        return NOT_INVERTIBLE;
    if (rows != cols)
        return NOT_INVERTIBLE;
    resMatrix.resize(rows);
    std::vector<std::vector<Type>> E(rows);
    for (std::size_t i = 0; i < rows; i++){
        E[i].resize(cols, 0);
        resMatrix[i].resize(cols);
    }
    for (std::size_t i = 0; i < rows; i++){
        E[i][i] = 1;
    }
    std::vector<Type> solution(rows);
    std::vector<std::vector<Type>> tempMatrix(rows);
    SOLUTION_FLAG flag;
    for (std::size_t i = 0; i < rows; i++){
        tempMatrix = inputMatrix;
        flag = (*method)(tempMatrix, E[i], solution, 1e-6);
        if (flag == NO_SOLUTION){
            for (std::size_t i = 0; i < rows; i++)
                resMatrix[i].clear();
            resMatrix.clear();
            return NOT_INVERTIBLE;
        }
        for (std::size_t k = 0; k < rows; k++)
            resMatrix[k][i] = solution[k];
    }
    return IS_INVERTIBLE;
}

template<typename Type>
Type norm1OfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type sum = 0;
    for (std::size_t i = 0; i < vector.size(); i++)
        sum += std::abs(vector[i]);
    return sum;
}

template<typename Type>
Type normInfOfVector(const std::vector<Type> &vector){
    if (!vector.size())
        return NAN;
    Type max = std::abs(vector[0]);
    for (std::size_t i = 1; i < vector.size(); i++)
        if (std::abs(vector[i]) > max)
            max = std::abs(vector[i]);
    return max;
}

template<typename Type>
Type findLowerBoundOfcond1(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1, Type delta2, Type delta3){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<Type> rCoefs1(rows), rCoefs2(rows), rCoefs3(rows);
    for (std::size_t i = 0; i < rows; i++){
        rCoefs1[i] = rCoefs[i] + delta1;
        rCoefs2[i] = rCoefs[i] + delta2;
        rCoefs3[i] = rCoefs[i] + delta3;
    }
    std::vector<std::vector<Type>> Q;
    findQMatrix(lCoefs, Q);
    transposeMatrix(Q);

    std::vector<Type> solution(rows), perturbSolution(rows), deltaSolution(rows);
    std::vector<Type> tempRCoefs;

    multiplyMatrix(Q, rCoefs, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, solution);

    multiplyMatrix(Q, rCoefs1, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs1(rows, delta1);
    Type lowerBound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs1));

    multiplyMatrix(Q, rCoefs2, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs2(rows, delta2);
    Type bound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs2));
    if (bound > lowerBound)
        lowerBound = bound;

    multiplyMatrix(Q, rCoefs3, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs3(rows, delta3);
    bound = (norm1OfVector(deltaSolution) * norm1OfVector(rCoefs))/(norm1OfVector(solution) * norm1OfVector(deltaRCoefs3));
    if (bound > lowerBound)
        lowerBound = bound;
    
    return lowerBound;
}

template<typename Type>
Type findLowerBoundOfcondInf(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs, Type delta1, Type delta2, Type delta3){
    std::size_t rows = lCoefs.size(); // Количество строк в СЛАУ
    std::size_t cols = 0;
    if (rows != 0)
        cols = lCoefs[0].size();
    else
        return NAN;
    if (rows != cols)
        return NAN;
    std::vector<Type> rCoefs1(rows), rCoefs2(rows), rCoefs3(rows);
    for (std::size_t i = 0; i < rows; i++){
        rCoefs1[i] = rCoefs[i] + delta1;
        rCoefs2[i] = rCoefs[i] + delta2;
        rCoefs3[i] = rCoefs[i] + delta3;
    }
    std::vector<std::vector<Type>> Q;
    findQMatrix(lCoefs, Q);
    transposeMatrix(Q);

    std::vector<Type> solution(rows), perturbSolution(rows), deltaSolution(rows);
    std::vector<Type> tempRCoefs;

    multiplyMatrix(Q, rCoefs, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, solution);

    multiplyMatrix(Q, rCoefs1, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs1(rows, delta1);
    Type lowerBound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs1));

    multiplyMatrix(Q, rCoefs2, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs2(rows, delta2);
    Type bound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs2));
    if (bound > lowerBound)
        lowerBound = bound;

    multiplyMatrix(Q, rCoefs3, tempRCoefs);
    qrMethod<Type>(lCoefs, tempRCoefs, perturbSolution);
    for (std::size_t i = 0; i < cols; i++){
        deltaSolution[i] = solution[i] - perturbSolution[i];
    }
    std::vector<Type> deltaRCoefs3(rows, delta3);
    bound = (normInfOfVector(deltaSolution) * normInfOfVector(rCoefs))/(normInfOfVector(solution) * normInfOfVector(deltaRCoefs3));
    if (bound > lowerBound)
        lowerBound = bound;
    
    return lowerBound;
}

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(const std::vector<std::vector<Type>> &matrix1, const std::vector<std::vector<Type>> &matrix2, std::vector<std::vector<Type>> &result){
    std::size_t rows1 = matrix1.size();
    std::size_t cols1 = 0;
    if (rows1 != 0)
        cols1 = matrix1[0].size();
    else
        return NOT_MULTIPLIED;
    std::size_t rows2 = matrix2.size();
    std::size_t cols2 = 0;
    if (rows2 != 0)
        cols2 = matrix2[0].size();
    else
        return NOT_MULTIPLIED;
    if (cols1 != rows2)
        return NOT_MULTIPLIED;
    result.resize(rows1);
    for (std::size_t i = 0; i < rows1; i++){
        result[i].resize(cols2);
    }
    std::size_t rows = result.size();
    std::size_t cols = result[0].size();
    for (size_t i = 0; i < rows; i++){
        for (size_t j = 0; j < cols; j++){
            Type sum = 0;
            for (size_t k = 0; k < cols1; k++){
                sum += matrix1[i][k] * matrix2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return IS_MULTIPLIED;
}

template<typename Type>
MULTIPLIED_FLAG multiplyMatrix(const std::vector<std::vector<Type>> &matrix, const std::vector<Type> &vec, std::vector<Type> &result){
    std::size_t rows1 = matrix.size();
    std::size_t cols = 0;
    if (rows1 != 0)
        cols = matrix[0].size();
    else
        return NOT_MULTIPLIED;
    std::size_t rows2 = vec.size();
    if (cols != rows2)
        return NOT_MULTIPLIED;
    result.resize(rows1);
    std::size_t rows = result.size();
    for (size_t i = 0; i < rows; i++){
        Type sum = 0;
        for (size_t k = 0; k < cols; k++){
            sum += matrix[i][k] * vec[k];
        }
        result[i] = sum;
    }
    return IS_MULTIPLIED;
}

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<std::vector<Type>> &matrix){
    std::size_t rows = matrix.size();
    std::size_t cols = 0;
    if (rows != 0)
        cols = matrix[0].size();
    else{
        os << 0;
        return os;
    }
    for (std::size_t i = 0; i < rows - 1; i++){
        for (std::size_t j = 0; j < cols; j++){
            os << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
    for (std::size_t j = 0; j < cols; j++){
            os << matrix[rows - 1][j] << ' ';
        }
    return os;
}

template<typename Type>
std::ostream& operator<<(std::ostream &os, const std::vector<Type> &vector){
    std::size_t rows = vector.size();
    if (!rows){
        os << 0;
        return os;
    }  
    os << "{ ";
    for (std::size_t i = 0; i < rows - 1; i++)
        os << vector[i] << ", ";
    os << vector[rows - 1] << ' ';
    os << '}';
    return os;
}