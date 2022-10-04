#include "ioData.h"

template<typename Type>
FILE_FLAG generateRandomTest(size_t equationDim, Type minValue, Type maxValue, const std::string& FILE_PATH){
    std::ofstream file;
    file.open(FILE_PATH);
    if (!file.is_open())
        exit(NOT_OPEN);
    file << equationDim << '\n';
    PRNG generator;
    for (std::size_t i = 0; i < equationDim; i++){
        for (std::size_t j = 0; j < equationDim + 1; j++)
            file << getRandomNumber(generator, minValue, maxValue) << '\t';
        file << '\n';
    }
    file.close();
    return IS_CLOSED;
}

template<typename Type>
FILE_FLAG readData(std::vector<std::vector<Type>> &lCoefs, std::vector<Type> &rCoefs,  const std::string& IN_FILE_PATH) {
	std::ifstream file;
	file.open(IN_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
	std::size_t n;
	file >> n;
    std::vector<Type> hVec;
    hVec.reserve(n);
    Type hValue = 0;
    lCoefs.clear();
    rCoefs.clear();
    for (std::size_t i = 0; i < n; i++){
        for (std::size_t j = 0; j < n; j++){ 
            file >> hValue;
            hVec.push_back(hValue);    
        }
        file >> hValue;
        rCoefs.push_back(hValue);
        lCoefs.push_back(hVec);
        hVec.clear();
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeData(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, SOLUTION_FLAG FLAG){
	std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (FLAG == NO_SOLUTION){
        file << "The solution:" << '\n';
	    file << "No solution.";
        file.close();
        return IS_CLOSED;
    }
    file << "The solution:" << '\n';
	file << "X = " << solution;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeQRMatrix(const std::vector<std::vector<Type>> &Q, const std::vector<std::vector<Type>> &R, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t dimMatrix = Q.size(); 
    file << '\n' << '\n';
    file << "Q matrix:" << '\n';
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            file << Q[i][j] << '\t';    
        }
        file << '\n';
    }
    file << '\n';
    file << "R matrix:" << '\n';
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            file << R[i][j] << '\t';    
        }
        file << '\n';
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeMatrixMultiply(const std::vector<std::vector<Type>> &B, const std::vector<std::vector<Type>> &A, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    std::size_t dimMatrix = A.size();
    file << '\n' << '\n';
    file << "Multiplication of invert A and original A:" << '\n';
    for (std::size_t i = 0; i < dimMatrix; i++){
        for (std::size_t j = 0; j < dimMatrix; j++){
            Type sum = 0;
            for (std::size_t k = 0; k < dimMatrix; k++){
                sum += B[i][k]*A[k][j];
            }
            file << sum << '\t';
        }
        file << '\n';
    }
    file.close();
    return IS_CLOSED;
} 

template<typename Type>
FILE_FLAG writeResidual(Type residual, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << '\n' << '\n';
    file << "Residual = " << residual;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG addPerturbation(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type perturbation, SOLUTION_FLAG FLAG){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (FLAG == NO_SOLUTION){
        file << '\n' << '\n';
        file << "The solution:" << '\n';
	    file << "No solution.";
        file.close();
        return IS_CLOSED;
    }
    file << '\n' << '\n';
    file << "The solution after perturbation = " << perturbation << '\n';
	file << "X1 = " << "{ ";
    for (std::size_t i = 0; i < solution.size() - 1; i++)
        file << solution[i] << ", ";
    file << solution[solution.size() - 1] << ' ';
    file << '}';
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeConds(Type cond_1, Type cond_inf, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    file << '\n' << '\n';
    file << "cond_1 A = " << cond_1;
    file << '\n' << '\n';
    file << "cond_inf A = " << cond_inf;
    file.close();
    return IS_CLOSED;
} 

template<typename Type>
FILE_FLAG writeLowerBounds(Type lowerBound1, Type lowerBoundInf, const std::string& OUT_FILE_PATH){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    file << '\n' << '\n';
    file << lowerBound1 << " ≤ cond_1 A";
    file << '\n' << '\n';
    file << lowerBoundInf << " ≤ cond_inf A";
    file.close();
    return IS_CLOSED;
}