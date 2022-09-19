#include "ioData.h"

template<typename Type>
FILE_FLAG generateRandomTest(size_t equationDim, Type minValue, Type maxValue, const std::string& FILE_PATH){
    std::ofstream file;
    file.open(FILE_PATH);
    if (!file.is_open())
        exit(NOT_OPEN);
    file << equationDim << '\n';
    PRNG generator;
    for (size_t i = 0; i < equationDim; i++){
        for (size_t j = 0; j < equationDim + 1; j++)
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
	size_t n;
	file >> n;
    std::vector<Type> hVec;
    hVec.reserve(n);
    Type hValue = 0;
    lCoefs.clear();
    rCoefs.clear();
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < n; j++){ 
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
	file << "X = " << "{ ";
    for (size_t i = 0; i < solution.size() - 1; i++)
        file << solution[i] << ", ";
    file << solution[solution.size() - 1] << ' ';
    file << '}';
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeQRMatrix(const std::vector<std::vector<Type>> &Q, const std::vector<std::vector<Type>> &R, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    size_t dimMatrix = Q.size(); 
    file << '\n' << '\n';
    file << "Q matrix:" << '\n';
    for (size_t i = 0; i < dimMatrix; i++){
        for (size_t j = 0; j < dimMatrix; j++){
            file << Q[i][j] << '\t';    
        }
        file << '\n';
    }
    file << '\n';
    file << "R matrix:" << '\n';
    for (size_t i = 0; i < dimMatrix; i++){
        for (size_t j = 0; j < dimMatrix; j++){
            file << R[i][j] << '\t';    
        }
        file << '\n';
    }
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeDiscrepancy(const std::vector<Type> &discrepancyVec, Type discrepancy, const std::string& OUT_FILE_PATH){
	std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN); 
    file << '\n' << '\n';
    file << "Discrepancy vector:" << '\n';
    file << "b - b1 = " << "{ ";
    for (size_t i = 0; i < discrepancyVec.size() - 1; i++)
        file << discrepancyVec[i] << ", ";
    file << discrepancyVec[discrepancyVec.size() - 1] << ' ';
    file << '}';
    file << '\n' <<'\n';
    file << "Discrepancy = " << discrepancy;
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG addDisturbance(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, Type disturbance, SOLUTION_FLAG FLAG){
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
    file << "The solution after disturbance = " << disturbance << '\n';
	file << "X1 = " << "{ ";
    for (size_t i = 0; i < solution.size() - 1; i++)
        file << solution[i] << ", ";
    file << solution[solution.size() - 1] << ' ';
    file << '}';
	file.close();
	return IS_CLOSED;
}

template<typename Type>
FILE_FLAG writeConds(Type cond_1, Type cond_inf, const std::string& OUT_FILE_PATH, SOLUTION_FLAG FLAG){
    std::ofstream file;
	file.open(OUT_FILE_PATH, std::ios::app);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (FLAG == NO_SOLUTION){
        file << '\n' << '\n';
        file << "cond A = inf";
        return IS_CLOSED;
    }
    file << '\n' << '\n';
    file << "cond_1 A = " << cond_1;
    file << '\n' << '\n';
    file << "cond_inf A = " << cond_inf;
    file.close();
    return IS_CLOSED;
} 