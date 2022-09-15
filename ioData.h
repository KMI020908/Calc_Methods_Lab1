#ifndef IO_DATA_H
#define IO_DATA_H

#include<fstream> 
#include<iomanip> 
#include<vector>

enum SOLUTION_FLAG{
    HAS_SOLUTION,   // 0
    NO_SOLUTION     // 1 
};

enum FILE_FLAG{
    NOT_OPEN, // 0
    IS_CLOSED // 1
};

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
FILE_FLAG writeData(const std::vector<Type> &solution, const std::string& OUT_FILE_PATH, SOLUTION_FLAG FLAG = HAS_SOLUTION){
	std::ofstream file;
	file.open(OUT_FILE_PATH);
	if (!file.is_open())
		exit(NOT_OPEN);
    if (FLAG == NO_SOLUTION){
        file << "The solution:" << '\n';
	    file << "Нет решений.";
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

#endif