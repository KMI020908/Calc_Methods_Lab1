#include<vector>
#include"ioData.h"
#include"methods.cpp"
#include"filePath.h"

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys, const std::string IN_FILE_PATH, 
const std::string G_OUT_FILE_PATH, const std::string QR_OUT_FILE_PATH){
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    gaussMethod<Type>(lCoefSys, rCoefSys, G_OUT_FILE_PATH);
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    qrMethod<Type>(lCoefSys, rCoefSys, QR_OUT_FILE_PATH);
}

int main(){

    std::vector<std::vector<double>> lCoefSys; // Матрица левых коэффициентов
    std::vector<double> rCoefSys; // Вектор правых коэффициентов

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_1, G_OUT_FILE_PATH_1, QR_OUT_FILE_PATH_1);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_2, G_OUT_FILE_PATH_2, QR_OUT_FILE_PATH_2);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_3, G_OUT_FILE_PATH_3, QR_OUT_FILE_PATH_3);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_4, G_OUT_FILE_PATH_4, QR_OUT_FILE_PATH_4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_5, G_OUT_FILE_PATH_5, QR_OUT_FILE_PATH_5);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_6, G_OUT_FILE_PATH_6, QR_OUT_FILE_PATH_6);

    return 0;
}