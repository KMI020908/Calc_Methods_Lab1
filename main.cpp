#include<vector>
#include"ioData.h"
#include"methods.cpp"
#include"filePath.h"
int main()
{
    std::vector<std::vector<double>> lCoefSys1; // Матрица левых коэффициентов
    std::vector<double> rCoefSys1; // Вектор правых коэффициентов

    std::vector<std::vector<double>> lCoefSys2; // Матрица левых коэффициентов
    std::vector<double> rCoefSys2; // Вектор правых коэффициентов

    std::vector<std::vector<double>> lCoefSys3; // Матрица левых коэффициентов
    std::vector<double> rCoefSys3; // Вектор правых коэффициентов

    std::vector<std::vector<double>> lCoefSys4; // Матрица левых коэффициентов
    std::vector<double> rCoefSys4; // Вектор правых коэффициентов

    std::vector<std::vector<double>> lCoefSys5; // Матрица левых коэффициентов
    std::vector<double> rCoefSys5; // Вектор правых коэффициентов

    std::vector<std::vector<double>> lCoefSys6; // Матрица левых коэффициентов
    std::vector<double> rCoefSys6; // Вектор правых коэффициентов

    readData<double>(lCoefSys1, rCoefSys1, IN_FILE_PATH_1);
    gaussMethod<double>(lCoefSys1, rCoefSys1, G_OUT_FILE_PATH_1);
    readData<double>(lCoefSys1, rCoefSys1, IN_FILE_PATH_1);
    qrMethod<double>(lCoefSys1, rCoefSys1, QR_OUT_FILE_PATH_1);

    readData<double>(lCoefSys2, rCoefSys2, IN_FILE_PATH_2);
    gaussMethod<double>(lCoefSys2, rCoefSys2, G_OUT_FILE_PATH_2);
    readData<double>(lCoefSys2, rCoefSys2, IN_FILE_PATH_2);
    qrMethod<double>(lCoefSys2, rCoefSys2, QR_OUT_FILE_PATH_2);

    readData<double>(lCoefSys3, rCoefSys3, IN_FILE_PATH_3);
    gaussMethod<double>(lCoefSys3, rCoefSys3, G_OUT_FILE_PATH_3);
    readData<double>(lCoefSys3, rCoefSys3, IN_FILE_PATH_3);
    qrMethod<double>(lCoefSys3, rCoefSys3, QR_OUT_FILE_PATH_3);

    readData(lCoefSys4, rCoefSys4, IN_FILE_PATH_4);
    gaussMethod(lCoefSys4, rCoefSys4, G_OUT_FILE_PATH_4);
    readData(lCoefSys4, rCoefSys4, IN_FILE_PATH_4);
    qrMethod(lCoefSys4, rCoefSys4, QR_OUT_FILE_PATH_4);

    readData<double>(lCoefSys5, rCoefSys5, IN_FILE_PATH_5);
    gaussMethod<double>(lCoefSys5, rCoefSys5, G_OUT_FILE_PATH_5);
    readData<double>(lCoefSys5, rCoefSys5, IN_FILE_PATH_5);
    qrMethod<double>(lCoefSys5, rCoefSys5, QR_OUT_FILE_PATH_5);

    generateRandomTest<double>(4, 1.0, 10.0, IN_FILE_PATH_6);
    readData<double>(lCoefSys6, rCoefSys6, IN_FILE_PATH_6);
    gaussMethod<double>(lCoefSys6, rCoefSys6, G_OUT_FILE_PATH_6);
    readData<double>(lCoefSys6, rCoefSys6, IN_FILE_PATH_6);
    qrMethod<double>(lCoefSys6, rCoefSys6, QR_OUT_FILE_PATH_6);

    return 0;
}