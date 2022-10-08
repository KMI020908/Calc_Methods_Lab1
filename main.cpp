#include<vector>
#include"methods.cpp"
#include"ioData.cpp"
#include"filePath.h"

// Процедура проверки алгоритмов
template<typename Type>
void checkTest(std::vector<std::vector<Type>> &lCoefSys, std::vector<Type> &rCoefSys,
const std::string &IN_FILE_PATH, const std::string &G_OUT_FILE_PATH, const std::string &QR_OUT_FILE_PATH, Type perturbation = 0.01){
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    std::vector<Type> solution;
    SOLUTION_FLAG flag = gaussMethodFull<Type>(lCoefSys, rCoefSys, solution);
    if (flag == HAS_SOLUTION){
        writeData<Type>(solution, G_OUT_FILE_PATH);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH); 
        Type cond_1 = findCond_1<Type>(lCoefSys);
        Type cond_inf = findCond_inf<Type>(lCoefSys);
        writeConds<Type>(cond_1, cond_inf, G_OUT_FILE_PATH);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type lowerBound1 = findLowerBoundOfcond1(lCoefSys, rCoefSys);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type lowerBoundInf = findLowerBoundOfcondInf(lCoefSys, rCoefSys);
        writeLowerBounds<Type>(lowerBound1, lowerBoundInf, G_OUT_FILE_PATH);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        std::vector<std::vector<Type>> invertlCoefSysGauss;
        std::vector<std::vector<Type>> invertlCoefSysQR;
        invertMatrix<Type>(lCoefSys, invertlCoefSysGauss, gaussMethod);
        invertMatrix<Type>(lCoefSys, invertlCoefSysQR, qrMethod);
        writeMatrixMultiplyInvA<Type>(invertlCoefSysGauss, lCoefSys, G_OUT_FILE_PATH, "with Gauss method");
        writeMatrixMultiplyInvA<Type>(invertlCoefSysQR, lCoefSys, G_OUT_FILE_PATH, "with QR method");
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type residual = findResidual<Type>(lCoefSys, rCoefSys, solution);
        writeResidual<Type>(residual, G_OUT_FILE_PATH);
        // Возбуждение = perturbation
        for (std::size_t i = 0; i < lCoefSys.size(); i++){
            rCoefSys[i] += perturbation;
        }
        flag = gaussMethod<Type>(lCoefSys, rCoefSys, solution);
        addPerturbation<Type>(solution, G_OUT_FILE_PATH, perturbation, flag);
        // Возбуждение = -perturbation
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        for (std::size_t i = 0; i < lCoefSys.size(); i++){
            rCoefSys[i] -= perturbation;
        }
        flag = gaussMethod<Type>(lCoefSys, rCoefSys, solution);
        addPerturbation<Type>(solution, G_OUT_FILE_PATH, -perturbation, flag);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    }
    else{
        writeData<Type>(solution, G_OUT_FILE_PATH, NO_SOLUTION);    
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH); 
        Type cond_1 = findCond_1<Type>(lCoefSys);
        Type cond_inf = findCond_inf<Type>(lCoefSys);
        writeConds<Type>(cond_1, cond_inf, G_OUT_FILE_PATH);
    }

    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    std::vector<std::vector<Type>> Q;
    findQMatrix<Type>(lCoefSys, Q);
    readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
    flag = qrMethod<Type>(lCoefSys, rCoefSys, solution); 
    if (flag == HAS_SOLUTION){
        writeData<Type>(solution, QR_OUT_FILE_PATH);
        writeQRMatrix<Type>(Q, lCoefSys, QR_OUT_FILE_PATH);
        writeMatrixMultiplyQR(Q, lCoefSys, QR_OUT_FILE_PATH);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type cond_1 = findCond_1<Type>(lCoefSys);
        Type cond_inf = findCond_inf<Type>(lCoefSys);
        writeConds<Type>(cond_1, cond_inf, QR_OUT_FILE_PATH);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type lowerBound1 = findLowerBoundOfcond1(lCoefSys, rCoefSys);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type lowerBoundInf = findLowerBoundOfcondInf(lCoefSys, rCoefSys);
        writeLowerBounds<Type>(lowerBound1, lowerBoundInf, QR_OUT_FILE_PATH);
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        std::vector<std::vector<Type>> invertlCoefSysGauss;
        std::vector<std::vector<Type>> invertlCoefSysQR;
        invertMatrix<Type>(lCoefSys, invertlCoefSysGauss, gaussMethod);
        invertMatrix<Type>(lCoefSys, invertlCoefSysQR, qrMethod);
        writeMatrixMultiplyInvA<Type>(invertlCoefSysGauss, lCoefSys, QR_OUT_FILE_PATH, "with Gauss method");
        writeMatrixMultiplyInvA<Type>(invertlCoefSysQR, lCoefSys, QR_OUT_FILE_PATH, "with QR method");
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        Type residual = findResidual<Type>(lCoefSys, rCoefSys, solution);
        writeResidual<Type>(residual, QR_OUT_FILE_PATH);
        // Возбуждение = perturbation
        for (std::size_t i = 0; i < lCoefSys.size(); i++){
            rCoefSys[i] += perturbation;
        }
        flag = qrMethod<Type>(lCoefSys, rCoefSys, solution);
        addPerturbation<Type>(solution, QR_OUT_FILE_PATH, perturbation, flag);
        // Возбуждение = -perturbation
        readData<Type>(lCoefSys, rCoefSys, IN_FILE_PATH);
        for (std::size_t i = 0; i < lCoefSys.size(); i++){
            rCoefSys[i] -= perturbation;
        }
        flag = qrMethod<Type>(lCoefSys, rCoefSys, solution);
        addPerturbation<Type>(solution, QR_OUT_FILE_PATH, -perturbation, flag);
    }
    else{
        writeData<Type>(solution, QR_OUT_FILE_PATH, NO_SOLUTION);    
        Type cond_1 = findCond_1<Type>(lCoefSys);
        Type cond_inf = findCond_inf<Type>(lCoefSys);
        writeConds<Type>(cond_1, cond_inf, QR_OUT_FILE_PATH);
    }
}

template<typename Type>
void temp_main(){

    std::vector<std::vector<Type>> lCoefSys; // Матрица левых коэффициентов
    std::vector<Type> rCoefSys; // Вектор правых коэффициентов

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_1, G_OUT_FILE_PATH_1, QR_OUT_FILE_PATH_1);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_2, G_OUT_FILE_PATH_2, QR_OUT_FILE_PATH_2);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_3, G_OUT_FILE_PATH_3, QR_OUT_FILE_PATH_3);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_4, G_OUT_FILE_PATH_4, QR_OUT_FILE_PATH_4);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_5, G_OUT_FILE_PATH_5, QR_OUT_FILE_PATH_5);

    generateRandomTest<Type>(4, 1.0, 20.0, IN_FILE_PATH_6);
    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_6, G_OUT_FILE_PATH_6, QR_OUT_FILE_PATH_6);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_7, G_OUT_FILE_PATH_7, QR_OUT_FILE_PATH_7);

    checkTest(lCoefSys, rCoefSys, IN_FILE_PATH_8, G_OUT_FILE_PATH_8, QR_OUT_FILE_PATH_8);
}

int main(){
    temp_main<double>();

    std::vector<std::vector<float>> A; // Матрица левых коэффициентов
    std::vector<float> b; // Вектор правых коэффициентов
    std::size_t dim = 8;
    std::vector<float> tempVec;
    for (size_t i = 1; i <= dim; i++){
        for (size_t j = 1; j <= dim; j++){
            tempVec.push_back(1.0/(i + j - 1));
        }
        A.push_back(tempVec);
        tempVec.clear();
    }
    for (size_t i = 1; i <= dim; i++){
        b.push_back(1.0 / i);
    }

    std::cout << "Числа обусловленности для плохого теста: " << findCond_1(A) << '\t' << findCond_inf(A) << '\n';
    std::vector<float> X; // Решение
    SOLUTION_FLAG flag = gaussMethod<float>(A, b, X, 1e-7);
    if (flag == NO_SOLUTION)
        std::cout << "Нет решения плохого теста" << '\n';
    else
        for (size_t i = 0; i < dim; i++){
           std::cout << X[i] << '\n';
        }

    return 0;
}