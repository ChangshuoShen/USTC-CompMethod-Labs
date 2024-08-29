#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <cassert>
#include <iomanip> //设置输出流精度

// 代表矩阵和向量的类型
typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

// 定义出无穷范数，LU分解，解线性方程，反幂法
int max_elem_index(const Vector &vec);
void LU_decomposition(const Matrix &A, Matrix &L, Matrix &U);
Vector solveEquation(const Matrix &L, const Matrix &U, const Vector &b);
Vector inversePowerMethod(const Matrix &A, const Vector &initialGauss, double tolerance);
void printMatrix(Matrix &matrix);
void printVector(Vector &vector);

int main() {
    Matrix A = {{1.0/9, 1.0/8, 1.0/7, 1.0/6, 1.0/5},
                {1.0/8, 1.0/7, 1.0/6, 1.0/5, 1.0/4},
                {1.0/7, 1.0/6, 1.0/5, 1.0/4, 1.0/3},
                {1.0/6, 1.0/5, 1.0/4, 1.0/3, 1.0/2},
                {1.0/5, 1.0/4, 1.0/3, 1.0/2, 1.0/1}};

    printMatrix(A);
    
    std::cout<<std::endl;
    Matrix B = {{4.0, -1.0, 1.0, 3.0},
                {16.0, -2.0, -2.0, 5.0},
                {16.0, -3.0, -1.0, 7.0},
                {6.0, -4.0, 2.0, 9.0}};
    printMatrix(B);
    std::cout<<std::endl;
    Vector initialVector_A = {1, 1, 1, 1, 1};
    Vector initialVector_B = {1, 1, 1, 1};
    // 允许误差
    double tolerance = 1e-6;

    Vector eigenVector_A = inversePowerMethod(A, initialVector_A, tolerance);
    Vector eigenVector_B = inversePowerMethod(B, initialVector_B, tolerance);

    std::cout<<"A'a eigen vector: ";
    printVector(eigenVector_A);
    std::cout<<"B's eigen vector: ";
    printVector(eigenVector_B);
    return 0;
}

// 反幂法
Vector inversePowerMethod(const Matrix &A, const Vector &initialGuess, double tolerance) {
    assert(A.size() == A[0].size() && "Matrix A must be square");
    assert(A.size() == initialGuess.size() && "Initial guess vector must have the same size as A");
    int iteration_cnt = 0;
    Vector eigenVector = initialGuess;
    double lambda = 0.0;
    double prevLambda = 0.0;
    Matrix L, U;
    LU_decomposition(A, L, U);
    do 
    {
        prevLambda = lambda;

        // 计算Ax_{k+1} = x_k，也就是LUx = x
        Vector nextEigenVector = solveEquation(L, U, eigenVector);
        std::cout<<"next eigen vector:\t";
        printVector(nextEigenVector);
        // 使用无穷范数规范化
        int max_index = max_elem_index(nextEigenVector);
        // 规范化之前先记录按模最大的元素，在迭代后期就是特征值lambda
        lambda = nextEigenVector[max_index];
        // double norm = std::abs(nextEigenVector[max_index]);
        for (double &element : nextEigenVector)
        {
            element /= lambda;
        }
        std::cout<<"normalized vector:";
        printVector(nextEigenVector);
        std::cout<<"Lambda - prevLambda:\t"<<lambda - prevLambda<<std::endl;
        std::cout<<"lambda this loop(before processed): "<<lambda<<"iteration times so far: "<<++iteration_cnt<<std::endl;

        eigenVector = nextEigenVector;

    } while (std::abs(lambda - prevLambda) > tolerance && iteration_cnt <= 1000);

    lambda = 1.0 / lambda;

    std::cout<<"Real Lambda: "<<lambda<<std::endl<<"eigenVector:";
    printVector(eigenVector);
    
    return eigenVector;
}

int max_elem_index(const Vector &vec) {
    // 返回按模最大元素的下标
    // 先不定义无穷范数，可能会有问题
    int n = vec.size();
    double max_elem = std::abs(vec[0]);
    int max_index = 0;
    for(int i = 1; i < n; i++)
    {
        if(std::abs(vec[i]) > max_elem)
        {
            max_elem = std::abs(vec[i]);
            max_index = i;
        }
    }
    return max_index;
}


// 执行 LU 分解
void LU_decomposition(const Matrix &A, Matrix &L, Matrix &U) {
    int n = A.size();
    L = Matrix(n, Vector(n, 0));
    U = Matrix(n, Vector(n, 0));

    // 先初始化L，对角元都是1
    for (int i = 0; i < n; ++i)
    {
        L[i][i] = 1.0;
    }
        
    // Perform LU decomposition
    for (int i = 0; i < n; ++i) 
    {
        // Compute U matrix
        for (int k = i; k < n; ++k) 
        {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        // Compute L matrix
        for (int k = i + 1; k < n; ++k) 
        {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[k][j] * U[j][i];
            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }
    std::cout<<"LU decomposition finished"<<std::endl;
    std::cout<<"L: "<<std::endl;
    printMatrix(L);
    std::cout<<"U: "<<std::endl;
    printMatrix(U);
}

// 使用分解好的LU计算Ax = b，防止没必要的多次分解
Vector solveEquation(const Matrix &L, const Matrix &U, const Vector &b) {
    int n = L.size();

    // Solve Ly = b
    Vector y(n);
    for (int i = 0; i < n; ++i) 
    {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
            sum += L[i][j] * y[j];
        y[i] = b[i] - sum;
    }

    // Solve Ux = y
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) 
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j)
            sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

void printMatrix(Matrix &matrix) {
    int n = matrix.size();
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            std::cout << std::setprecision(10) << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void printVector(Vector &vector) {
    int n = vector.size();
    std::cout << "(";
    for(int i = 0; i < n; i++) {
        std::cout << std::setprecision(10) << vector[i] << "\t";
    }
    std::cout << ")" << std::endl;
}
