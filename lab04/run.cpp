#include<iostream>
#include<vector>
#include<ctime>
#include<cmath>
#include<cstddef>
#include<random>
#include<fstream>
#include<sstream>
#include<string>

typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

Matrix matrixMultiplication(const Matrix& A, const Matrix& B);
Matrix transpose(const Matrix& A);
Matrix generateRandomMatrix(size_t rows, size_t cols);
void printMatrix(const Matrix& matrix);
void normalize_matrix(Matrix& matrix);
Matrix reshapeMatrix(const Matrix& input, size_t targetRows, size_t targetCols);
void testJacobi();
Matrix read_iris();
void write_data(const Matrix& input);
Matrix removeLastColumn(const Matrix& matrix);
void addLastColumn(Matrix& matrix1, const Matrix& matrix2);
bool isclose(const Matrix& A, const Matrix& B, double tolerance);




class SymmetricMatrix
{
// 对称矩阵类定义
private:
    Matrix matrix;
    Matrix rotationMatrix;
    size_t size;

public:
    // 构造函数
    SymmetricMatrix(Matrix A) : matrix(A), size(A.size()) {
        // 初始化旋转矩阵为单位矩阵，之后对matrix进行操作的时候对这个矩阵做行变换，记录总的旋转矩阵
        rotationMatrix = Matrix(size, Vector(size, 0));
        for (size_t i = 0; i < size; i++) 
        {
            rotationMatrix[i][i] = 1.0;
        }
    }

    double getSquareSumOfNonDiagonalElements() const 
    {
        double squareSum = 0.0;
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = i + 1; j < size; ++j) { // Only iterate over the upper triangular part
                squareSum += 2 * matrix[i][j] * matrix[i][j];
            }
        }
        return squareSum;
    }

    // Jacobi 变换
    void JacobiTransform() 
    {
        while (true) 
        {
            size_t p = 0, q = 1;
            // 找到最大的非对角元素的下标p, q
            for (size_t i = 0; i < size; ++i) 
            {
                for (size_t j = i + 1; j < size; ++j) 
                {
                    if (std::abs(matrix[i][j]) > std::abs(matrix[p][q])) 
                    {
                        p = i;
                        q = j;
                    }
                }
            }

            // 如果最大的非对角元素小于精度阈值，则停止迭代
            std::cout<<"square sum of non-diagonal elems:\t"<<getSquareSumOfNonDiagonalElements()<<std::endl;
            if (std::abs(matrix[p][q]) < 1e-6) 
            {
                break;
            }

            // 计算旋转角度
            double theta;
            if (matrix[p][p] == matrix[q][q]) 
            {
                theta = M_PI / 4;
            } 
            else 
            {
                theta = 0.5 * std::atan(2 * matrix[p][q] / (matrix[q][q] - matrix[p][p]));
            }

            // 计算旋转矩阵的元素
            double c = std::cos(theta);
            double s = std::sin(theta);

            // 更新旋转矩阵的第p行和第q行，将前面计算得到的旋转变换作用到旋转矩阵上
            // A = JDJ^T上，最终的旋转矩阵就是J
            for (size_t i = 0; i < size; i++) 
            {
                double temp1 = rotationMatrix[p][i];
                double temp2 = rotationMatrix[q][i];
                rotationMatrix[p][i] = c * temp1 - s * temp2;
                rotationMatrix[q][i] = s * temp1 + c * temp2;
            }

            // 更新矩阵元素
            for (size_t i = 0; i < size; i++) 
            {
                if (i != p && i != q) 
                {
                    double aip = matrix[i][p];
                    double aiq = matrix[i][q];
                    matrix[i][p] = c * aip - s * aiq;
                    matrix[i][q] = s * aip + c * aiq;
                    matrix[p][i] = matrix[i][p];
                    matrix[q][i] = matrix[i][q];
                }
            }

            double app = matrix[p][p];
            double aqq = matrix[q][q];
            double apq = matrix[p][q];
            matrix[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
            matrix[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
            matrix[p][q] = 0; // 对称矩阵的非对元素变为0
            matrix[q][p] = 0;
        }
        // 循环结束之后，按照对角元素大小重排矩阵D使之DESC排列，此时需要同时做行列变换，同时对旋转矩阵做相同的行变换
         // 循环结束之后，按照对角元素大小重排矩阵D使之DESC排列
        for (size_t i = 0; i < size; i++) 
        {
            for (size_t j = i + 1; j < size; j++) 
            {
                if (matrix[i][i] < matrix[j][j]) 
                {
                    // 同时交换旋转矩阵的行和列
                    for (size_t k = 0; k < size; k++) 
                    {
                        std::swap(matrix[k][i], matrix[k][j]);
                    }
                    for (size_t k = 0; k < size; k++)
                    {
                        std::swap(matrix[i][k], matrix[j][k]);

                        std::swap(rotationMatrix[i][k], rotationMatrix[j][k]);
                    }
                }
            }
        }
        rotationMatrix = transpose(rotationMatrix);
    }

    Matrix getRotationMatrix()
    {
        return rotationMatrix;
    }

    Matrix getDiagram()
    {
        return matrix;
    }

    Matrix checkIt()
    {
        Matrix result = matrixMultiplication(matrixMultiplication(rotationMatrix, matrix), transpose(rotationMatrix));
        return result;
    }


    // 打印矩阵
    void printInformation() 
    {
        std::cout<<"transformed metrix:"<<std::endl;
        printMatrix(matrix);
        
        std::cout<<"Rotation Matrix:"<<std::endl;
        printMatrix(rotationMatrix);
        
        std::cout<<"check if sufficient"<<std::endl;
        Matrix output = checkIt();
        printMatrix(output);
    }
};


class SingularValueDecomposition
{
private:
    Matrix A; // 原始矩阵 A
    Matrix U;
    Matrix V;
    Matrix sigma;

public:
    SingularValueDecomposition(Matrix input_A): A(input_A) {}

    void computeSVD() 
    {
        // 记录A矩阵的行数和列数，便于后面调用
        size_t rows = A.size();
        size_t cols = A[0].size();

        // 在这里实现 SVD 分解的步骤
        // 可以直接使用成员变量 A 来访问原始矩阵

        // 1. 计算 A * A^T, A^T * A
        Matrix AAT = matrixMultiplication(A, transpose(A));
        Matrix ATA = matrixMultiplication(transpose(A), A);

        std::cout<<"A * A^T:"<<std::endl;
        printMatrix(AAT);
        std::cout<<"A^T * A:"<<std::endl;
        printMatrix(ATA);

        // 3. 对以上两个矩阵进行 Jacobi 变换，得到 U 和 V
        // 先创建两个对象
        SymmetricMatrix symAAT(AAT);
        SymmetricMatrix symATA(ATA);

        // 使用Jacobi计算U,V
        symAAT.JacobiTransform();
        symATA.JacobiTransform();

        U = symAAT.getRotationMatrix();
        V = symATA.getRotationMatrix();

        std::cout<<"for AAT:"<<std::endl;
        symAAT.printInformation();
        std::cout<<"for ATA:"<<std::endl;
        symATA.printInformation();

        // 4. 从对角化后的矩阵中提取奇异值，构建 sigma 矩阵
        Matrix tmp_sigma = symAAT.getDiagram();
        std::cout<<"tmp sigma:"<<std::endl;
        printMatrix(tmp_sigma);
        for(size_t i = 0; i < rows; i++)
        {   
            tmp_sigma[i][i] = std::sqrt(tmp_sigma[i][i]);
            for(size_t j = i+1; j < rows; j++)
            {
                tmp_sigma[i][j] = tmp_sigma[j][i] = 0;
            }
        }
        std::cout<<std::endl<<std::endl;
        printMatrix(tmp_sigma);
        sigma = reshapeMatrix(tmp_sigma, rows, cols);
        std::cout<<"sigma for SVD:"<<std::endl;
        printMatrix(sigma);
        
    }

    void checkIt()
    {
        std::cout<<"U:"<<std::endl;
        printMatrix(U);
        std::cout<<"V:"<<std::endl;
        printMatrix(V);
        std::cout<<"sigma:"<<std::endl;
        printMatrix(sigma);
        std::cout<<"U * sigma * V^T:"<<std::endl;
        Matrix USV = matrixMultiplication(matrixMultiplication(U, sigma), transpose(V));
        printMatrix(USV);
        std::cout<<"if the USV == A:\t"<<isclose(A, USV, 1e-5)<<std::endl;
    }

};

int main()
{
    std::cout << "~~~Lab 04 Executing~~~" << std::endl;
    testJacobi();

    std::cout<<"SVD:"<<std::endl<<std::endl;

    //生成随机矩阵
    Matrix A = generateRandomMatrix(4, 3);
    std::cout << "Original matrix:" << std::endl;
    printMatrix(A);
    SingularValueDecomposition SVD_A(A);
    SVD_A.computeSVD();
    SVD_A.checkIt();
    
    std::cout<<"~~~SVD over~~~"<<std::endl;



    // Matrix original_data = read_iris();
    // // std::cout<<"original data"<<std::endl;
    // // printMatrix(original_data);
    // Matrix processed_data = removeLastColumn(original_data);
    // // normalize_matrix(processed_data);
    // // std::cout<<"remove the last column & normalize it :"<<std::endl;
    // // printMatrix(processed_data);
    // Matrix DTD = matrixMultiplication(transpose(processed_data), processed_data);
    // std::cout<<"DTD:"<<std::endl;
    // printMatrix(DTD);
    
    // SymmetricMatrix iris_pca(DTD);

    // iris_pca.JacobiTransform();
    // iris_pca.printInformation();
    // Matrix rotate_M = iris_pca.getRotationMatrix();
    // printMatrix(rotate_M);
    // // printMatrix(matrixMultiplication(rotate_M, transpose(rotate_M))); // 检查旋转矩阵是否是单位正交阵
    // Matrix diag = iris_pca.getDiagram();
    // std::cout<<"Eign Matrix:"<<std::endl;
    // printMatrix(rotate_M);
    // std::cout<<"diagram:"<<std::endl;

    // std::cout<<"transformed matrix"<<std::endl;
    // Matrix transformed_matrix = matrixMultiplication(processed_data, rotate_M);
    // // normalize_matrix(transformed_matrix);
    // // printMatrix(transformed_matrix);
    // addLastColumn(transformed_matrix, original_data);
    // printMatrix(transformed_matrix);
    // write_data(transformed_matrix);
    // std::cout<<"~~~Ending~~~"<<std::endl;
    
    return 0;
}


void testJacobi()
{
    // 生成随机矩阵
    Matrix A = generateRandomMatrix(4, 3);
    std::cout<<"~~~begin testing jacobi~~~"<<std::endl;
    std::cout << "Original matrix:" << std::endl;
    printMatrix(A);

    // 计算 A * A^T 得到对称矩阵
    Matrix AAT = matrixMultiplication(A, transpose(A));
    std::cout<<"AA^T:"<<std::endl;
    printMatrix(AAT);

    // 创建对称矩阵对象，并进行 Jacobi 变换
    SymmetricMatrix symMatrix(AAT);
    symMatrix.JacobiTransform();

    std::cout << "Transformed matrix:" << std::endl;
    symMatrix.printInformation(); // 打印变换后的矩阵
}


Matrix generateRandomMatrix(size_t rows, size_t cols) 
{
    std::mt19937 rng(std::time(nullptr)); // 以当前时间作为随机数种子
    std::uniform_real_distribution<double> dist(0.0, 1.0); // 生成 [0,1) 之间的均匀分布随机数
    Matrix matrix(rows, Vector(cols));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = dist(rng); // 使用随机数生成器生成随机数
        }
    }
    return matrix;
}


void normalize_matrix(Matrix& matrix)
{
    // 遍历矩阵的每一列
    for (size_t j = 0; j < matrix[0].size(); ++j)
    {
        double column_mean = 0.0;

        // 计算当前列的均值
        for (size_t i = 0; i < matrix.size(); ++i) 
        {
            column_mean += matrix[i][j];
        }
        column_mean /= matrix.size();

        // 减去当前列的均值
        for (size_t i = 0; i < matrix.size(); ++i) 
        {
            matrix[i][j] -= column_mean;
        }
    }
}


void printMatrix(const Matrix& matrix) 
{
    for (const auto& row : matrix) {
        for (double element : row) {
            std::cout << element << "\t\t";
        }
        std::cout << std::endl;
    }
}

Matrix transpose(const Matrix& A) 
{
    // 矩阵转置函数
    size_t m = A.size();
    size_t n = A[0].size();
    Matrix result(n, Vector(m, 0));

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[j][i] = A[i][j];
        }
    }

    return result;
}


Matrix matrixMultiplication(const Matrix& A, const Matrix& B) 
{
    // 矩阵乘法函数
    size_t m = A.size();
    size_t n = B[0].size();
    size_t p = B.size();
    Matrix result(m, Vector(n, 0));

    for (size_t i = 0; i < m; ++i) 
    {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < p; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

Matrix reshapeMatrix(const Matrix& input, size_t targetRows, size_t targetCols) {
    size_t rows = input.size();
    size_t cols = input[0].size();

    Matrix result(targetRows, Vector(targetCols, 0.0)); // 初始化目标矩阵

    // 遍历输入矩阵的行和列，将元素复制到目标矩阵中
    for (size_t i = 0; i < std::min(rows, targetRows); ++i) {
        for (size_t j = 0; j < std::min(cols, targetCols); ++j) {
            result[i][j] = input[i][j];
        }
    }

    return result;
}

Matrix read_iris()
{
     // 打开文件
    std::ifstream file("iris.txt");
    
    // 检查文件是否成功打开
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
    }

    Matrix data; // 存储读取的数据

    std::string line;
    // 逐行读取文件内容
    while (std::getline(file, line)) 
    {
        Vector row; // 存储当前行的数据
        std::istringstream iss(line);
        std::string token;
        // 从当前行中逐个读取数字，并添加到当前行的vector中
        while (std::getline(iss, token, ',')) 
        {
            double num = std::stoi(token); // 将分割得到的字符串转换为整数
            row.push_back(num);
        }
        // 将当前行的vector添加到整个数据的vector中
        data.push_back(row);
    }
    // 关闭文件
    file.close();
    return data;
}

void write_data(const Matrix& input)
{
    // 打开文件流
    std::ofstream outFile("transformed_data.txt");

    // 检查文件是否成功打开
    if (!outFile.is_open()) 
    {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    // 写入数据到文件
    for (const auto& row : input) 
    {
        for (double value : row) 
        {
            outFile << value << " ";
        }
        outFile << "\n"; // 写入换行符
    }

    // 关闭文件流
    outFile.close();
}

Matrix removeLastColumn(const Matrix& matrix) 
{
    Matrix newMatrix;
    for (const auto& row : matrix) {
        Vector newRow(row.begin(), row.end() - 1); // 复制当前行的数据，但不包含最后一列
        newMatrix.push_back(newRow);
    }
    return newMatrix;
}

void addLastColumn(Matrix& matrix1, const Matrix& matrix2)
{
    size_t cols = matrix1[0].size();
    for(size_t i = 0; i < matrix1.size(); ++i)
    {
        matrix1[i].push_back(matrix2[i][cols]);
    }
}

bool is_num_close(double a, double b, double tolerance) {
    return std::abs(a - b) <= tolerance;
}

bool isclose(const Matrix& A, const Matrix& B, double tolerance) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        return false; // Matrices must have the same dimensions
    }

    size_t rows = A.size();
    size_t cols = A[0].size();

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (!is_num_close(A[i][j], B[i][j], tolerance)) {
                return false;
            }
        }
    }

    return true;
}