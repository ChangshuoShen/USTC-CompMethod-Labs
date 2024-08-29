#include<iostream>
#include<vector>
#include<cmath>
#include<fstream>
// #include<cstdlib>

void printMatrix(const std::vector<std::vector<double>>& matrix);
std::vector<double> gaussElimination(std::vector<std::vector<double>>& matrix);
double norm(const std::vector<double>& x);
std::vector<double> gaussSeidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double epsilon, int maxIterations);
double analyticalSolution(double x, double epsilon, double a)
{
    return ((1 - a) / (1 - exp(-1 / epsilon))) * (1 - exp(-x / epsilon)) + a * x;
}
void test();


int main()
{
    std::cout<<"Linear Equations?   启动！！！"<<std::endl;
    test();
    double epsilon = 1.0;
    double a = 0.5;
    int n = 100;
    // 计算步长h
    double h = 1.0 / n;
    // 构建系数矩阵和右侧向量
    std::vector<std::vector<double>> A(n - 1, std::vector<double>(n - 1, 0.0));
    std::vector<double> b(n - 1, 0.0);
    // 设置这个三对角矩阵
    for (int i = 0; i < n - 1; i++) {
        if (i > 0) 
        {
            A[i][i - 1] = epsilon;
        }
        A[i][i] = -(2 * epsilon + h);

        if (i < n - 2) 
        {
            A[i][i + 1] = epsilon + h;
        }
        // 右侧矩阵
        // b[i] = (i == 0) ? a * h * h: a * h * h - epsilon;
        b[i] = a * h * h;
    }
    b[n-2] -= (epsilon + h);

    // 这里加一个增广矩阵
    std::vector<std::vector<double>> EA = A;
    for(int i = 0; i < n - 1; i++)
    {
        EA[i].push_back(b[i]);
    }
    // printMatrix(EA);
    // 使用列主元 Gauss 消元法求解
    std::vector<double> solution_gauss = gaussElimination(EA);
    // 使用 Gauss-Seidel 迭代法求解
    std::vector<double> solution_seidel = gaussSeidel(A, b, 5e-5, 10000);
    // 计算解析解和计算解的误差
    // 使用几个vector存储x_1 ~ x_99的各个数据，y_0和y_n是准确解，所以不再考虑
    std::vector<double> x(n-1, 0.0);
    std::vector<double> analytical_y(n-1, 0.0);
    std::vector<double> error_gauss(n-1, 0.0);
    std::vector<double> error_seidel(n-1, 0.0);
    for (int i = 0; i < n - 1; ++i) {
        x[i] = (i + 1) * h;
        analytical_y[i] = analyticalSolution(x[i], epsilon, a);
        error_gauss[i] = std::abs(solution_gauss[i] - analytical_y[i]);
        error_seidel[i] = std::abs(solution_seidel[i] - analytical_y[i]);
        std::cout<<"x:"<<x[i]<<"\tanalytical_y"<<analytical_y[i]<<"\tgauss y: "<<solution_gauss[i]<<"\tseidel: "<<solution_seidel[i]<<std::endl;
    }
    // 打开data.txt写入数据
    std::ofstream file("data.txt");
    for(int i = 0; i < n - 1; i++)
    {
        file << x[i] << " " << analytical_y[i] << " " << solution_gauss[i] << " " << solution_seidel[i] << " " << error_gauss[i] << " " << error_seidel[i] << std::endl;
    }
    file.close();

    return 0;
}


std::vector<double> gaussElimination(std::vector<std::vector<double>>& matrix)
{
    /* 
        这里实现Gauss列主元消元法
        不需要记录顺序，计算完再换回来，这是我最开始的误区，但是代码不删除，选择注释，之后时刻提醒自己
    */ 
    int rows = matrix.size();
    int cols = matrix[0].size() - 1; // 最后一列为向量b
    // 我们做的A都是方阵
    // std::vector<int> swap_rec(rows);
    // for(int i = 0; i < rows; i++)
    // {
    //     swap_rec[i] = i;
    //     /*
    //         专门开一个数组用于存储后面交换的两行
    //         初始化每个索引和对应值相等，后续增广矩阵交换两行时，该数组对应着更换对应位置的值，用于记录
    //         保证解出向量后对应元素可以换回来
    //     */
    // }

    // 化简为上三角形矩阵
    for (int i = 0; i < rows - 1; i++) {
        // 找到当前列绝对值最大的行，并将其交换至当前行
        int maxRow = i;
        double maxVal = std::fabs(matrix[i][i]);
        for (int k = i + 1; k < rows; ++k) {
            if (std::fabs(matrix[k][i]) > maxVal) {
                maxVal = std::fabs(matrix[k][i]);
                maxRow = k;
            }
        }
        // 将主元换到对角位置(若非对角元素)，同时将swap_rec中两对应位置元素互换
        if (maxRow != i) {
            std::swap(matrix[i], matrix[maxRow]);
            // std::swap(swap_rec[i], swap_rec[maxRow]);
        }

        // 消元
        for (int j = i + 1; j < rows; j++) {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k < cols + 1; k++) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }

    // 回代
    std::vector<double> tmp_solution(rows);
    for (int i = rows - 1; i >= 0; i--) {
        tmp_solution[i] = matrix[i][cols];
        for (int j = i + 1; j < cols; j++) {
            tmp_solution[i] -= matrix[i][j] * tmp_solution[j];
        }
        tmp_solution[i] /= matrix[i][i];
    }
    return tmp_solution;
    // // 恢复顺序
    // std::vector<double> solution(rows);
    // // 如果根据rec数组记录的更换顺序进行更换这个代码设计也太复杂了，鉴于只是一个向量，我们直接新开一个solution作解
    // for(int i = 0; i < rows; i++)
    // {
    //     solution[swap_rec[i]] = tmp_solution[i];
    // }
    // return solution;
}


double norm(const std::vector<double>& x)
{
    // 求解欧几里德范数
    double result = 0.0;
    for (double val : x)
    {
        result += val * val;
    }
    return sqrt(result);
}


// 这里实现Gauss-Seidel迭代法求解线性方程组
std::vector<double> gaussSeidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double epsilon, int maxIterations)
{
    // A是一个方阵
    int n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    for(int it = 0; it < maxIterations; it++)
    {
        for(int i = 0; i < n; i++)
        {
            double sigma = 0.0; // sigma用于存储暂时的计算结果
            for (int j = 0; j < i; j++) 
            {
                sigma += A[i][j] * x_new[j];
            }
            for (int j = i+1; j < n; j++)
            {
                sigma += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sigma) / A[i][i];
        }

        // 检查收敛条件啦
        std::vector<double> residual(n, 0.0);
        for(int i = 0; i < n; i++)
        {
            residual[i] = x_new[i] - x[i];
        }
        // 用2范数判断
        double error = norm(residual);
        if(error < epsilon)
        {
            std::cout<<"converged after"<< it + 1<<"iterations"<<std::endl;
            return x_new;
        }
        x = x_new;
    }

    // 制定迭代次数还没收敛到指定的精度，简单进行一个报错
    std::cerr << "Did not converge after " << maxIterations << " iterations." << std::endl;
    return x;
}

void printMatrix(const std::vector<std::vector<double>>& matrix)
{
    // 写一个打印矩阵的函数，万一用上了呢
    int rows = matrix.size();
    int cols = matrix[0].size();

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < cols; j++)
        {
            std::cout<<matrix[i][j]<<"\t";
        }
        std::cout<<std::endl;
    }
}

void test()
{
    // 测试环节
    std::vector<std::vector<double>> matrix = {
        {2, 1, 1, 8},
        {1, -1, 0, 3},
        {3, -1, 2, 11}
    };

    std::vector<double> solution = gaussElimination(matrix);

    std::cout << "Solution:" << std::endl;

    for (double x : solution) {
        std::cout << x << " ";
    }
    std::cout << std::endl;



    // 这里测试Gauss-Seidel
    std::vector<std::vector<double>> A = {
        {2, 1, 1},
        {1, -1, 0},
        {3, -1, 2}
    };

    std::vector<double> b = {8, 3, 11};

    std::vector<double> sol2 = gaussSeidel(A, b, 1e-5, 50);

    std::cout << "Solution:" << std::endl;

    for (double x : sol2) {
        std::cout << x << " ";
    }
    std::cout << std::endl;
}
