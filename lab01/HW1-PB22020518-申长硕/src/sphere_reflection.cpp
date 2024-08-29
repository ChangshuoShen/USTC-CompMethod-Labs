#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<cmath>


int sgn(double x);
double* coefficient(double uvw[3]);
double func(double x, double coe[5]);
double bisection(double left, double right, double epsilon, double coe[5]);
double calc_beta(double x, double uvw[3]);
double *calc_uvw(double x_p, double x_q, double y_q);
double newton(double x, double epsilon, double coe[5]);

// 计算Q点对称点Q'
// double* calc_reflection(double x_t, double y_t, double x_q, double y_q);
double* calc_reflection(double x_t, double y_t, double x_p, double x_q, double y_q);



int main()
{
    // 不通过主函数传参，打开input.txt文件传入一系列x_p, x_q, y_q
    std::ifstream file("./input.txt");
    std::vector<std::vector<double>> data;
    std::vector<std::vector<double>> ans;

    if(file.is_open())
    {
        std::string line;
        while(std::getline(file, line))
        {
            // 逐行读取嘛
            std::vector<double> row;
            std::stringstream ss(line);
            double value;

            for(int i = 0; i < 3; i++)
            {
                if(ss >> value)
                {
                    row.push_back(value);
                    // std::cout << value << "\t";
                }
                else{
                    std::cerr << "Expected 3 numbers for each line" << std::endl;
                    return 1;
                }
            }
            data.push_back(row);
        }
    }
    else
    {
        std::cerr << "Unable to open the file" << std::endl;
        return 1;
    }
    // 现在我们已经将数据存入<vector> data,
    for(const auto& row: data)
    {
        // 输出检查一番
        // for(const auto& value : row)
        // {
        //     std::cout << value <<"\t";
        // }
        // std::cout << std::endl;
        // 先计算u, v, w
        double* uvw = calc_uvw(row[0], row[1], row[2]);

        // 之后计算该四次函数的系数：
        double* coe = coefficient(uvw);

        // 计算coefficients之后使用二分法计算方程的根

        // 从0, |1 / x_p|开始迭代，因为x_t在-1～0之间，其余无意义
        double alpha = bisection(0, -1.0 / row[0], 1e-8, coe);

        // 这里简单尝试一下牛顿法
        // double alpha = newton(-1.0 / row[0], 1e-10, coe);


        // 之后用alpha带入公式计算beta
        double beta = calc_beta(alpha, uvw);

        double x_t = alpha * row[0] + beta * row[1];
        double y_t = beta * row[2];

        // 这里计算Q',q_reflected
        double* q_r = calc_reflection(x_t, y_t, row[0], row[1], row[2]); 
        
        std::vector single_ans = {x_t, y_t, q_r[0], q_r[1]};
        ans.push_back(single_ans);
        // std::cout << "(" << x_t << ", " << y_t << ")"<< std::endl;
        // 简单进行一个内存的释放
        delete[] uvw;
        delete[] coe;
        delete[] q_r;
    }

    // 输出检查
    std::cout<<"----------answer summary:-----------"<<std::endl;
    for(const auto& single_ans : ans)
    {
        std::cout<<"T: ("<<single_ans[0] <<", "<<single_ans[1]<<")"<<"Q': ("<<single_ans[2]<<", "<<single_ans[3]<<")"<<std::endl;
    }
    return 0;

}


int sgn(double x)
{
    if(x >= 0){
        return 1;
    }
    else{
        return -1;
    }
}

double* calc_uvw(double x_p, double x_q, double y_q)
{
    // 用P,Q坐标计算出需要解的x的多项式系数
    // 按报告中公式求出u, v, w
    double* uvw = new double[3];
    uvw[0] = x_p * x_p;
    uvw[1] = x_p * x_q;
    uvw[2] = x_q * x_q + y_q * y_q;
    // std::cout<<"u, w, v:"<<uvw[0]<<"\t"<<uvw[1]<<"\t"<<uvw[2]<<std::endl;
    return uvw;
}


double* coefficient(double uvw[3])
{
   // 后利用u, v, w计算出多项式系数
   /*
        计算公式：
        4u(uw- v^2) x^4 
        + 4(v^2 - uw)x^3 
        + (u + 2v + w - 4uw)x^2 
        + 2(w - v)x 
        + (w - 1) 
        = 0
        依此计算各个系数
   */
    double* coe = new double[5];
    coe[0] = 4 * uvw[0] * (uvw[0] * uvw[2] - uvw[1] * uvw[1]);
    coe[1] = 4 * (uvw[1] * uvw[1] - uvw[0] * uvw[2]);
    coe[2] = uvw[0] + 2 * uvw[1] + uvw[2] - 4 * uvw[0] * uvw[2];
    coe[3] = 2 * (uvw[2] - uvw[1]);
    coe[4] = uvw[2] - 1;
    return coe;
}


double func(double x, double coe[5])
{
    // 使用coe传入五个系数，提前算出来防止重复计算影响效率
    return (coe[0] * std::pow(x, 4)
        + coe[1] * std::pow(x, 3)
        + coe[2] * std::pow(x, 2)
        + coe[3] * x
        + coe[4]
    );
}


// 二分法
double bisection(double left, double right, double epsilon, double coe[5])
{
    std::cout<<std::endl<<"~~~bisection begin~~~"<<std::endl;
    // std::cout<<"coe数组"<<std::endl;
    // for(int i = 0; i < 5; i++)
    // {
    //     std::cout<<coe[i]<<std::endl;
    // }
    double mid;
    int run_times = 0;
    while(fabs(right - left) > epsilon)
    {
        mid = (left + right) / 2.0;
        double fun_value = func(mid, coe);
        // std::cout<<"value this loop"<<"\t mid:"<<mid<<"func_value"<<fun_value<<std::endl;
        // std::cout<<"f(left)"<<func(left, coe)<<std::endl;
        // std::cout<<"f(right)"<<func(right, coe)<<std::endl;
        if(fabs(fun_value) < epsilon)
        {
            return mid;
        }
        else if(sgn(fun_value) * sgn(func(left, coe)) < 0)
        {
            right = mid;
        }
        else
        {
            left = mid;
        }
        run_times++;
    }
    std::cout<<"running times: "<<run_times<<std::endl;
    return (left + right) / 2.0;
}


// 简单进行一个导数的求
double deriv(double x, double coe[5])
{
    // derivative
    return (4 * coe[0] * std::pow(x, 3)
        + 3 * coe[1] * std::pow(x, 2)
        + 2 * coe[2] * x
        + coe[3]
    );
}


/*
    使用Newton迭代法：
    迭代格式：x_{n+1} = x_n - f(x_n) \over f'(x_n)
*/ 
double newton(double x, double epsilon, double coe[5])
{
    
    // 使用牛顿迭代
    double f = func(x, coe);
    while(fabs(f) > epsilon)
    {
        x = x - func(x, coe) / deriv(x, coe);
        std::cout<<"x:"<<x<<std::endl;
        f = func(x, coe);
        std::cout<<"f(x):"<<f;
    }
    std::cout<<"牛顿迭代结果:"<<x<<std::endl;

    return x;
}


// 得到alpha之后按照公式beta
double calc_beta(double x, double uvw[3])
{
    // y = {{1 + x - 2ux^2} \over {1 + 2vx}}
    double beta = (1.0 + x - 2.0 * uvw[0] * x * x) / (1.0 + 2.0 * uvw[1] * x);
    // std::cout<<"beta:"<<beta<<std::endl;
    return beta;
}

// double* calc_reflection(double x_t, double y_t, double x_q, double y_q)
// {
//     // 计算Q点的对称点Q'
//     /*
//         Q  = sT^{\perp} + tT
//         Q' = sT^{\perp} - tT
//         其中：
//         t = \vec{TQ} \cdot \vec{OT}\\
//         s = \vec{TQ} \cdot \vec{T^{\perp}}

//         根据T(x_t, y_t)得到T_perp是(y_t, -x_t)
//     */
//     double x_q_tmp = x_q - x_t;
//     double y_q_tmp = y_q - y_t;

//     double t = (x_q_tmp - x_t) * x_t + (y_q_tmp - y_t) * y_t;
//     double s = (x_q_tmp - y_t) * y_t - (y_q_tmp + x_t) * x_t;

//     // Q' = sT^{\perp} - tT
//     double x_qr = s * y_t - t * y_t + x_t;
//     double y_qr = s * y_t + t * x_t + y_t;
//     std::cout<<"Q\': ("<<x_qr<<", "<<y_qr<<")"<<std::endl;
//     double* coo_q_r = new double[2];
//     coo_q_r[0] = x_qr;
//     coo_q_r[1] = y_qr;
//     return coo_q_r;
// }

double* calc_reflection(double x_t, double y_t, double x_p, double x_q, double y_q)
{
    // lambda = |TQ| / |PT|
    double lambda = std::pow(((x_t - x_q) * (x_t - x_q) + (y_t - y_q) * (y_t - y_q)), 0.5) / std::pow(((x_t - x_p) * (x_t - x_p) + y_t * y_t), 0.5);
    // std::cout<<"lambda:"<<lambda<<std::endl;

    double x_qr = x_t + lambda * (x_t - x_p);
    double y_qr = y_t + lambda * y_t;   // y_p = 0

    double* coordinates = new double[2];
    coordinates[0] = x_qr;
    coordinates[1] = y_qr;

    // std::cout<<"Q': ("<<x_qr<<", "<<y_qr<<")"<<std::endl;

    return coordinates;
}

