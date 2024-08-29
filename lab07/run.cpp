#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

double ax(double t) {return std::sin(t) / (std::sqrt(t) + 1);}

double ay(double t) {return std::log(t + 1) / (t + 1);}

double Romberg(double a, double b, double e, int M, std::function<double(double)> f, int &hit_count, int &total_count) {
    std::vector<std::vector<double>> R(M, std::vector<double>(M, 0.0));
    R[0][0] = (f(a) + f(b)) * (b - a) / 2;
    auto h = [a, b](int x) { return (b - a) / std::pow(2, x - 1); };
    // 每次进入循环，记录总调用次数
    total_count++;
    int result_index = 0;
    for (int k = 1; k < M; k++) 
    {
        result_index = k;
        double sum = 0.0;
        for (int i = 1; i <= std::pow(2, k - 1); i++) {
            sum += f(a + (2 * i - 1) * h(k + 1));
        }
        R[k][0] = (R[k - 1][0] + h(k) * sum) / 2;

        for (int j = 1; j <= k; j++) {
            R[k][j] = R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1]) / (std::pow(4, j) - 1);
        }

        if (std::abs(R[k][k] - R[k - 1][k - 1]) < e) {
            hit_count++;
            break;
        }
    }

    return R[result_index][result_index];
}

double vx(double t, double a, double e, int M, std::function<double(double)> ax, int &hit_count, int &total_count) {
    return Romberg(a, t, e, M, ax, hit_count, total_count);
}

double vy(double t, double a, double e, int M, std::function<double(double)> ay, int &hit_count, int &total_count) {
    return Romberg(a, t, e, M, ay, hit_count, total_count);
}

void run_experiment(int M, double e) {
    std::ofstream file("points_M" + std::to_string(M) + ".txt");
    std::vector<double> time_area(100);
    for (int i = 0; i < 100; i++) {
        time_area[i] = (i + 1) / 10.0;
    }

    int hit_count = 0;
    int total_count = 0;
    double a = 0;
    std::vector<double> x;
    std::vector<double> y;

    for (double t : time_area) {
        double b = t;
        auto vx_lambda = [a, e, M, &hit_count, &total_count](double s) { 
            return vx(s, a, e, M, ax, hit_count, total_count); 
        };

        auto vy_lambda = [a, e, M, &hit_count, &total_count](double s) { 
            return vy(s, a, e, M, ay, hit_count, total_count); 
        };

        x.push_back(Romberg(a, b, e, M, vx_lambda, hit_count, total_count));
        y.push_back(Romberg(a, b, e, M, vy_lambda, hit_count, total_count));
    }

    for (size_t i = 0; i < x.size(); i++) {
        file << x[i] << " " << y[i] << std::endl;
    }

    file.close();
    double ratio = static_cast<double>(hit_count) / total_count;
    std::cout << "M = " << M << ", 达到误差要求的比例: " << ratio * 100 << "%" << std::endl;
}

int main() {
    double e = 1e-6;
    std::vector<int> M_values = {4, 8, 12, 16, 20};

    for (int M : M_values) {
        run_experiment(M, e);
    }

    return 0;
}