import numpy as np
import matplotlib.pyplot as plt
import copy

with open('./point.txt', 'r') as file:
    lines = file.readlines()

    # 将数据转换为NumPy数组
    data = np.zeros((len(lines), 2))
    for i, line in enumerate(lines):
        elements = line.strip().split()
        data[i, 0] = float(elements[0])
        data[i, 1] = float(elements[1])
    # data = np.array([elements for line in lines for elements in line.strip().split()], dtype=float)
print(data)

class Spline:
    def __init__(self, points):
        self.points = points
        self.size = len(self.points)
        self.h = [self.points[i + 1][0] - self.points[i][0] for i in range(self.size-1)]
        self.lambda_li = [(self.h[i+1] / (self.h[i] + self.h[i + 1])) for i in range(self.size-2)]
        self.miu = [1 - self.lambda_li[i] for i in range(self.size-2)]
        self.d =[(6 / (self.h[i] + self.h[i - 1]) * ((self.points[i + 1][1] - self.points[i][1]) / self.h[i] - (self.points[i][1] - self.points[i - 1][1]) / self.h[i - 1])) for i in range(1, self.size-1)]
        self.func_points = [] # 这个用于存储在各个点上的坐标

        
    def thomas_algorithm(self):
        a = [2.0 for _ in range(self.size - 2)]
        b = self.lambda_li[0:self.size-3]
        c = self.miu[1:self.size-2]
        d = self.d[0:self.size-1]
        n = len(d)
        
        c_star = np.zeros(n-1)
        d_star = np.zeros(n)
        
        c_star[0] = c[0] / a[0]
        d_star[0] = d[0] / a[0]
        
        for i in range(1, n):
            m = 1 / (a[i] - b[i-1] * c_star[i-1])
            if i < n - 1:
                c_star[i] = c[i] * m
            d_star[i] = (d[i] - d_star[i-1] * b[i-1]) * m
        
        x = np.zeros(n)
        x[-1] = d_star[-1]
        
        for i in range(n - 2, -1, -1):
            x[i] = d_star[i] - c_star[i] * x[i + 1]
        
        return x
    
    def calc_M(self):
        self.M = [0] + list(self.thomas_algorithm()) + [0]
        self.formatted_M = [round(val, 2) for val in self.M]
    
    def special_func_for_spline(self, i):
        # 这个函数根据index i返回一个区间内的1000个函数及函数值，依据上面的根据M表达的插值公式
        # 注意 i 表示的是第i + 1个
        
        x_start = self.points[i, 0]
        x_end = self.points[i + 1, 0]
        h = self.h[i]
        y1 = self.points[i, 1]
        y2 = self.points[i + 1, 1]
        m1 = self.M[i]
        m2 = self.M[i + 1]
        return_points = np.linspace(x_start, x_end, 1000, endpoint=False)
        # print(x_start, return_points)
        def func(m1, m2, y1, y2, x_start, x_end, h, x):
            # 这是一个使用M_i, M_i+1两个参数和x_start, x_end计算对应点x函数指的函数
            term1 = ((x_end - x)**3 * m1 + (x - x_start)**3 * m2) / (6 * h)
            term2 = (y1 * (x_end - x) + y2 * (x - x_start)) / h
            term3 = (h * ((m1 * (x_end - x)) + (m2 * (x - x_start)))) / 6
            return term1 + term2 - term3
        
        return_vals = np.array([func(m1, m2, y1, y2, x_start, x_end, h, x) for x in return_points])
        return np.column_stack((return_points, return_vals))

    
    def solve_and_plot(self):
        self.calc_M()
        
        for i in range(self.size - 1):
            points_this_set = self.special_func_for_spline(i)
            self.func_points.extend(points_this_set.tolist())
            plt.plot(points_this_set[:, 0], points_this_set[:, 1])
            
        plt.xlabel('x')
        plt.ylabel('S(x)')
        plt.title('Spline?')
        plt.grid(True)
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plt.show()
        
    def __sub__(self, other):
        self_points_array = np.array(self.func_points)
        other_points_array = np.array(other.func_points)
        
        sub_points = self_points_array[:, 1] - other_points_array[:, 1]  # 计算第二列的差值
        x_values = self_points_array[:, 0]  # x 值取自 func_points 的第一列
        
        plt.plot(x_values, sub_points)
        plt.xlabel('x')
        plt.ylabel('Sub(x)')
        plt.title('Error')
        plt.grid(True)
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plt.show()
        
        sub_points = np.column_stack((x_values, sub_points))
        for i in range(0, len(sub_points), 100):
            print(sub_points[i])
        
        

sample1 = Spline(data)
sample1.solve_and_plot()
print(sample1.M)
# sample1.func_points

data1 = copy.deepcopy(data)
data1[9, 1] = 10.
# print(data1)
sample2 = Spline(data1)
sample2.solve_and_plot()
print(sample2.M)

sample1 - sample2