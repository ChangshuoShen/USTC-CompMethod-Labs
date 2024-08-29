import numpy as np
from math import sin, log, sqrt
import matplotlib.pyplot as plt

# 定义 Romberg 算法
def Romberg(a, b, e, M, f):
    R = np.zeros((M, M))
    global hit_count, total_count  # 使用全局变量
    R[0][0] = (f(a) + f(b)) * (b - a) / 2
    h = lambda x: (b - a) / 2 ** (x - 1)
    total_count += 1  # 每次进入循环,总调用次数加 1
    for k in range(1, M):
        result_index = k
        R[k][0] = (R[k - 1][0] + h(k) * sum([f(a + (2 * i - 1) * h(k + 1)) for i in range(1, 2 ** (k - 1) + 1)])) / 2
        for j in range(1, k + 1):
            R[k][j] = R[k][j - 1] + (R[k][j - 1] - R[k - 1][j - 1]) / (4 ** j - 1)
        if abs(R[k][k] - R[k - 1][k - 1]) < e:
            hit_count += 1  # 达到误差要求,命中次数加 1
            break
    
    return R[result_index][result_index]

# 定义加速度函数 ax(t), ay(t)
def ax(t):
    return sin(t) / (sqrt(t) + 1)

def ay(t):
    return log(t + 1) / (t + 1)

# 对速度应用 Romberg 算法
time_area = [i/10 for i in range(1, 101)]
n, M, e = 1, 12, 10e-6
x = []
y = []
hit_count = 0
total_count = 0

for t in time_area:
    a, b = 0, t
    vx = lambda s: Romberg(a, s, e, M, ax)
    vy = lambda s: Romberg(a, s, e, M, ay)
    x.append(Romberg(a, b, e, M, vx))
    y.append(Romberg(a, b, e, M, vy))

ratio = hit_count / total_count
print(f"达到误差要求的比例: {ratio * 100:.2f}%")

plt.plot(x, y, marker='o', linestyle='-', color='r', label='Curve', markersize=2)
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Trajectory')
plt.legend()
plt.grid(True)

# Display the plot
plt.show()