import numpy as np

data = []

with open('./point.txt', encoding='utf8') as file:
    for line in file:
        line = line.strip()
        values = line.split('\t')
        point = [float(values[0]), float(values[1])]
        data.append(point)
print('data', data)


n = len(data)

h = [data[i + 1][0] - data[i][0] for i in range(n-1)]
lambda_li = [(h[i+1] / (h[i] + h[i + 1])) for i in range(n-2)]
miu = [1 - lambda_li[i] for i in range(n-2)]
d =[(6 / (h[i] + h[i - 1]) * ((data[i + 1][1] - data[i][1]) / h[i] - (data[i][1] - data[i - 1][1]) / h[i - 1])) for i in range(1, n-1)]

def thomas_algorithm(a, b, c, d):
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

a = [2.0 for _ in range(n-2)]
b = lambda_li[0:n-3]
c = miu[1:n-2]
d = d[0:n-1]

M = thomas_algorithm(a, b, c, d)
M = [0.0] + list(M) + [0.0]

formatted_M = [round(value, 2) for value in M]

print(formatted_M)