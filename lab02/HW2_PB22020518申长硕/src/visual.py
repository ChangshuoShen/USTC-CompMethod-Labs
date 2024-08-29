'''
这是一个简单的实现可视化的程序，接受cpp传来的数值，对其进行可视化
同时将error进行可视化
'''
import numpy as np
import matplotlib.pyplot as plt

def data_visualization():
    with open("./data.txt", "r") as file:
        lines = file.readlines()
        data = [line.strip().split(' ') for line in lines]
    data_array = np.array(data, dtype=float)

    # 打印结果
    print("Data as NumPy array:")
    print(data_array)

    fig, axes = plt.subplots(2, 1, figsize=(6, 6))
    
    axes[0].plot(data_array[:, 0], data_array[:, 1], color='black', label='analytical ans')
    axes[0].plot(data_array[:, 0], data_array[:, 2], color='blue', label='gauss ans')
    axes[0].plot(data_array[:, 0], data_array[:, 3], color='red', label='seidel ans')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].set_title('epsilon = 0.001')
    axes[0].legend()

    axes[1].plot(data_array[:, 0], data_array[:, 4], color='blue', label='gauss abs error')
    axes[1].plot(data_array[:, 0], data_array[:, 5], color='red', label='seidel abs error')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('error')
    axes[1].set_title('error of different methods')
    axes[1].legend()

    plt.savefig('visualization.png')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    data_visualization()