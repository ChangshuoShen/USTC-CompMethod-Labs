import numpy as np
import matplotlib.pyplot as plt


def read_file(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()
    data = [[float(num) for num in line.replace(',', ' ').strip().split()] for line in lines]
    data_array = np.array(data)
    # print(data_array)
    return data_array

def set_axes_limits(ax_list):
    # 找到四个子图中 x 和 y 轴的范围的最大值和最小值
    x_min = min(ax.get_xlim()[0] for ax in ax_list)
    x_max = max(ax.get_xlim()[1] for ax in ax_list)
    y_min = min(ax.get_ylim()[0] for ax in ax_list)
    y_max = max(ax.get_ylim()[1] for ax in ax_list)
    min_ = min(x_min, y_min)
    max_ = max(x_max, y_max)
    # 设置所有子图的 x 和 y 轴范围
    for ax in ax_list:
        ax.set_xlim(min_, max_)
        ax.set_ylim(min_, max_)


if __name__ == '__main__':
    data = read_file('iris.txt')
    transformed_data = read_file('transformed_data.txt')
    # 提取最后一列作为颜色分类
    colors = data[:, -1]
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    # 绘制散点图，并根据最后一列的值选择颜色
    axs[0, 0].scatter(data[:, 0], data[:, 1], c=colors, cmap=plt.cm.get_cmap('viridis', 3), s=8)
    axs[0, 1].scatter(data[:, 2], data[:, 3], c=colors, cmap=plt.cm.get_cmap('viridis', 3), s=8)
    axs[1, 0].scatter(transformed_data[:, 0], transformed_data[:, 1], c=colors, cmap=plt.cm.get_cmap('viridis', 3), s=8)
    axs[1, 1].scatter(transformed_data[:, 2], transformed_data[:, 3], c=colors, cmap=plt.cm.get_cmap('viridis', 3), s=8)
    # 获取所有的子图对象列表
    all_axes = fig.axes
    # 统一设置所有子图的 x 和 y 轴范围
    set_axes_limits(all_axes)
    plt.tight_layout()
    plt.show()

