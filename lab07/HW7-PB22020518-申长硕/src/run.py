import matplotlib.pyplot as plt

def plot_results(M_values):
    fig, axs = plt.subplots(1, len(M_values), figsize=(len(M_values) * 6, 15))
    fig.suptitle('Trajectory for different M values')
    
    for i, M in enumerate(M_values):
        filename = f'points_M{M}.txt'
        x = []
        y = []
        with open(filename, 'r') as file:
            for line in file:
                point = line.split()
                x.append(float(point[0]))
                y.append(float(point[1]))
        
        axs[i].plot(x, y, marker='o', linestyle='-', markersize=2, label=f'M = {M}')
        axs[i].set_xlabel('X-axis')
        axs[i].set_ylabel('Y-axis')
        axs[i].legend()
        axs[i].grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Adjust layout to make room for the title
    plt.show()


if __name__ == "__main__":
    
    M_values = [4, 8, 12, 16, 20]
    plot_results(M_values)
