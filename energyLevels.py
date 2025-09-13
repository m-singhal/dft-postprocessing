import numpy as np 
import matplotlib.pyplot as plt

def energyLevels(Emin = None, Emax = None):
    with open("EIGENVAL", "r") as f:
        lines = f.readlines()
        n = int(lines[5].strip().split()[-1])
        data = np.loadtxt(lines[8:8+n])
        data = data[:,1:]
        spin = data.shape[1]//2
    
    if spin == 1:
        for y in range(n):
            if int(data[y,1]) == 1:
                plt.hlines(data[y,0], xmin=-0.5, xmax=0.5, color="red", linewidth=2, alpha=0.5)
            else:
                plt.hlines(data[y,0], xmin=-0.5, xmax=0.5, color="black", linewidth=2,  alpha=0.5)
        plt.ylabel("Energy Levels (eV)")
        plt.xticks([])
        plt.ylim(Emin, Emax)
        plt.xlim(-1,1)

    else:
        fig, axes = plt.subplots(1, 2, figsize=(10, 6), tight_layout=True)
        for y in range(n):
            if int(data[y,2]) == 1:
                axes[0].hlines(data[y,0], xmin=-0.5, xmax=0.5, color="red", linewidth=2, alpha=0.5)
            elif int(data[y,2]) == 0:
                axes[0].hlines(data[y,0], xmin=-0.5, xmax=0.5, color="black", linewidth=2, alpha=0.5)
            
            if int(data[y,3]) == 1:
                axes[1].hlines(data[y,1], xmin=-0.5, xmax=0.5, color="red", linewidth=2, alpha=0.5)
            elif int(data[y,3]) == 0:
                axes[1].hlines(data[y,1], xmin=-0.5, xmax=0.5, color="black", linewidth=2, alpha=0.5)
        axes[0].set_title("Spin Up")
        axes[1].set_title("Spin Down")
        axes[0].set_ylabel("Energy Levels (eV)")
        axes[0].set_xticks([])
        axes[1].set_xticks([])
        axes[0].set_ylim(Emin, Emax)
        axes[1].set_ylim(Emin, Emax)
        axes[0].set_xlim(-1,1)
        axes[1].set_xlim(-1,1)
        
    
    
    
    plt.show()

energyLevels(-3,7)
