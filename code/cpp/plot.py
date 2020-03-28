import matplotlib.pyplot as plt
import numpy as np
import os

if os.name == 'posix':
    path = "/home/ben/Documents/uni/760/assignment1/output/"


def main():
    obj = np.genfromtxt(path + "objTest")

    _, ax = plt.subplots(1)
    ax.set_yscale("log")
    ax.plot(obj)
    plt.show()

if __name__ == "__main__":
    main()