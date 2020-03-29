import matplotlib.pyplot as plt
import numpy as np
import os

if os.name == 'posix':
    path = "/home/ben/Documents/uni/760/assignment1/output/"


def main():
    bestObj = np.genfromtxt(path + "bestObj")
    allObj = np.genfromtxt(path + "allObj") / 53843.4
    bestObjIdx = np.genfromtxt(path + "bestObjIdx")
    idx = np.genfromtxt(path + "indices", dtype=int)

    _, ax = plt.subplots(1)
    ax.set_yscale("log")
    ax.scatter(np.arange(allObj.shape[0]), allObj, marker='x', color='k', s=0.1)
    ax.plot(bestObjIdx[0:idx[0]], bestObj[0:idx[0]])
    ax.plot(bestObjIdx[idx[0]:idx[1]], bestObj[idx[0]:idx[1]])
    ax.plot(bestObjIdx[idx[1]:], bestObj[idx[1]:])
    ax.set_ylim(3e-5, 10)

    ax.set_title("Next Descent with 3 Random Restarts")
    ax.set_ylabel("Objective")
    ax.set_xlabel("Objective function Evaluations")

    plt.show()

if __name__ == "__main__":
    main()