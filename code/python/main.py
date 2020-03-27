import numpy as np
import os

if(os.name == 'posix'):
    path = "/home/ben/Documents/uni/760/assignment1/data/"

def runTests():
    pass

def swap(arr, a, b):
    s = np.copy(arr)
    s[a], s[b] = np.copy(s[b]), np.copy(s[a])
    return s

def objective(dx, dy):
    return 5 * abs(dy) + abs(dx)

def readData(problemFileName = "ProbA.txt", n=120):
    # path = "/home/ben/Documents/uni/760/assignment1/data/"
    postionsFileName = "Positions.txt"

    weights = np.insert(np.genfromtxt(path+problemFileName, skip_header=1), 0, 0) # prepend 0
    
    # weights = np.pad(weights, (0,n-len(weights)), 'constant', constant_values=(0,0))
    s, x, y = np. genfromtxt(path+postionsFileName, skip_header=1).T

    s[weights.shape[0]-1:] = 0
    s = np.array(s, dtype=int)

    return s, weights, x, y

def computeDxDy(s, weights, x, y, massTotal):
    dx = 0
    dy = 0

    for i, sI in enumerate(s):
        # x of position, weight of container
        dx += x[i] * weights[sI]
        dy += y[i] * weights[sI]

    return dx / massTotal, dy / massTotal

def computeObjective(s, weights, x, y, massTotal):
    dx, dy = computeDxDy(s, weights, x, y, massTotal)
    return objective(dx, dy)

def nextDescentIteration(s, weights, x, y, massTotal, obj):
    for i in range(s.shape[0]):
        for j in range(i):
            if((s[i] == 0 and s[j] == 0) or (j == i + 60)): continue
            # compute obj of swapping containers i and j. Swap if obj (i, j) less than current obj
            tempS = swap(s, i, j)
            nextObj = computeObjective(tempS, weights, x, y, massTotal)
            if(nextObj < obj): 
                return tempS, nextObj, 0

    return None, None, 1

def runNextDescent(s, weights, x, y, massTotal, obj):
    its = 0

    try:
        while(True):
            s, obj, status = nextDescentIteration(s, weights, x, y, massTotal, obj)
            if(status == 1): break
            print("Current obj: ", obj)
            its += 1

    except(KeyboardInterrupt):
        pass

    return s, obj, its

def saveSolution(fname, solution):
    np.savetxt("/home/ben/Documents/uni/760/assignment1")


def main():
    s, weights, x, y = readData()
    # print(s)
    # print(s.shape)
    # print(weights.shape)

    massTotal = np.sum(weights)

    if(True):
        np.random.seed(0)
        np.random.shuffle(s) # works in place for some reason

        obj = computeObjective(s, weights, x, y, massTotal)

        s, obj, iterations = runNextDescent(s, weights, x, y, massTotal, obj)

        # print(s)
        print("Objective", obj, "found in", iterations, "iterations")

if __name__ == "__main__":
    main()
