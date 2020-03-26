import numpy as np
import os

if(os.name == 'posix'):
    path = "/home/ben/Documents/uni/760/assignment1/data/"

def runTests():
    pass

def swap(arr, a, b):
    # assert(len(s)==120, "s must be 120 but is {}".format(len(s)))
    assert(a < b, "a must be less than b, a={}, b={}".format(a, b))
    s = np.copy(arr)
    s[a], s[b] = s[b], s[a]
    return s

def objective(dx, dy):
    return 5 * abs(dy) + abs(dx)

def readData(problemFileName = "ProbA.txt", n=120):
    # path = "/home/ben/Documents/uni/760/assignment1/data/"
    postionsFileName = "Positions.txt"

    containers = np.genfromtxt(path+problemFileName, skip_header=1)
    containers = np.pad(containers, (0,n-len(containers)), 'constant', constant_values=(0,0))
    _, x, y = np. genfromtxt(path+postionsFileName, skip_header=1).T
    containers = np.vstack((np.arange(120), containers)).T

    return containers, x, y

def computeDxDy(containers, x, y, massTotal):
    dx = 0
    dy = 0

    for container in containers:
        i = int(container[0])
        dx += x[i] * container[1]
        dy += y[i] * container[1]

    return dx / massTotal, dy / massTotal

def computeObjective(containers, x, y, massTotal):
    dx, dy = computeDxDy(containers, x, y, massTotal)
    return objective(dx, dy)

def nextDescentIteration(containers, x, y, massTotal, obj):
    for i in range(containers.shape[0]):
        for j in range(i):
            if((containers[i][1] == 0 and containers[j][1] == 0) or (j == i + 60)): continue
            # compute obj of swapping containers i and j. Swap if obj (i, j) less than current obj
            temp = swap(containers, i, j)
            nextObj = computeObjective(temp, x, y, massTotal)

            if(nextObj < obj): 
                containers = temp
                return containers, nextObj, 0

    return None, None, 1

def runNextDescent(containers, x, y, massTotal, obj):
    its = 0

    try:
        while(True):
            tempContainers, tempObj, status = nextDescentIteration(containers, x, y, massTotal, obj)
            if(status == 1): break
            containers = tempContainers
            obj = tempObj
            print("Current obj: ", obj)
            its += 1

    except(KeyboardInterrupt):
        pass

    return containers, obj, its

def saveSolution(fname, solution):
    np.savetxt("/home/ben/Documents/uni/760/assignment1")


def main():
    containers, x, y = readData()
    
    np.random.seed(0)
    np.random.shuffle(containers) # works in place for some reason

    massTotal = np.sum(containers, axis = 0)[1]
    obj = computeObjective(containers, x, y, massTotal)

    solution, slnObj, iterations = runNextDescent(containers, x, y, massTotal, obj)

    print(containers)
    print("Objective", obj, "found in", iterations, "iterations")

if __name__ == "__main__":
    main()
