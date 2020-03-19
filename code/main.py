import numpy as np

def runTests():
    pass

def swap(s, a, b):
    assert(len(s)==120, "s must be 120 but is {}".format(len(s)))
    assert(a<b, "a must be less than b, a={}, b={}".format(a, b))
    s[a], s[b] = s[b], s[a]
    return s

def readData(problemFileName = "ProbA.txt"):
    path = "/home/ben/Documents/uni/760/assignment1/data/data/"
    postionsFileName = "Positions.txt"

    containers = np.genfromtxt(path+problemFileName, skip_header=1)
    containers = np.pad(containers, (0,n-len(containers)), 'constant', constant_values=(0,0))
    _, x, y = np. genfromtxt(path+postionsFileName, skip_header=1).T

    return containers, x, y

if __name__=="__main__":

    containers, x, y = readData()
    n = 120
    np.random.shuffle(containers)

    