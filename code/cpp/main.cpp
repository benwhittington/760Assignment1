#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <time.h>

#define STORE_SIZE 1000000

void makeTokens(std::string s, double* out);

template<typename T>
void swap(T* s, int a, int b); // done

void readData(std::string fname, double* x, double* y); // done

void readData(std::string fname, double* weights); // done

int getProblemSize(std::string fname); // done

template<typename T>
void fillSolution(T* s, size_t n); // done

template<typename S, typename F>
F objective(const S* s, const F* weights, const F* x, const F* y);

template<typename S, typename F>
int nextDescentIteration(S* s, const F* weights, const F* x, const F* y, const F massTotal, F& obj, int& startI, int& startJ);

int runNextDescent(int* s, double* weights, double* x, double* y, double massTotal, double& obj);

template<typename T>
void copy(T* aBegin, T* aEnd, const T* bBegin, const T* bEnd);

template<typename T>
void print(const T* start, const T* end);

template<typename T>
void writeSolution(std::string fname, const T* begin, const T* end);

void makeTokens(std::string s, double* out) { // cannot template easily because of stof vs. stod.
    std::string delimiter = "  ";

    size_t pos = 0;
    std::string token;
    int i = 0;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        out[i] = std::stod(token);
        s.erase(0, pos + delimiter.length());
        ++i;
    }
    out[i] = std::stof(s);
}

template<typename T>
void swap(T* s, int a, int b) {
    T temp = s[a];
    s[a] = s[b];
    s[b] = temp;
}

template<typename T>
void readData(std::string path, std::string fname, T* x, T* y) {
    // std::string path = "/home/ben/Documents/uni/760/assignment1/data/";

    std::ifstream positionsStream;
    positionsStream.open(path + "Positions.txt"); //put your program together with thsi file in the same folder.
    double* temp = new double[3];

    if(positionsStream.is_open()){
        std::string line;
        getline(positionsStream, line);

        const int n = std::stoi(line);

        for(int i = 0;i < n; ++i) {
            getline(positionsStream, line); //read number
            makeTokens(line, temp);
            x[i] = temp[1];
            y[i] = temp[2];
        }        
    }
    delete temp;
}

void readData(std::string path, std::string fname, double* weights) {
    // std::string path = "/home/ben/Documents/uni/760/assignment1/data/";
    std::ifstream weightsStream;
    weightsStream.open(path+fname);
    weights[0] = 0;

    if(weightsStream.is_open()) {
        std::string line;
        getline(weightsStream, line);
        const int n = std::stoi(line);
        
        for(int i = 1; i <= n; ++i) {
            getline(weightsStream, line);
            weights[i] = std::stod(line);
        }
    }
}

int getProblemSize(std::string path, std::string fname) {
    // std::string path = "/home/ben/Documents/uni/760/assignment1/data/";
    std::ifstream weightsStream;
    weightsStream.open(path+fname);


    if(weightsStream.is_open()) {
        std::string line;
        getline(weightsStream, line);
        return std::stoi(line);
    }

    return -1;
}

template<typename T>
void fillSolution(T* s, size_t n) {
    unsigned int i = 1;
    for(;i <= n; ++i) {
        s[i-1] = i;
    }
    for(;i <= (unsigned)120; ++i) {
        s[i-1] = 0;
    }
}

template<typename S, typename F>
F objective(const S* s, const F*  weights, const F* x, const F* y) {
    F dx = 0.0;
    F dy = 0.0;

    for(int i = 0;i < 120; ++i) {
        dx += x[i] * weights[s[i]];
        dy += y[i] * weights[s[i]];
    }
    return fabs(dx) + 5 * fabs(dy);
}

template<typename S, typename F>
void computeDxDy(S* s, const F* weights, const F* x, const F* y, F& dx, F& dy) {
    dx = 0;
    dy = 0;

    for(int i = 0; i < 120; ++i) {
        dx += x[i] * weights[s[i]];
        dy += y[i] * weights[s[i]];
    }
}

template<typename S, typename F>
void centerOfMass(const S* s, const int a, const int b, const F* weights, const F* x, const F* y, F& dx, F& dy) {
    dx = dx - (x[a] * weights[s[a]] + x[b] * weights[s[b]]) + (x[a] * weights[s[b]] + x[b] * weights[s[a]]);
    dy = dy - (y[a] * weights[s[a]] + y[b] * weights[s[b]]) + (y[a] * weights[s[b]] + y[b] * weights[s[a]]);
}

template<typename S, typename F>
void kicks(S* s, const F* weights, const F* x, const F* y, const int n, F* bestObjStore, int* bestObjIdx, int& objIdx, int& bestIdx, const F massTotal) {
    int c = 0;
    int i, j;
    while(c < n) {
        i = (int)(119 * (double(rand()) / RAND_MAX));
        j = (int)(119 * (double(rand()) / RAND_MAX));
        if(i == j || (s[i] == 0 && s[j] == 0)) { continue; }
        
        swap(s, i, j);
        ++c;

        if(bestIdx != -1) { 
            bestObjStore[bestIdx] = objective(s, weights, x, y) / massTotal; 
            bestObjIdx[bestIdx] = objIdx;
            ++bestIdx;
            ++objIdx;
        }
    }
}

int countValidNeighbors(int* begin, int* end, const int n, int itrCount) {
    int* p = begin;
    int res = 0;
    while(p < end) {
        res += (itrCount - *p < n) ? 0 : 1;
        // std::cout << "tabulist elem: " << *p << " itrCount: " << itrCount << " diff: " << itrCount - *p << std::endl;
        ++p;
    }
    // exit(0);
    return res;
}

template<typename S, typename F>
void tabuIteration(S* s, S* tabuList, const int tabuListLength, const F* weights, const F* x, const F* y, const F massTotal, F& obj, F& dx, F& dy, F* objStore, int& objIdx, const int itrCount) {
    F nextObj, tempDx, tempDy;
    F bestObj = INFINITY;
    F bestDx, bestDy;
    int bestI, bestJ;
    int i, j;

    for(i = 0;i < 120; ++i) {
        for(j = 0; j < i; ++j) {
            if((s[i] == 0 && s[j] == 0) || (j + 60 == i)) { continue; }

            if(itrCount - tabuList[s[i]] < tabuListLength) {
                ++i;
                if(i > 119) { goto endTabuLoop; }
                else { j = 0; }
            }
            if(itrCount - tabuList[s[j]] < tabuListLength) {
                continue;
            }

            tempDx = dx;
            tempDy = dy;

            centerOfMass(s, i, j, weights, x, y, tempDx, tempDy);
            nextObj = fabs(tempDx) + 5 * fabs(tempDy);

            // for plotting
            if(objIdx != -1) { 
                objStore[objIdx] = nextObj / massTotal;
                ++objIdx;
                if(objIdx > STORE_SIZE) { throw std::out_of_range("Store index exceeded " + std::to_string(STORE_SIZE)); }
            }

            if(nextObj < bestObj) { 
                bestObj = nextObj;
                bestI = i;
                bestJ = j;
                bestDx = tempDx;
                bestDy = tempDy;
            }
        }
    }

    endTabuLoop:
    dx = bestDx;
    dy = bestDy;
    obj = bestObj;

    swap(s, bestI, bestJ);
    tabuList[s[bestI]] = itrCount;
    tabuList[s[bestJ]] = itrCount;
}

template<typename S, typename F>
int nextDescentIteration(S* s, const F* weights, const F* x, const F* y, const F massTotal, F& obj, int& startI, int& startJ, F& dx, F& dy, F* objStore, int& objIdx) {
    int i, j;
    F nextObj, tempDx, tempDy;

    for(int itr = 0;itr < 120; ++itr) {
        for(int jtr = 0; jtr < itr; ++jtr) {
            i = (itr + startI) % 120;
            j = (jtr + startJ) % 120;

            if((s[i] == 0 && s[j] == 0) || (j + 60 == i)) { continue; }

            tempDx = dx;
            tempDy = dy;

            centerOfMass(s, i, j, weights, x, y, tempDx, tempDy);
            nextObj = fabs(tempDx) + 5 * fabs(tempDy);

            // for plotting
            if(objIdx != -1) { 
                objStore[objIdx] = nextObj / massTotal; 
                ++objIdx;
            }

            if(nextObj < obj) { 
                obj = nextObj;
                startI = i;
                startJ = j;
                dx = tempDx;
                dy = tempDy;
                swap(s, i, j);
                return 0;

            }
        }
    }
    return 1;
}

template<typename S, typename F>
void runTabu(S* s, const int tabuListLength, const F* weights, const F* x, const F* y, const F massTotal, const bool store, const bool randRestarts, F& obj) {
    double bestObj;
    double currentObj;

    double dx, dy;

    int noDescents = 0;
    int bestIdx = store ? 0: -1;    // index for storing best objective (plotting)
    int indexStore[3];              // store indices of restarts (plotting)

    const int noRestarts = 1;
    int restartNo = 0;
    int itrCount;

    int objIdx = 0;
    double nextObj;
    double tabuObj, bestTabuObj;
    double tempDx, tempDy;
    double bestDx, bestDy;
    int bestI, bestJ;
    int i, j;

    int* tabuList = new int[120];
    int* bestS = new int[120];

    std::fill(tabuList, tabuList + 120, -tabuListLength); // initially un-ban all containers

    double* bestObjStore; // for storing data to write to file
    double* objStore;
    int* bestObjIndex;

    if(store) {
        bestObjStore = new double[1000];
        objStore = new double[STORE_SIZE];
        bestObjIndex = new int[1000];
    }

    for(; restartNo < noRestarts; ) {
        computeDxDy(s, weights, x, y, dx, dy);

        nextObj = objective(s, weights, x, y);

        if(store) {
            bestObjIndex[bestIdx] = objIdx;
            bestObjStore[bestIdx] = nextObj / massTotal;
            ++bestIdx;
        }

        currentObj = bestObj = nextObj;
        itrCount = 0;

        // while(res != 1) {
        for(int k = 0; k < 100; ++k) {

            // res = nextDescentIteration(s, weights, x, y, massTotal, obj, startI, startJ, dx, dy, objStore, objIdx);
            // tabuIteration(s, tabuList, tabuLength, weights, x, y, massTotal, nextObj, dx, dy, objStore, objIdx, itrCount);
            bestTabuObj = INFINITY;

            for(i = 0;i < 120; ++i) {
                for(j = 0; j < i; ++j) {
                    if((s[i] == 0 && s[j] == 0) || (j + 60 == i)) { continue; }

                    if(itrCount - tabuList[s[i]] < tabuListLength) {
                        ++i;
                        if(i > 119) { goto endTabuLoop; }
                        else { j = 0; }
                    }
                    if(itrCount - tabuList[s[j]] < tabuListLength) {
                        continue;
                    }

                    tempDx = dx;
                    tempDy = dy;

                    centerOfMass(s, i, j, weights, x, y, tempDx, tempDy);
                    tabuObj = fabs(tempDx) + 5 * fabs(tempDy);

                    if(store) { 
                        objStore[objIdx] = tabuObj;// / massTotal;
                        ++objIdx;
                        if(objIdx > STORE_SIZE) { throw std::out_of_range("Store index exceeded " + std::to_string(STORE_SIZE)); }
                    }

                    if(tabuObj < bestTabuObj) { 
                        bestTabuObj = tabuObj;
                        bestI = i;
                        bestJ = j;
                        bestDx = tempDx;
                        bestDy = tempDy;
                    }
                }
            }

            endTabuLoop:
            dx = bestDx;
            dy = bestDy;
            nextObj = bestTabuObj;

            swap(s, bestI, bestJ);
            tabuList[s[bestI]] = itrCount;
            tabuList[s[bestJ]] = itrCount;

            if(store) {
                bestObjIndex[bestIdx] = objIdx;
                bestObjStore[bestIdx] = currentObj / massTotal;
                ++bestIdx;
                bestObjIndex[bestIdx] = objIdx;
                bestObjStore[bestIdx] = nextObj / massTotal;                
                ++bestIdx;
            }

            currentObj = nextObj;

            if(currentObj < bestObj) { // keep s if obj is best so far
                copy(&bestS[0], &bestS[120], &s[0], &s[120]);
                bestObj = currentObj;
            }
            ++noDescents;
            ++itrCount;
        }

        if(store) {
            indexStore[restartNo] = bestIdx;
        }

        if(randRestarts) {
            std::random_shuffle(&s[0], &s[120]);
        } else {
            kicks(s, weights, x, y, 15, bestObjStore, bestObjIndex, objIdx, bestIdx, massTotal);
        }
        ++restartNo;
    }

    writeSolution("sln", bestS, bestS + 120);

    std::cout << "Best objective: " << currentObj / massTotal << std::endl;

    if(store) {
        writeSolution("allObj", &objStore[0], &objStore[objIdx]);
        writeSolution("bestObj", &bestObjStore[0], &bestObjStore[bestIdx]);
        writeSolution("bestObjIdx", &bestObjIndex[0], &bestObjIndex[bestIdx]);
        writeSolution("indices", &indexStore[0], &indexStore[3]);
    } 

    delete tabuList;
    delete bestS;

    if(store) {
        delete objStore;    
        delete bestObjStore;            
        delete bestObjIndex;      
    }
}

template<typename S, typename F>
void runNextDescent(S* s, const int tabuListLength, const F* weights, const F* x, const F* y, const F massTotal, const bool store, const bool randRestarts, F& obj) {
    double bestObj;
    double currentObj;

    double dx, dy;

    int noDescents = 0;
    int bestIdx = store ? 0: -1;    // index for storing best objective (plotting)
    int indexStore[3];              // store indices of restarts (plotting)

    const int noRestarts = 1;
    int restartNo = 0;
    int itrCount;
    int res = 0;

    int objIdx = 0;
    double nextObj;
    double tempDx, tempDy;
    int i, j;
    int startI = 0, startJ = 0;

    int* bestS = new int[120];

    double* bestObjStore; // for storing data to write to file
    double* objStore;
    int* bestObjIndex;

    if(store) {
        bestObjStore = new double[1000];
        objStore = new double[STORE_SIZE];
        bestObjIndex = new int[1000];
    }

    for(; restartNo < noRestarts; ) {
        computeDxDy(s, weights, x, y, dx, dy);

        nextObj = objective(s, weights, x, y);

        if(store) {
            bestObjIndex[bestIdx] = objIdx;
            bestObjStore[bestIdx] = nextObj / massTotal;
            ++bestIdx;
        }

        currentObj = bestObj = nextObj;
        itrCount = 0;

        while(res != 1) {
            double descentObj; // replace nextObj
            // res = nextDescentIteration(s, weights, x, y, massTotal, obj, startI, startJ, dx, dy, objStore, objIdx);
            // tabuIteration(s, tabuList, tabuLength, weights, x, y, massTotal, nextObj, dx, dy, objStore, objIdx, itrCount);
            
            for(int itr = 0;itr < 120; ++itr) {
                for(int jtr = 0; jtr < itr; ++jtr) {
                    i = (itr + startI) % 120;
                    j = (jtr + startJ) % 120;

                    if((s[i] == 0 && s[j] == 0) || (j + 60 == i)) { continue; }

                    tempDx = dx;
                    tempDy = dy;

                    centerOfMass(s, i, j, weights, x, y, tempDx, tempDy);
                    descentObj = fabs(tempDx) + 5 * fabs(tempDy);

                    // for plotting
                    if(objIdx != -1) { 
                        objStore[objIdx] = descentObj / massTotal; 
                        ++objIdx;
                    }

                    if(nextObj < obj) { 
                        nextObj = descentObj;
                        startI = i;
                        startJ = j;
                        dx = tempDx;
                        dy = tempDy;
                        swap(s, i, j);
                        res = 0;
                        goto endDescentLoop;

                    }
                }
            }
            res = 1;
            endDescentLoop:

            if(store) {
                bestObjIndex[bestIdx] = objIdx;
                bestObjStore[bestIdx] = currentObj / massTotal;
                ++bestIdx;
                bestObjIndex[bestIdx] = objIdx;
                bestObjStore[bestIdx] = nextObj / massTotal;                
                ++bestIdx;
            }

            currentObj = nextObj;

            if(currentObj < bestObj) { // keep s if obj is best so far
                copy(&bestS[0], &bestS[120], &s[0], &s[120]);
                bestObj = currentObj;
            }
            ++noDescents;
            ++itrCount;
        }

        if(store) {
            indexStore[restartNo] = bestIdx;
        }

        if(randRestarts) {
            std::random_shuffle(&s[0], &s[120]);
        } else {
            kicks(s, weights, x, y, 15, bestObjStore, bestObjIndex, objIdx, bestIdx, massTotal);
        }
        ++restartNo;
    }

    writeSolution("sln", bestS, bestS + 120);

    std::cout << "Best objective: " << currentObj / massTotal << std::endl;

    if(store) {
        writeSolution("allObj", &objStore[0], &objStore[objIdx]);
        writeSolution("bestObj", &bestObjStore[0], &bestObjStore[bestIdx]);
        writeSolution("bestObjIdx", &bestObjIndex[0], &bestObjIndex[bestIdx]);
        writeSolution("indices", &indexStore[0], &indexStore[3]);
    } 

    delete bestS;

    if(store) {
        delete objStore;    
        delete bestObjStore;            
        delete bestObjIndex;      
    }
}

template<typename T>
void print(const T* start, const T* end) {
    const T* p = start;
    std::cout << "[ ";
    while(p < end - 1) {
        std::cout << *p << ", ";
        ++p;
    }
    std::cout << *p << " ]" << std::endl;
}

template<typename T>
T sum(T* begin, T* end) {
    T* p = begin;
    T total = 0;
    while(p < end) {
        total += *p;
        ++p;
    }
    return total;
}

template<typename T>
void writeSolution(std::string fname, const T* begin, const T* end) {
    const std::string path = "/home/ben/Documents/uni/760/assignment1/output/";
    std::ofstream outStream;
    outStream.open(path + fname);
    const T* p = begin;

    while(p < end -1) {
        outStream << *p << "\n";
        ++p;
    }
    outStream << *p << std::endl;
}

template<typename T>
void copy(T* aBegin, T* aEnd, const T* bBegin, const T* bEnd) {
    T* a = aBegin;
    const T* b = bBegin;
    while(a < aEnd) {
        *a = *b;
        ++a;
        ++b;
    }
}

void run() {

    // read data
    // std::string path  = "/home/benwhittington/Documents/uni/760/760Assignment1/data/";
    std::string path = "/home/ben/Documents/uni/760/assignment1/data/";
    std::string fname = "Positions.txt";
    std::string problem = "ProbA.txt";
    const int n = getProblemSize(path, problem);
    const int noContainers = 120;

    int* s = new int[noContainers]; // alloc arrays. Could be done on stack but ¯\_(ツ)_/¯
    int* bestS = new int[noContainers];
    double* x = new double[noContainers];
    double* y = new double[noContainers];
    double* weights = new double[n+1];

    const bool randRestarts = false; // false for kicks, true for rand restarts
    const bool store = true;

    double* bestObjStore; // for storing data to write to file
    double* objStore;
    int* bestObjIndex;

    readData(path, fname, x, y); // read in problem and coords of container positions
    readData(path, problem, weights);

    fillSolution(s, n); // fills positions with containers in problem and empty containers for the rest
    srand(1000000);
    std::random_shuffle(&s[0], &s[noContainers]);

    const int tabuLength = std::min(50, int(n/3));
    double massTotal = sum(&weights[0], &weights[n+1]);


    double nextObj;
    double dx = 0;                  // for storing com for fast obj calc
    double dy = 0;
    int objIdx = store ? 0 : -1;    // index for storing all computed objs (-1 if not plotting)
    int itrCount;

    // runTabu(s, tabuLength, weights, x, y, massTotal, store, randRestarts, nextObj);
    runNextDescent(s, tabuLength, weights, x, y, massTotal, store, randRestarts, nextObj);

    /*
    double bestObj;
    double currentObj;

    int startI, startJ;             // for continuing through neighborhood pick up where left off
    int res;                        // status of descent

    int noDescents = 0;
    int bestIdx = store ? 0: -1;    // index for storing best objective (plotting)
    int indexStore[3];              // store indices of restarts (plotting)

    const int noRestarts = 1;
    int restartNo = 0;

    double testDx, testDy;

    for(; restartNo < noRestarts; ) {
        computeDxDy(s, weights, x, y, dx, dy);

        nextObj = objective(s, weights, x, y);

        if(store) {
            bestObjIndex[bestIdx] = objIdx;
            bestObjStore[bestIdx] = nextObj / massTotal;
            ++bestIdx;
        }

        currentObj = bestObj = nextObj;
        startJ = startJ = 0;
        itrCount = 0;
        res = 0;

        // while(res != 1) {
        for(int k = 0; k < 100; ++k) {

            // res = nextDescentIteration(s, weights, x, y, massTotal, obj, startI, startJ, dx, dy, objStore, objIdx);
            tabuIteration(s, tabuList, tabuLength, weights, x, y, massTotal, nextObj, dx, dy, objStore, objIdx, itrCount);

            if(store) {
                bestObjIndex[bestIdx] = objIdx;
                bestObjStore[bestIdx] = currentObj / massTotal;
                ++bestIdx;
                bestObjIndex[bestIdx] = objIdx;
                bestObjStore[bestIdx] = nextObj / massTotal;                
                ++bestIdx;
            }

            currentObj = nextObj;

            if(currentObj < bestObj) { // keep s if obj is best so far
                copy(&bestS[0], &bestS[120], &s[0], &s[120]);
                bestObj = currentObj;
            }
            ++noDescents;
            ++itrCount;
        }

        if(store) {
            indexStore[restartNo] = bestIdx;
        }

        if(randRestarts) {
            std::random_shuffle(&s[0], &s[120]);
        } else {
            kicks(s, weights, x, y, 15, bestObjStore, bestObjIndex, objIdx, bestIdx, massTotal);
        }

        ++restartNo;
    }

    writeSolution("sln", bestS, bestS + 120);

    std::cout << "Best objective: " << currentObj / massTotal << std::endl;

    if(store) {
        writeSolution("allObj", &objStore[0], &objStore[objIdx]);
        writeSolution("bestObj", &bestObjStore[0], &bestObjStore[bestIdx]);
        writeSolution("bestObjIdx", &bestObjIndex[0], &bestObjIndex[bestIdx]);
        writeSolution("indices", &indexStore[0], &indexStore[3]);
    }

    */

    delete x;
    delete y;
    delete s;
    delete weights;
    delete bestS;

    // if(store) {
    //     delete objStore;    
    //     delete bestObjStore;            
    //     delete bestObjIndex;      
    // }

}

int main(void) {
    run();
}