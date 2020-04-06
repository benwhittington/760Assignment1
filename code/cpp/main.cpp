#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <time.h>

void makeTokens(std::string s, double* out);

template<typename T>
void swap(T* s, int a, int b); // done

void readData(std::string fname, float* x, float* y); // done

void readData(std::string fname, float* weights); // done

int getProblemSize(std::string fname); // done

template<typename T>
void fillSolution(T* s, size_t n); // done

template<typename S, typename F>
F objective(S* s, const F* weights, const F* x, const F* y);

template<typename S, typename F>
int nextDescentIteration(S* s, const F* weights, const F* x, const F* y, const F massTotal, F& obj, int& startI, int& startJ);

int runNextDescent(int* s, double* weights, double* x, double* y, double massTotal, double& obj);

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
F objective(S* s, const F*  weights, const F* x, const F* y) {
    F dx = 0;
    F dy = 0;

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
    dx += - (x[a] * weights[s[a]] + x[b] * weights[s[b]]) + (x[a] * weights[s[b]] + x[b] * weights[s[a]]);
    dy += - (y[a] * weights[s[a]] + y[b] * weights[s[b]]) + (y[a] * weights[s[b]] + y[b] * weights[s[a]]);
}

template<typename S, typename F>
int nextDescentIteration(S* s, F* weights, const F* x, const F* y, const F massTotal, F& obj, int& startI, int& startJ, F& dx, F& dy, F* objStore, int& objIdx) {
    int i;
    int j;
    F nextObj;
    F testObj;
    F tempDx, tempDy;
    F testDx, testDy;

    std::cout << "Inside obj: " << obj << std::endl;

    for(int itr = 0;itr < 120; ++itr) {
        for(int jtr = 0; jtr < itr; ++jtr) {
            i = (itr + startI) % 120;
            j = (jtr + startJ) % 120;
            // if(j + 60 == i) { std::cout << j << " + 60 = " << i << std::endl; }
            if((s[i] == 0 && s[j] == 0) || (j + 60 == i)) { continue; }

            swap(s, i, j);
            computeDxDy(s, weights, x, y, testDx, testDy);
            testObj = objective(s, weights, x, y);
            swap(s, i, j); // swap i and j back

            tempDx = dx;
            tempDy = dy;

            centerOfMass(s, i, j, weights, x, y, tempDx, tempDy);
            // std::cout << tempDx << std::endl;
            nextObj = fabs(tempDx) + 5 * fabs(tempDy);
            // std::cout << testObj << (testObj == nextObj ? " == " : " != ") << nextObj << std::endl; 
            // std::cout << nextObj << ", " << obj << std::endl;
            std::cout << tempDx << ", " << testDx << std::endl;
            std::cout << tempDy << ", " << testDy << std::endl;
            std::cout << "-----------------------" << std::endl;

            if(objIdx != -1) { 
                objStore[objIdx] = nextObj; 
                ++objIdx;
            }

            if(nextObj < obj) { 
                obj = nextObj;
                startI = i;
                startJ = j;
                dx = tempDx;
                dy = tempDy;
                swap(s, i, j);
                std::cout << "done" << std::endl;
                return 0;
            }
            // swap(s, i, j); // swap i and j back
        }
    }
    return 1;
}

int runNextDescent(int* s, double* weights, double* x, double* y, double massTotal, double& obj) {
    return 0;
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
    std::string path  = "/home/benwhittington/Documents/uni/760/760Assignment1/data/";
    std::string fname = "Positions.txt";
    std::string problem = "ProbA.txt";
    const int n = getProblemSize(path, problem);

    // alloc arrays. Could be done on stack but ¯\_(ツ)_/¯
    int* s = new int[120];
    int* bestS = new int[120];
    double* x = new double[120];
    double* y = new double[120];
    double* weights = new double[n+1];

    double* bestObjStore;
    double* objStore;
    int* bestObjIndex;

    const bool store = false;

    if(store) {
        bestObjStore = new double[400];
        objStore = new double[100000];
        bestObjIndex = new int[400];
    }


    // read in problem and coords of container positions
    readData(path, fname, x, y);
    readData(path, problem, weights);

    double massTotal = sum(&weights[0], &weights[n+1]);

    // fills positions with containers in problem and empty containers for the rest
    fillSolution(s, n);
    srand(1000000);

    // initial COM location for fast obj calc
    double dx = 0;
    double dy = 0;

    std::random_shuffle(&s[0], &s[120]);

    computeDxDy(s, weights, x, y, dx, dy);
    
    std::cout << "dx: " << dx << " | ";
    std::cout << "dy: " << dy << std::endl;

    int a = 5, b = 5;

    centerOfMass(s, a, b, weights, x, y, dx, dy);
    std::cout << "dx: " << dx << " | ";
    std::cout << "dy: " << dy << std::endl;

    swap(s, a, b);
    computeDxDy(s, weights, x, y, dx, dy);
    std::cout << "dx: " << dx << " | ";
    std::cout << "dy: " << dy << std::endl;

    double obj = objective(s, weights, x, y);
    std::cout << obj << std::endl;
    swap(s, a, b);
    obj = objective(s, weights, x, y);
    std::cout << obj << std::endl;    
    swap(s, a, b);
    obj = objective(s, weights, x, y);
    std::cout << obj << std::endl;
    /**/

    /*

    double obj;
    double bestObj = 10000000;

    int startI, startJ;             // for continuing through neighborhood pick up where left off
    int res;                        // status of descent

    int noDescents = 0;
    int bestIdx = 0;                // index for storing best objective (plotting)
    int objIdx = store ? 0 : -1;    // index for storing all computed objs (-1 if not plotting)
    int indexStore[3];              // store indices of restarts (plotting)


    const int noRestarts = 1;
    int restartNo = 0;

    const double runTime = 40;      // mins
    clock_t start, end;
    double elapsed = 0;
    start = clock();
    
    double testDx, testDy;

    for(; restartNo < noRestarts; ) {
    // while(elapsed < runTime * 60) {
        std::random_shuffle(&s[0], &s[120]);
        obj = objective(s, weights, x, y);
        
        res = 0;
        startJ = 0;
        startI = 0;
        
        while(res != 1) {

            // std::cout << obj << std::endl;
            std::cout << "Outside obj: " << obj << std::endl;
            res = nextDescentIteration(s, weights, x, y, massTotal, obj, startI, startJ, dx, dy, objStore, objIdx);
            computeDxDy(s, weights, x, y, testDx, testDy);
            // std::cout << testDx << ", " << testDy << std::endl;

            ++noDescents;

            if(obj < bestObj) { 
                copy(&bestS[0], &bestS[120], &s[0], &s[120]);
                bestObj = obj;
                if(store) {
                    bestObjIndex[bestIdx] = objIdx;
                    bestObjStore[bestIdx] = bestObj / massTotal;
                    ++bestIdx;
                }
            }
        }

        if(store) {
            bestObjStore[bestIdx] = bestObj;
            bestObjIndex[bestIdx] = objIdx;
            ++bestIdx;
            indexStore[restartNo] = bestIdx;
        }

        if(restartNo % 1000 == 0) {
            end = clock();
            elapsed = (double)(end - start)/(double)CLOCKS_PER_SEC;

            std::cout 
                << "best objective: " << bestObj
                << "\nnumber of descents: " << noDescents
                << "\nrestart: " << restartNo
                << "\ntime elapsed: " << elapsed / 60 << "mins" 
                << "\n------------------------------------------"
                << std::endl;
        }
        ++restartNo;
    }

    writeSolution("sln", bestS, bestS + 120);

    std::cout << "Best objective: " << bestObj / massTotal << std::endl;

    if(store) {
        writeSolution("allObj", &objStore[0], &objStore[objIdx]);
        writeSolution("bestObj", &bestObjStore[0], &bestObjStore[bestIdx]);
        writeSolution("bestObjIdx", &bestObjIndex[0], &bestObjIndex[bestIdx]);
        writeSolution("indices", &indexStore[0], &indexStore[3]);
    }

    /**/

    delete x;
    delete y;
    delete s;
    delete weights;
    delete bestS;

    if(store) {
        delete bestObjIndex;      
        delete objStore;        
        delete bestObjStore;            
    }

}

int main(void) {
    run();
}