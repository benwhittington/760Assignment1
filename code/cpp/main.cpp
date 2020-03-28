#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <math.h>

void makeTokens(std::string s, double* out);

template<typename T>
void swap(T* s, int a, int b); // done

// template<typename T>
// void shuffle(T* s, size_t size); // done

void readData(std::string fname, float* x, float* y); // done

void readData(std::string fname, float* weights); // done

int getProblemSize(std::string fname); // done

template<typename T>
void fillSolution(T* s, size_t n); // done

template<typename S, typename F>
F objective(S* s, const F* weights, const F* x, const F* y, const F massTotal, F obj);

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
void readData(std::string fname, T* x, T* y) {
    std::string path = "/home/ben/Documents/uni/760/assignment1/data/";

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

void readData(std::string fname, double* weights) {
    std::string path = "/home/ben/Documents/uni/760/assignment1/data/";
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

int getProblemSize(std::string fname) {
    std::string path = "/home/ben/Documents/uni/760/assignment1/data/";
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
F objective(S* s, const F*  weights, const F* x, const F* y, const F massTotal, F obj) {
    F dx = 0;
    F dy = 0;

    for(int i = 0;i < 120; ++i) {
        dx += x[i] * weights[s[i]];
        dy += y[i] * weights[s[i]];
    }
    return fabs(dx) + 5 * fabs(dy);
}

template<typename S, typename F>
int nextDescentIteration(S* s, F* weights, const F* x, const F* y, const F massTotal, F& obj, int& startI, int& startJ) {
    int i;
    int j;
    F nextObj;

    for(int itr = 0;itr < 120; ++itr) {
        for(int jtr = 0; jtr < itr; ++jtr) {
            i = (itr + startI) % 120;
            j = (jtr + startJ) % 120;

            if((s[i] == 0 && s[j] == 0) || (j + 60 == i)) { continue; }

            swap(s, i, j);
            nextObj = objective(s, weights, x, y, massTotal, obj);

            if(nextObj < obj) { 
                obj = nextObj;
                startI = i;
                startJ = j;
                return 0;
            }
            swap(s, i, j); // swap i and j back
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
        std::string fname = "Positions.txt";
    std::string problem = "ProbA.txt";
    const int n = getProblemSize(problem);

    int* s = new int[120];
    int* bestS = new int[120];
    double* x = new double[120];
    double* y = new double[120];
    double* weights = new double[n+1];

    // double* bestObjStore = new double[3000];
    // std::fill(&bestObjStore[0], &bestObjStore[3000], -1);

    readData(fname, x, y);
    readData(problem, weights);

    double massTotal = sum(&weights[0], &weights[n+1]);

    fillSolution(s, n);

    srand(1000000);
    std::random_shuffle(&s[0], &s[120]);

    double obj;
    double bestObj = 1000;
    int startI;
    int startJ;
    int res;
    int i = 0;
    const int noRestarts = 40000;

    for(int j = 0; j < noRestarts; ++j) {
        std::random_shuffle(&s[0], &s[120]);
        obj = objective(s, weights, x, y, massTotal, obj);
        res = 0;
        startJ = 0;
        startI = 0;

        while(res != 1) {
            res = nextDescentIteration(s, weights, x, y, massTotal, obj, startI, startJ);
            ++i;
        }

        if(obj < bestObj) { 
            copy(&bestS[0], &bestS[120], &s[0], &s[120]); 
            bestObj = obj;
        }
        if(j % 1000 == 0) { std::cout << "best objective: " << bestObj/massTotal << " iteration: " << i << " restart: " << j << std::endl; }
    }

    std::cout << "Best objective: " << bestObj/massTotal << std::endl;
    // print(&s[0], &s[120]);
    writeSolution("testSln", &bestS[0], &bestS[120]);
    // writeSolution("objTest", &bestObjStore[0], &bestObjStore[i]);

    delete x;
    delete y;
    delete s;
    delete weights;
    delete bestS;
    // delete bestObjStore;
}

int main(void) {
    // std::cout << "done" << std::endl;
    // run();

    // double* a = new double[10];
    // double* b = new double[10];

    // std::fill(&a[0], &a[10], -1);

    // copy(&b[0], &b[10], &a[0], &a[10]);

    // print(&b[0], &b[10]);
    

}