#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

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
float objective(S* s, F*  weights, F* x, F* y, F*, F massTotal, F obj, F startI, F startJ, size_t n);

int nextDescentIteration(int* s, double* weights, double* x, double* y, double massTotal, double obj, int startI, int startJ);

int runNextDescent(int* s, double* weights, double* x, double* y, double massTotal, double& obj);

void makeTokens(std::string s, double* out) {
    std::string delimiter = "  ";

    size_t pos = 0;
    std::string token;
    int i = 0;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        out[i] = std::stof(token);
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

// template<typename T>
// void shuffle(T* s, size_t size) {
//     unsigned int j;
//     for(unsigned int i = 0; i < 120; ++i) {
//         j = i + (int)(rand()/(1.0 + RAND_MAX) * (119 - i));
//         // std::cout << i << ", " << j << std::endl;
//         std::cout << i << ", " << j << std::endl;
//         swap(s, i, j);
//     }
// }

void readData(std::string fname, float* x, float* y) {
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

void readData(std::string fname, float* weights) {
    std::string path = "/home/ben/Documents/uni/760/assignment1/data/";
    std::ifstream weightsStream;
    weightsStream.open(path+fname);


    if(weightsStream.is_open()) {
        std::string line;
        getline(weightsStream, line);
        const int n = std::stoi(line);
        
        for(int i = 0; i < n; ++i) {
            getline(weightsStream, line);
            weights[i] = std::stof(line);
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
    for(;i <= 120; ++i) {
        s[i-1] = 0;
    }
}

template<typename S, typename F>
float objective(S* s, F*  weights, F* x, F* y, F*, F massTotal, F obj, F startI, F startJ, size_t n) {
    F dx = 0;
    F dy = 0;

    for(unsigned int i = 0;i < n; ++i) {
        dx += x[i] * weights[s[i]];
        dy += y[i] * weights[s[i]];
    }
    return abs(dx) + 5 * abs(dy);
}

int nextDescentIteration(int* s, double* weights, double* x, double* y, double massTotal, double obj, int startI, int startJ) {

}

int runNextDescent(int* s, double* weights, double* x, double* y, double massTotal, double& obj) {

}

template<typename T>
void print(T* s) {
    for(unsigned int i = 0;i < 120; ++i) {
        std::cout << s[i] << std::endl;
    }
}

int main(void) {
    std::string fname = "Positions.txt";
    std::string problem = "ProbA.txt";
    const unsigned int n = getProblemSize(problem);

    int* s = new int[120];
    float* x = new float[120];
    float* y = new float[120];
    float* weights = new float[n];

    readData(fname, x, y);
    readData(problem, weights);

    fillSolution(s, n);

    srand(5);

    std::random_shuffle(&s[0], &s[120]);

    print(s);

    delete x;
    delete y;
    delete s;
    delete weights;

    return 0;
}