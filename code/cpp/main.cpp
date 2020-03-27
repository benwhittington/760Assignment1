#include <iostream>
#include <string>
#include <fstream>

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

void readData(std::string fname, int n);

template<typename T>
void swap(T* s, int a, int b) {
    T temp = s[a];
    s[a] = s[b];
    s[b] = temp;
}

void shuffle(double* s, size_t size) {
    int j;
    for(unsigned int i = 0; i < size; ++i) {
        j = i + 1 + (int)(rand()/(1.0 + RAND_MAX) * (120 - i));
        swap(s, i, j);
    }
}


inline float objective(float dx, float dy);

int nextDescentIteration(double* s, double* weights, double* x, double* y, double massTotal, double obj, int startI, int startJ);

int runNextDescent(double* s, double* weights, double* x, double* y, double massTotal, double& obj);

void readData(std::string fname, double* x, double* y) {
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

void fillSolution(double* s, size_t n) {
    unsigned int i = 1;
    for(;i <= n; ++i) {
        s[i] = i;
    }
    for(;i < 120; ++i) {
        s[i] = 0;
    }
}

int main(void) {
    double* x = new double[120];
    double* y = new double[120];
    double* weights;
    double* s = new double[120];

    std::string fname = "Positions.txt";
    std::string problem = "ProbA.txt";
    readData(fname, x, y);
    const unsigned int n = getProblemSize(problem);
    weights = new double[n];
    readData(problem, weights);

    fillSolution(s, n);
    shuffle(s, n);

    for(int i = 0;i < 120; ++i) {
        std::cout << s[i] << std::endl;
    }

    delete x;
    delete y;
    delete weights;

    return 0;
}