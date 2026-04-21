#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <limits>

// Define types with std:: namespace
typedef std::vector<double> VectorDouble;
typedef std::vector<VectorDouble> MatrixDouble;
typedef std::vector<int> VectorInt;
typedef std::vector<VectorInt> MatrixInt;

// Easy functions
int test();
double calc_LA(VectorDouble x, VectorDouble y, VectorDouble z);

//// LSA section
class LSA_Data {
public:
    int max_shift;
    VectorDouble X;
    VectorDouble Y;
    LSA_Data() : max_shift(std::numeric_limits<int>::infinity()) { }
    LSA_Data(int shift, VectorDouble x, VectorDouble y) 
        : max_shift(shift), X(x), Y(y) { }
    void assign(int shift, VectorDouble x, VectorDouble y) {
        max_shift = shift;
        X = x;
        Y = y;
    }
};

class LSA_Result {
public:
    double score;
    MatrixInt trace;
    LSA_Result() : score(0.0), trace() { }
};

LSA_Result DP_lsa(const LSA_Data& data, bool keep_trace = true);

//// LLA section
class LLA_Data {
public:
    int max_shift;
    VectorDouble X;
    VectorDouble Y;
    VectorDouble Z;
    LLA_Data() : max_shift(std::numeric_limits<int>::infinity()) { }
    LLA_Data(int shift, VectorDouble x, VectorDouble y, VectorDouble z)
        : max_shift(shift), X(x), Y(y), Z(z) { }
    void assign(int shift, VectorDouble x, VectorDouble y, VectorDouble z) {
        max_shift = shift;
        X = x;
        Y = y;
        Z = z;
    }
};

class LLA_Result {
public:
    double score;
    MatrixInt trace;
    LLA_Result() : score(0.0), trace() { }
};

LLA_Result DP_lla(const LLA_Data& data, bool keep_trace = true);
