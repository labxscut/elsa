#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <limits>
//#include <numeric>

using namespace std;

//Customized types
typedef vector<double> VectorDouble;
typedef vector<VectorDouble> MatrixDouble;
typedef vector<int> VectorInt;
typedef vector<VectorInt> MatrixInt;

//// LSA section

class LSA_Data {
public:
  int max_shift;
	VectorDouble X;
	VectorDouble Y;
	LSA_Data(){ VectorDouble X; VectorDouble Y; max_shift=std::numeric_limits<int>::infinity(); };
	LSA_Data(int shift, VectorDouble x, VectorDouble y): max_shift(shift),X(x),Y(y){ };
  void assign(int, VectorDouble, VectorDouble);
};

class LSA_Result {
public:
  double score;
  MatrixInt trace;
};

LSA_Result DP_lsa( const LSA_Data&, bool ); 


//// LLA and LA data types
/*
class LLA_Data {
public:
  int max_shift;
	VectorDouble X;
	VectorDouble Y;
	VectorDouble Z;
	LLA_Data(){ VectorDouble X; VectorDouble Y; VectorDouble Z; max_shift=std::numeric_limits<int>::infinity(); };
	LLA_Data(int shift, VectorDouble x, VectorDouble y, VectorDouble z): max_shift(shift),X(x),Y(y),Z(z){ };
	//inline int random_shuffle();
};
*/

/*
int LLA_Data::random_shuffle(){
	std::random_shuffle(X.begin(),X.end());
	//cout<<"X="<<X[0]<<","<<X[1]<<","<<X[2]<<","<<X[3]<<endl;
	std::random_shuffle(Y.begin(),Y.end());
	//cout<<"Y="<<Y[0]<<","<<Y[1]<<","<<Y[2]<<","<<Y[3]<<endl;
	std::random_shuffle(Z.begin(),Z.end());
	//cout<<"Z="<<Z[0]<<","<<Z[1]<<","<<Z[2]<<","<<Z[3]<<endl;
	return 0;
};
*/

/*
class LLA_Result {
public:
	double score;
	MatrixInt trace;
};

class LA_Result {
public:
  double score;
};
*/


//// Permutation test template
/*
class PT_Return {
public:
	VectorDouble scores;
	double pvalue;
};
*/

//Declaration of functions
//LLA_Result DP_lla( const LLA_Data& ); //const: passing the reference not for modifying
//LA_Result ST_la( const LLA_Data& );


/*
template <class DataType, class ResultType>
    PT_Return PT( DataType data, const ResultType& result, int pn, ResultType (*test_func)( const DataType& ) )
{
  PT_Return pt_return;
  double score;
  int count=0;
  for( int i=0; i<pn; i++ ){
    data.random_shuffle_x();
    score = (test_func(data)).score;
    //cout << i << "-th: " << score << endl;
    pt_return.scores.push_back(score);
    if ( abs(score) > abs(result.score) ) count++; // for two tailed p-value
  }
  //cout << count <<endl;
  pt_return.pvalue = (double)count/(double)pn;
  return pt_return;
};
*/

//// LSA with replicates data types
/* currently not used
class LSA_Rep_Data{
public:
  int max_shift;
  MatrixDouble X;
  MatrixDouble Y;
  LSA_Rep_Data(){ MatrixDouble X; MatrixDouble Y; max_shift=std::numeric_limits<int>::infinity(); };
  LSA_Rep_Data(int shift, MatrixDouble x, MatrixDouble y): max_shift(shift), X(x), Y(y){ };
  inline int random_shuffle();
};

int LSA_Rep_Data::random_shuffle(){  //shuffle rows, we want each column is a time sequence, do proper transformation at input
  std::random_shuffle(X.begin(),X.end());
  std::random_shuffle(Y.begin(),Y.end());
  return 0;
};
*/
