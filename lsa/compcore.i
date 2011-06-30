%module compcore
%include "std_vector.i"
%{
#include "compcore.h"
%}
namespace std{
    %template(VectorDouble) vector<double>;
    %template(VectorInt) vector<int>;
    %template(MatrixDouble) vector<vector<double> >;
    %template(MatrixInt) vector<vector<int> >;
};
%include "compcore.h"  

/* %constant is like #define */
//%constant LSA_Result (*LSA_test)(const LSA_Data&, bool) = DP_lsa; // call LSA_test
//%constant LSA_Result (*LSA_test)(const LSA_Data&, bool) = LSA_DP; // call LSA_test

//%constant LLA_Result (*LLA_test)( const LLA_Data& ) = DP_lla; // call LLA_test
//%constant LA_Result (*LA_test)(const LLA_Data& ) = ST_la; // call LA_test
//%template(LLA_PT) PT<LLA_Data, LLA_Result>;
//%template(LA_PT) PT<LLA_Data, LA_Result>;
//%template(LSA_PT) PT<LSA_Data, LSA_Result>;
