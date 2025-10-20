#include "compcore.hpp"
using namespace std;

/* Test section */
int test(){
  cout<<"tested";
  return 0;
}

/* LSA section */
/*
void LSA_Data::assign(int shift, VectorDouble x, VectorDouble y){
  max_shift = shift;
  X.assign(x.begin(),x.end());
  Y.assign(y.begin(),y.end());
}*/

LSA_Result DP_lsa( const LSA_Data& data, bool keep_trace ){  //python does not support default value
  LSA_Result lsa_result;
  int max_p[2]={0}; int porn=0;
  double max_s=-std::numeric_limits<double>::infinity(); //initialize to negative infinite
  //MatrixDouble psm; MatrixDouble nsm;
  MatrixDouble psm = vector<vector<double> >(data.X.size()+1, vector<double>(data.Y.size()+1));
  MatrixDouble nsm = vector<vector<double> >(data.X.size()+1, vector<double>(data.Y.size()+1));
  //double psm[data.X.size()+1][data.Y.size()+1]; // positive score matrix
  //cout << "here" <<endl;
  //double nsm[data.X.size()+1][data.Y.size()+1]; // negative score matrix
  //cout << data.X.size() <<endl;
  //cout << data.X.size() <<endl;
  //for( unsigned int i=0; i<=data.X.size(); i++ )
  //  for( unsigned int j=0; j<=data.Y.size(); j++ )
  //    { psm[i][j]=0.; nsm[i][j]=0.; } //initialize score matrix elements to 0
  //int tm[data.X.size()+1][data.Y.size()+1]; // trace matrix
  //for( unsigned int i=0; i<=data.X.size(); i++ )
  //  for( unsigned int j=0; j<=data.Y.size(); j++ )
  //    tm[i][j]=0; //initialize trace matrix elements to 0
                        //0x0??? for the moving direction, no gap <=> only 0x0000 and 0x0111
  //cout << tm[3][4][5] << endl;
  for( unsigned int i=1; i<=data.X.size(); i++)
    for( unsigned int j=std::max(1,(int)i-data.max_shift); (int)j<=std::min((int)data.Y.size(),(int)i+data.max_shift); j++ ){
        //if( abs((int)(i-j))>data.max_shift ) continue; // alignment not allowed
        double s1=data.X[i-1]*data.Y[j-1];
        //double s0=sm[i-1][j-1];
        psm[i][j] = std::max(0., psm[i-1][j-1]+s1);
        nsm[i][j] = std::max(0., nsm[i-1][j-1]-s1);
        if ( psm[i][j] >= max_s ) {
          max_p[0] = i; max_p[1] = j; max_s = psm[i][j]; porn=1;
        }
        if ( nsm[i][j] >= max_s ) {
          max_p[0] = i; max_p[1] = j; max_s = nsm[i][j]; porn=-1;
        }
      }
  //calculation of the score matrix
  //cout << "max_s=" << max_s << ";max_p="<< max_p[0] <<","<<max_p[1]<<endl;
  int length=0; vector<int> step; step.push_back(max_p[0]); step.push_back(max_p[1]);
  if (porn == -1) {
    lsa_result.score=-1*nsm[max_p[0]][max_p[1]]/data.X.size();      //though it is not enforced, assert(len(X)==len(Y))
    //cout<<"porn="<<porn<<":score="<<lsa_result.score<<endl;
    while(nsm[max_p[0]-length][max_p[1]-length]!=0. && keep_trace == true) { 
      length++; lsa_result.trace.push_back(step); step.clear(); 
      step.push_back(max_p[0]-length); step.push_back(max_p[1]-length); }
  }
  else {
    lsa_result.score=psm[max_p[0]][max_p[1]]/data.X.size();         //though it is not enforced, assert(len(X)==len(Y))
    //cout<<"porn="<<porn<<":score="<<lsa_result.score<<endl;
    while(psm[max_p[0]-length][max_p[1]-length]!=0. && keep_trace == true) { 
      length++; lsa_result.trace.push_back(step); step.clear(); 
      step.push_back(max_p[0]-length); step.push_back(max_p[1]-length); }
  }
  //cout << "length="<<lsa_result.trace.size()<< endl;
  return lsa_result;
}

/* LLA section */
//definition of functions
LLA_Result DP_lla(const LLA_Data& data, bool keep_trace) {
    // Initialize result
    LLA_Result lla_result;
    
    // Input validation
    if (data.X.empty() || data.Y.empty() || data.Z.empty()) {
        lla_result.score = 0;
        return lla_result;
    }
    if (data.X.size() != data.Y.size() || data.Y.size() != data.Z.size()) {
        lla_result.score = 0;
        return lla_result;
    }

    // Initialize tracking variables
    const size_t n = data.X.size();
    int max_p[3] = {0};
    double max_s = 0.0; // non-negative best absolute score
    int best_len = 0;    // length of best segment (for tie-breaking: prefer shorter)
    int best_sign = 1;   // 1 for positive, -1 for negative winner

    // Positive/Negative score matrices and lengths
    vector<vector<vector<double>>> psm(n + 1,
        vector<vector<double>>(n + 1,
            vector<double>(n + 1, 0.0)));
    vector<vector<vector<double>>> nsm(n + 1,
        vector<vector<double>>(n + 1,
            vector<double>(n + 1, 0.0)));
    vector<vector<vector<int>>> lpsm(n + 1,
        vector<vector<int>>(n + 1,
            vector<int>(n + 1, 0)));
    vector<vector<vector<int>>> lnsm(n + 1,
        vector<vector<int>>(n + 1,
            vector<int>(n + 1, 0)));

    // Direction matrices for backtracking (only when keep_trace is true)
    vector<vector<vector<bool>>> p_ext, n_ext;
    if (keep_trace) {
        p_ext.resize(n + 1, vector<vector<bool>>(n + 1, vector<bool>(n + 1, false)));
        n_ext.resize(n + 1, vector<vector<bool>>(n + 1, vector<bool>(n + 1, false)));
        lla_result.trace.reserve(n);
    }

    // Main DP loop
    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= n; j++) {
            for (size_t k = 1; k <= n; k++) {
                // Shift constraint
                if (data.max_shift != std::numeric_limits<int>::max()) {
                    if (abs((int)i - (int)j) > data.max_shift || 
                        abs((int)i - (int)k) > data.max_shift || 
                        abs((int)j - (int)k) > data.max_shift) {
                        continue;
                    }
                }

                double s1 = data.X[i - 1] * data.Y[j - 1] * data.Z[k - 1];

                // Positive accumulation
                double prev_p = psm[i - 1][j - 1][k - 1];
                double cand_p = prev_p + s1;
                if (cand_p > 0) {
                    psm[i][j][k] = cand_p;
                    lpsm[i][j][k] = lpsm[i - 1][j - 1][k - 1] + 1;
                    if (keep_trace) p_ext[i][j][k] = (lpsm[i - 1][j - 1][k - 1] > 0);
                } else {
                    psm[i][j][k] = 0.0;
                    lpsm[i][j][k] = 0;
                    if (keep_trace) p_ext[i][j][k] = false;
                }

                // Negative accumulation
                double prev_n = nsm[i - 1][j - 1][k - 1];
                double cand_n = prev_n - s1; // accumulate negative
                if (cand_n > 0) {
                    nsm[i][j][k] = cand_n;
                    lnsm[i][j][k] = lnsm[i - 1][j - 1][k - 1] + 1;
                    if (keep_trace) n_ext[i][j][k] = (lnsm[i - 1][j - 1][k - 1] > 0);
                } else {
                    nsm[i][j][k] = 0.0;
                    lnsm[i][j][k] = 0;
                    if (keep_trace) n_ext[i][j][k] = false;
                }

                // Decide current best at (i,j,k)
                double cur_abs = 0.0; int cur_len = 0; int cur_sign = 1;
                if (psm[i][j][k] >= nsm[i][j][k]) {
                    cur_abs = psm[i][j][k];
                    cur_len = lpsm[i][j][k];
                    cur_sign = 1;
                } else {
                    cur_abs = nsm[i][j][k];
                    cur_len = lnsm[i][j][k];
                    cur_sign = -1;
                }

                if (cur_abs > 0) {
                    bool better = false;
                    if (cur_abs > max_s) {
                        better = true;
                    } else if (cur_abs == max_s) {
                        // tie-break: prefer shorter length
                        if (cur_len > 0 && (best_len == 0 || cur_len < best_len)) {
                            better = true;
                        } else if (cur_len == best_len && best_len > 0) {
                            // then lexicographically smallest end index
                            if ( (int)i < max_p[0] ||
                                 ( (int)i == max_p[0] && ( (int)j < max_p[1] || ( (int)j == max_p[1] && (int)k < max_p[2] ) ) ) ) {
                                better = true;
                            }
                        }
                    }
                    if (better) {
                        max_s = cur_abs;
                        best_len = cur_len;
                        best_sign = cur_sign;
                        max_p[0] = (int)i;
                        max_p[1] = (int)j;
                        max_p[2] = (int)k;
                    }
                }
            }
        }
    }

    if (max_s == 0.0) {
        lla_result.score = 0.0;
        return lla_result;
    }

    // Store the maximum score in `lla_result` (normalize by length of X for consistency)
    double signed_score = (best_sign == 1) ? psm[max_p[0]][max_p[1]][max_p[2]] : -nsm[max_p[0]][max_p[1]][max_p[2]];
    lla_result.score = signed_score / data.X.size();

    // Backtrace if requested
    if (keep_trace && best_len > 0) {
        int i = max_p[0], j = max_p[1], k = max_p[2];
        
        // Choose the correct direction matrix based on best_sign
        vector<vector<vector<bool>>>* ext_matrix;
        if (best_sign == 1) {
            ext_matrix = &p_ext;
        } else {
            ext_matrix = &n_ext;
        }
        
        // Debug: print max position and constraints
        // std::cout << "Debug: max_p = [" << i << ", " << j << ", " << k << "], best_sign = " << best_sign << std::endl;
        // std::cout << "Debug: delayLimit = " << data.max_shift << std::endl;
        
        // Backtrack along the path
        while (i > 0 && j > 0 && k > 0 && (*ext_matrix)[i][j][k]) {
            lla_result.trace.push_back({i, j, k});
            i--; j--; k--;
            
            // Ensure we still satisfy delay constraints after decrement
            if (data.max_shift != std::numeric_limits<int>::max()) {
                if (abs(i - j) > data.max_shift || 
                    abs(i - k) > data.max_shift || 
                    abs(j - k) > data.max_shift) {
                    // std::cout << "Debug: Breaking due to constraint violation at [" << i << ", " << j << ", " << k << "]" << std::endl;
                    break;
                }
            }
        }
        // Add the starting position
        lla_result.trace.push_back({i, j, k});
        
    // Debug: print final trace before reverse (only when DEBUG_COMPCORE is defined)
#ifdef DEBUG_COMPCORE
    std::cout << "Debug: Final trace before reverse: ";
    for (const auto& pos : lla_result.trace) {
        std::cout << "[" << pos[0] << ", " << pos[1] << ", " << pos[2] << "] ";
    }
    std::cout << std::endl;
#endif
        
    // Keep trace in end->start order (consistent with DP_lsa and Python expectations)
    }

    return lla_result;
}

/* Static calc LA functions */
double calc_LA(VectorDouble x, VectorDouble y, VectorDouble z) {
    if (x.size() != y.size() || y.size() != z.size()) {
        throw std::runtime_error("All input vectors must have the same length");
    }
    
    double sum = 0.0;
    int n = x.size();
    
    for (int i = 0; i < n; i++) {
        sum += (x[i] * y[i] * z[i]);
    }
    
    return sum / n;
}

/*
LLA_Result ST_lla( const LLA_Data& lla_data ) {
  LLA_Result lla_result;
  double score=0;
  for( unsigned int i=0; i<lla_data.X.size(); i++ ){
    score=score+lla_data.X[i]*lla_data.Y[i]*lla_data.Z[i];
  }
  lla_result.score=score;
  return lla_result;
}
*/

/*
LSA_Result DP_rep_lsa( const LSA_Rep_Data& data ){
  LSA_Result lsa_result;
  int max_p[2]={0}; int porn=0;
  double max_s=-std::numeric_limits<double>::infinity(); //initialize to negative infinite
  double psm[data.X.size()+1][data.Y.size()+1]; // positive score matrix
  double nsm[data.X.size()+1][data.Y.size()+1]; // negative score matrix
  for( unsigned int i=0; i<=data.X.size(); i++ )
    for( unsigned int j=0; j<=data.Y.size(); j++ )
      { psm[i][j]=0.; nsm[i][j]=0.; } //initialize score matrix elements to 0
  //cout << sm[3][4][5] << endl;
  int tm[data.X.size()+1][data.Y.size()+1]; // trace matrix
  for( unsigned int i=0; i<=data.X.size(); i++ )
    for( unsigned int j=0; j<=data.Y.size(); j++ )
      tm[i][j]=0; //initialize trace matrix elements to 0
                        //0x0??? for the moving direction, no gap <=> only 0x0000 and 0x0111
  //cout << tm[3][4][5] << endl;
  for( unsigned int i=1; i<=data.X.size(); i++)
    for( unsigned int j=std::max(1,(int)i-data.max_shift); (int)j<=std::min((int)data.Y.size(),(int)i+data.max_shift); j++ ){
        //if( abs((int)(i-j))>data.max_shift ) continue; // alignment not allowed
        //double s1=data.X[i-1]*data.Y[j-1];
        double s1=mcorr(X[i-1], Y[j-1]);  //Correlation Function
        //double s0=sm[i-1][j-1];
        psm[i][j] = std::max(0., psm[i-1][j-1]+s1);
        nsm[i][j] = std::max(0., nsm[i-1][j-1]-s1);
        //int t=-std::numeric_limits<int>::infinity();
        //if(s0>0)
        //    t=(s1>0)?0x0011:0x0000; //only move forward when strict larger than 0
        //else if(s0<0)
        //    t=(s1<0)?0x0011:0x0000;
        //else
        //  t=0x0000;
        //sm[i][j]=(t==0x0000)?s1:s0+s1;
        //tm[i][j]=t;
        //cout << "sm["<<i<<","<<j<<","<<k<<"]="<<sm[i][j][k]<<";s1="<< s1 << ";s0=" << s0 << ";t="<<hex<<tm[i][j][k]<<endl;
        if ( psm[i][j] > max_s ) {
          max_p[0] = i; max_p[1] = j; max_s = psm[i][j]; porn=1;
        }
        if ( nsm[i][j] > max_s ) {
          max_p[0] = i; max_p[1] = j; max_s = nsm[i][j]; porn=-1;
        }
      }
  //calculation of the score matrix
  //cout << "max_s=" << max_s << ";max_p="<< max_p[0] <<","<<max_p[1]<<endl;
  int length=0; vector<int> step; step.push_back(max_p[0]); step.push_back(max_p[1]);
  if (porn == -1) {
    lsa_result.score=-1*nsm[max_p[0]][max_p[1]];
    //cout<<"porn="<<porn<<":score="<<lsa_result.score<<endl;
    while(nsm[max_p[0]-length][max_p[1]-length]!=0.) { 
      length++; lsa_result.trace.push_back(step); step.clear(); 
      step.push_back(max_p[0]-length); step.push_back(max_p[1]-length); }
  }
  else {
    lsa_result.score=psm[max_p[0]][max_p[1]];
    //cout<<"porn="<<porn<<":score="<<lsa_result.score<<endl;
    while(psm[max_p[0]-length][max_p[1]-length]!=0.) { 
      length++; lsa_result.trace.push_back(step); step.clear(); 
      step.push_back(max_p[0]-length); step.push_back(max_p[1]-length); }
  }
  //for( int i=max_p[0], j=max_p[1]; ;
  //        i-=(int)((tm[i][j]&0x0001)>0), j-=(int)((tm[i][j]&0x0010)>0) ){
  //    vector<int> step;
  //  step.push_back(i); step.push_back(j);
  //  lsa_result.trace.push_back(step);
  //  if (tm[i][j] == 0x0000) break;
  //}
  //cout << "length="<<lsa_result.trace.size()<< endl;
  return lsa_result;
*/

/*
template <class DataType, class ResultType>
  PT_Return PT( DataType data, const ResultType& result, int pn,  ResultType (*test_func)( const DataType& ) ){
      PT_Return pt_return;
      double score;
      int count=0;
      for( int i=0; i<pn; i++ ){
        data.random_shuffle();
        score = (test_func(data)).score;
        cout << score << endl;
        pt_return.scores.push_back(score);
        if (score < result.score) count++;
      }
      //cout << count <<endl;
      pt_return.pvalue = (double)count/(double)pn;
      return pt_return;
}
*/
