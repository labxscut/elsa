#include "compcore.h"

//// LSA section
void LSA_Data::assign(int shift, VectorDouble x, VectorDouble y){
  max_shift = shift;
  X.assign(x.begin(),x.end());
  Y.assign(y.begin(),y.end());
}

LSA_Result DP_lsa( const LSA_Data& data, bool keep_trace ){  //python does not support default value
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
				double s1=data.X[i-1]*data.Y[j-1];
        //double s0=sm[i-1][j-1];
        psm[i][j] = std::max(0., psm[i-1][j-1]+s1);
        nsm[i][j] = std::max(0., nsm[i-1][j-1]-s1);
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


//definitions of functions
/*
LLA_Result DP_lla( const LLA_Data& data ){
	LLA_Result lla_result;
	int max_p[3]={0};
	double max_s=-std::numeric_limits<double>::infinity(); //initialize to negative infinite
	double sm[data.X.size()+1][data.Y.size()+1][data.Z.size()+1]; // score matrix
	for( unsigned int i=0; i<=data.X.size(); i++ )
		for( unsigned int j=0; j<=data.Y.size(); j++ )
			for( unsigned int k=0; k<=data.Z.size(); k++ )
				sm[i][j][k]=0.; //initialize score matrix elements to 0
	//cout << sm[3][4][5] << endl;
	int tm[data.X.size()+1][data.Y.size()+1][data.Z.size()+1]; // trace matrix
	for( unsigned int i=0; i<=data.X.size(); i++ )
		for( unsigned int j=0; j<=data.Y.size(); j++ )
			for( unsigned int k=0; k<=data.Z.size(); k++ )
				tm[i][j][k]=0; //initialize trace matrix elements to 0
												//0x0??? for the moving direction, no gap <=> only 0x0000 and 0x0111
	//cout << tm[3][4][5] << endl;
	for( unsigned int i=1; i<=data.X.size(); i++)
		for( unsigned int j=1; j<=data.Y.size(); j++ )
			for( unsigned int k=1; k<=data.Z.size(); k++ ){
				double s1=data.X[i-1]*data.Y[j-1]*data.Z[k-1];
				double s0=sm[i-1][j-1][k-1];
				int t=-std::numeric_limits<int>::infinity();
				if(s0>0)
					t=(s1>0)?0x0111:0x0000; //only move forward when strict larger than 0
				else if(s0<0)
					t=(s1<0)?0x0111:0x0000;
				else
					t=0x0000;
				sm[i][j][k]=(t==0x0000)?s1:s0+s1;
				tm[i][j][k]=t;
				//cout << "sm["<<i<<","<<j<<","<<k<<"]="<<sm[i][j][k]<<";s1="<< s1 << ";s0=" << s0 << ";t="<<hex<<tm[i][j][k]<<endl;
				if ( abs(sm[i][j][k]) > max_s ) {
					max_p[0] = i; max_p[1] = j; max_p[2] = k;
					max_s = sm[i][j][k];
				}
			}
	//calculation of the score matrix
	//cout << "max_s=" << max_s << ";max_p="<< max_p[0] <<","<<max_p[1]<<","<<max_p[2]<<endl;
	lla_result.score=sm[max_p[0]][max_p[1]][max_p[2]];
	for( int i=max_p[0], j=max_p[1], k=max_p[2]; ;
				i-=(int)((tm[i][j][k]&0x0001)>0), j-=(int)((tm[i][j][k]&0x0010)>0), k-=(int)((tm[i][j][k]&0x0100)>0) ){
		vector<int> step;
		step.push_back(i); step.push_back(j); step.push_back(k);
		lla_result.trace.push_back(step);
		if (tm[i][j][k] == 0x0000) break;
	}
	//cout << "length="<<lla_result.trace.size()<< endl;
	return lla_result;
}

LA_Result ST_la( const LLA_Data& la_data ) {
  LA_Result la_result;
  double score=0;
  for( unsigned int i=0; i<la_data.X.size(); i++ ){
    score=score+la_data.X[i]*la_data.Y[i]*la_data.Z[i];
  }
  la_result.score=score;
  return la_result;
}

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
			  //		t=(s1>0)?0x0011:0x0000; //only move forward when strict larger than 0
				//else if(s0<0)
			  //		t=(s1<0)?0x0011:0x0000;
				//else
				//	t=0x0000;
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
  //				i-=(int)((tm[i][j]&0x0001)>0), j-=(int)((tm[i][j]&0x0010)>0) ){
  //		vector<int> step;
	//	step.push_back(i); step.push_back(j);
	//	lsa_result.trace.push_back(step);
	//	if (tm[i][j] == 0x0000) break;
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



