#LSA_rep_test
lsa_compute ../test/LSA_rep_test.txt ../test/LSA_rep_test.lsa -r 2 -s 4 -d 1
#LSA_rep_test
lsa_compute ../test/LSA_rep_test.txt ../test/LSA_rep_test.lsa -r 2 -s 4 -d 1 -n none
#testna, appearantly PCC of a zero vector is not defined, so nan in output
lsa_compute ../test/testna.txt ../test/testna.lsa -r 1 -s 4 -d 1 
#olddata
lsa_compute ../test/olddata.txt ../test/olddata.lsa -r 1 -s 35 -d 3
#query_olddata
lsa_query ../test/olddata.lsa ../test/olddata.lsa.entry -x ../test/olddata.lsa.xgmml
#newdata
#lsa_compute.py ../test/newdata.txt ../test/newdata.lsa -r 1 -s 35 -d 3
#debug
lsa_query ../debug/debug2.lsa ../debug/debug2.lsa.entry -x ../debug2.lsa.xgmml
