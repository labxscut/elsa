#LSA_rep_test
lsa_compute.py ../test/LSA_rep_test.txt ../test/LSA_rep_test.lsa -r 2 -s 4 -d 1
#testna, appearantly PCC of a zero vector is not defined, so nan in output
lsa_compute.py ../test/testna.txt ../test/testna.lsa -r 1 -s 4 -d 1
#olddata
lsa_compute.py ../test/olddata.txt ../test/olddata.lsa -r 1 -s 35 -d 3
#query_olddata
lsa_query.py ../test/olddata.lsa ../test/olddata.lsa.entry -s ../test/olddata.lsa.sif
#newdata
#lsa_compute.py ../test/newdata.txt ../test/newdata.lsa -r 1 -s 35 -d 3
