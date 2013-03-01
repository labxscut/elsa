# Current Tests
lsa_compute ../test/ARISA20.csv ../test/ARISA20.theo.lsa -r 1 -s 127 -d 3 -p theo -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 #72.951840
lsa_compute ../test/ARISA20.csv ../test/ARISA20.mix.lsa -r 1 -s 127 -d 3 -p mix -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 #3501.449665
lsa_compute ../test/ARISA20.csv ../test/ARISA20.perm.lsa -r 1 -s 127 -d 3 -p perm -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 #14737.259068

lsa_compute ../test/ARISA20.csv ../test/ARISA20.lsa -r 1 -s 127 -d 3 -p theo -x 1000 -f linear -n percentileZ -e ../test/ARISA20.csv -m 0
lsa_query ARISA20.lsa ARISA20.sig.lsaq -q '(!lsa$P>0.01)&(lsa$Q<0.01)' -x ARISA20.sig.xgmml -s ARISA20.sig.sif
la_compute ARISA20.csv ARISA20.sig.lsaq ARISA20.sig.la -s 127 -n percentileZ -m 0 -p perm -x 1000 -f linear
la_query ARISA20.sig.lsaq ARISA20.sig.la ARISA20.sig.laq -q '(!la$P>0.01)&(la$Q<0.01)' -x ARISA20.sig.xgmml

#### Local Trend Tests ####
to_trend -T 1 -S 119 ARISA20.txt ARISA20.trend.txt

#Outdated Tests
#LSA test na, appearantly PCC of a zero vector is not defined, so shown nan in output
#lsa_compute ../test/testna.txt ../test/testna.lsa -r 2 -s 4 -d 0 
#la_compute testna.txt testna.lsa testna.la -s 4 -r 2

#
#lsa_compute 090412_sv.txt 090412_sv.lsa -s 40 -r 1 -d 3 -f linear -p perm -x 2000 -m 0 -n percentile
#lsa_compute 090412_z.txt 090412_z.lsa -s 40 -r 1 -d 3 -f linear -p perm -x 2000 -m 0 -n percentile

#LSA test replicates simple
#lsa_compute ../test/testrep.txt ../test/testrep.simple.lsa -r 5 -s 20 -d 3 -n pnz -m 0 -x 1000 -p mix
#LSA test replicates SD
#lsa_compute ../test/testrep.txt ../test/testrep.SD.lsa -r 5 -s 20 -d 3 -n none -t SD
#LSA test replicates Med
#lsa_compute ../test/testrep.txt ../test/testrep.Med.lsa -r 5 -s 20 -d 3 -n none -t Med
#LSA test replicates MAD
#lsa_compute ../test/testrep.txt ../test/testrep.MAD.lsa -r 5 -s 20 -d 3 -n none -t MAD
#Some perculir data
#lsa_compute ../test/PreLSAspl3na.txt ../test/PreLSAspl3na.lsa -r 1 -s 43 -d 1 -p 100
#Some real olddata
#lsa_compute ../test/olddata.txt ../test/olddata.lsa -r 1 -s 35 -d 3
#lsa_compute ../test/olddata.txt ../test/olddata.1.lsa -r 1 -s 35 -d 3 -f linear -p perm -x 100 -m 0 -n percentile &
#
#lsa_compute ../test/ShortJAC.csv ../test/ShortJAC.lsa -r 1 -s 127 -d 3 -p -100
#
#lsa_compute ../test/JAC0503.txt ../test/JAC0503.lsa -r 1 -s 127 -d 3 -p -100000 -f none


#query_olddata
#lsa_query ../test/olddata.lsa ../test/olddata.lsa.entry -x ../test/olddata.lsa.xgmml
#newdata
#lsa_compute.py ../test/newdata.txt ../test/newdata.lsa -r 1 -s 35 -d 3
#debug
#lsa_query ../debug/debug2.lsa ../debug/debug2.lsa.entry -x ../debug2.lsa.xgmml
