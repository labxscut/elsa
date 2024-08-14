#activate elsa virtual v
conda activate elsa
python -c "from lsa import lsa_version; lsa_version.main()"

# Current LSA Tests
echo "ELSA Tests"
lsa_chkdat ../test/ARISA20.csv 1 127 # no output is expected if no data errors found. 
lsa_compute ../test/ARISA20.csv ../test/ARISA20.theo.lsa -r 1 -s 127 -d 3 -p theo -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 #72.951840
#lsa_compute ../test/ARISA20.csv ../test/ARISA20.mix.lsa -r 1 -s 127 -d 3 -p mix -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 #3501.449665
#lsa_compute ../test/ARISA20.csv ../test/ARISA20.perm.lsa -r 1 -s 127 -d 3 -p perm -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 #14737.259068

# Current LTA Tests 
echo "LTA Tests"
lsa_totrend -T 1 -S 119 ARISA20.txt ARISA20.trend.txt
lsa_compute ../test/ARISA20.csv ../test/ARISA20.theo.lta -r 1 -s 127 -d 3 -p theo -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 -a 0.5549 -T 1 #4.043180s
#lsa_compute ../test/ARISA20.csv ../test/ARISA20.mix.lta -r 1 -s 127 -d 3 -p mix -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 -a 0.5549 -T 1 #
#lsa_compute ../test/ARISA20.csv ../test/ARISA20.perm.lsa -r 1 -s 127 -d 3 -p perm -x 1000 -f none -n percentileZ -e ../test/ARISA20.csv -m 0 -a 0.5549 -T 1 #

# Current LA Tests
#echo "ELA Tests"
#lsa_compute ../test/ARISA20.csv ../test/ARISA20.lsa -r 1 -s 127 -d 3 -p theo -x 1000 -f linear -n percentileZ -e ../test/ARISA20.csv -m 0
#lsa_query ARISA20.lsa ARISA20.sig.lsaq -q '(!lsa$P>0.01)&(lsa$Q<0.01)' -x ARISA20.sig.xgmml -s ARISA20.sig.sif
#la_compute ARISA20.csv ARISA20.sig.lsaq ARISA20.sig.la -s 127 -n percentileZ -m 0 -p perm -x 1000 -f linear
#la_query ARISA20.sig.lsaq ARISA20.sig.la ARISA20.sig.laq -q '(!la$P>0.01)&(la$Q<0.01)' -x ARISA20.sig.xgmml


# Parallal (require user customization based on  hpc system, not supported) 
#echo "par_ana Tests"
#par_ana ARISA20.txt ARISA20.lsa 'lsa_compute %s %s -e ARISA20.txt -s 127 -r 1 -p theo' $PWD
#par_ana ARISA20.txt ARISA20.la 'la_compute %s ARISA20.laq %s -s 127 -r 1 -p 1000' $PWD

#activate elsa virtual v
conda deactivate
