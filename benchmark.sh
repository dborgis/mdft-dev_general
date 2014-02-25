red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

for i in 1
    do mkdir -p input output
    cp benchmark/dft.in.$i input/dft.in
    cp benchmark/solute.in.$i input/solute.in
    cp benchmark/solvent.in.$i input/solvent.in
    rm log
    ./mdft &> log
    TESTRES=`grep "TOTAL ENERGY" log | awk '{ print $4 }'`
    OBJECTIVERESULT="16.452397"
    ERROR=`echo "($RES-$OBJECTIVERESULT)/$OBJECTIVERESULT*100)" | bc -l`
    if [ "$ERROR" \< "0.1" ]
    then 
    echo -e "Test $i ${green}OK${NC}"
    else
    echo -e "Test $i ${red}failed${NC}"
    fi
done
