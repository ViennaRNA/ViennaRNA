echo "Testing RNAcofold (equilibrium probabilities):"

RETURN=0

function failed {
    RETURN=1
    echo " [ NOT OK ]"
}

function passed {
    echo " [ OK ]"
}

function testline {
  echo -en "...testing $1:\t\t"
}

# Test partition function folding (centroid, MEA, base pair- and stack probabilities)
testline "Partition function (centroid, MEA)"
RNAcofold --noPS -p --centroid --MEA < ${DATADIR}/rnacofold.small.seq > rnacofold_pf.fold
diff=$(${DIFF} -I frequency ${RNACOFOLD_RESULTSDIR}/rnacofold.small.pf.gold rnacofold_pf.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Partition function (centroid, dimers, monomers)"
RNAcofold --noPS -a < ${DATADIR}/rnacofold.small.seq > rnacofold_pf.fold
diff=$(${DIFF} -I frequency ${RNACOFOLD_RESULTSDIR}/rnacofold.small.pf.full.gold rnacofold_pf.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Partition function (base pair probabilities - Hetero-dimer)"
RNAcofold --noPS -a --auto-id --id-prefix="rnacofold_pf_test" < ${DATADIR}/rnacofold.small.seq > rnacofold_pf.fold
for file in ABrnacofold_pf_test_00*dp5.ps
do
  diff=$(${DIFF} -I CreationDate -I Creator ${RNACOFOLD_RESULTSDIR}/${file} ${file})
  if [ "x${diff}" != "x" ] ; then break; fi
done
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Partition function (base pair probabilities - Homo-dimers)"
for file in AArnacofold_pf_test_00*dp5.ps BBrnacofold_pf_test_00*dp5.ps
do
  diff=$(${DIFF} -I CreationDate -I Creator ${RNACOFOLD_RESULTSDIR}/${file} ${file})
  if [ "x${diff}" != "x" ] ; then break; fi
done
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Partition function (base pair probabilities - Monomers)"
for file in Arnacofold_pf_test_00*dp5.ps Brnacofold_pf_test_00*dp5.ps
do
  diff=$(${DIFF} -I CreationDate -I Creator ${RNACOFOLD_RESULTSDIR}/${file} ${file})
  if [ "x${diff}" != "x" ] ; then break; fi
done
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test concentrations
testline "Partition function (concentrations)"
RNAcofold --noPS -f ${DATADIR}/rnacofold.concentrations < ${DATADIR}/rnacofold.small.seq > rnacofold_concentrations.fold
diff=$(${DIFF} -I frequency ${RNACOFOLD_RESULTSDIR}/rnacofold.small.concentrations.gold rnacofold_concentrations.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# clean up
rm rnacofold_pf.fold ABdot5.ps AAdot5.ps BBdot5.ps Adot5.ps Bdot5.ps rnacofold_concentrations.fold
rm ABrnacofold_pf_test_00*dp5.ps AArnacofold_pf_test_00*dp5.ps BBrnacofold_pf_test_00*dp5.ps Arnacofold_pf_test_00*dp5.ps Brnacofold_pf_test_00*dp5.ps dot.ps

exit ${RETURN}
