echo "Testing RNAalifold (equilibrium probabilities):"

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
RNAalifold -q --noPS -p --MEA -r --auto-id --id-prefix="rnaalifold_pf_test" ${DATADIR}/rfam_seed_selected.stk > rnaalifold_pf.out
diff=$(${DIFF} -I frequency ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.pf.gold rnaalifold_pf.out)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Partition function (base pair probabilities)"
for file in rnaalifold_pf_test_00*_dp.ps rnaalifold_pf_test_00*_ali.out
do
  diff=$(${DIFF} -I CreationDate -I Creator ${RNAALIFOLD_RESULTSDIR}/${file} ${file})
  if [ "x${diff}" != "x" ] ; then break; fi
done
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# clean up
rm rnaalifold_pf.out rnaalifold_pf_test_00*_dp.ps rnaalifold_pf_test_00*_ali.out

exit ${RETURN}
