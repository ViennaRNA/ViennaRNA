echo "Testing RNAfold (longer sequences):"

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

# MFE and input without FASTA header
testline "MFE prediction (.seq file)"
RNAfold --noPS < ${DATADIR}/rnafold.seq > tmp.fold
diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.seq.mfe.gold tmp.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# MFE and input with FASTA header
testline "MFE prediction (.fasta file)"
RNAfold --noPS < ${DATADIR}/rnafold.fasta > tmp.fold
diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.fasta.mfe.gold tmp.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# clean up
rm tmp.fold

exit ${RETURN}
