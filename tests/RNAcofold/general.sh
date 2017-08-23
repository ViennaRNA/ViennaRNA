echo "Testing RNAcofold (MFE and general features):"

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

# simple version test
testline "correct version number"
rnacofold_version=$(RNAcofold --version)
rnacofold_version=${rnacofold_version/$CURRENT_VERSION/LATEST}
if [ "x${rnacofold_version}" != "xRNAcofold LATEST" ] ; then failed; else passed; fi

# Test dangle models (d0, d1, d2, d3)
for dangles in 0 1 2 3
do
  testline "MFE prediction (RNAcofold -d${dangles})"
  RNAcofold --noPS -d ${dangles} < ${DATADIR}/rnacofold.small.seq > rnacofold.fold
  diff=$(${DIFF} ${RNACOFOLD_RESULTSDIR}/rnacofold.small.d${dangles}.mfe.gold rnacofold.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test --noLP
testline "MFE prediction (RNAcofold --noLP)"
RNAcofold --noPS --noLP < ${DATADIR}/rnacofold.small.seq > rnacofold.fold
diff=$(${DIFF} ${RNACOFOLD_RESULTSDIR}/rnacofold.small.noLP.mfe.gold rnacofold.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test temperatures (25degC, 40degC)
for T in 25 40
do
  testline "MFE prediction (RNAcofold -T ${T})"
  RNAcofold --noPS -T ${T} < ${DATADIR}/rnacofold.small.seq > rnacofold.fold
  diff=$(${DIFF} ${RNACOFOLD_RESULTSDIR}/rnacofold.small.T${T}.mfe.gold rnacofold.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test maximum base pair span
testline "MFE prediction (RNAcofold --maxBPspan=30)"
RNAcofold --noPS --maxBPspan=30 < ${DATADIR}/rnacofold.small.seq > rnacofold.fold
diff=$(${DIFF} ${RNACOFOLD_RESULTSDIR}/rnacofold.small.span30.gold rnacofold.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test automatic ID generation
testline "MFE prediction and automatic sequence ID feature"
RNAcofold --noPS --auto-id --id-prefix="ballaballa" --id-start=40 < ${DATADIR}/rnacofold.small.seq > rnacofold.fold
diff=$(${DIFF} ${RNACOFOLD_RESULTSDIR}/rnacofold.small.id.mfe.gold rnacofold.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test parameter file loading
for par in dna_mathews1999.par dna_mathews2004.par rna_andronescu2007.par rna_turner1999.par
do
  testline "MFE prediction (-P ${par})"
  RNAcofold --noPS -P ${MISC_DIR}/${par} < ${DATADIR}/rnacofold.small.seq > rnacofold.fold
  diff=$(${DIFF} ${RNACOFOLD_RESULTSDIR}/rnacofold.small.${par}.mfe.gold rnacofold.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# clean up
rm rnacofold.fold

exit ${RETURN}
