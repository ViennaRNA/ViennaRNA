echo "Testing RNAfold (MFE and general features):"

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
rnafold_version=$(RNAfold --version)
rnafold_version=${rnafold_version/$CURRENT_VERSION/LATEST}
if [ "x${rnafold_version}" != "xRNAfold LATEST" ] ; then failed; else passed; fi

# Test dangle models (d0, d1, d2, d3)
for dangles in 0 1 2 3
do
  testline "MFE prediction (RNAfold -d${dangles})"
  RNAfold --noPS -d ${dangles} < ${DATADIR}/rnafold.small.seq > rnafold.fold
  diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.small.d${dangles}.mfe.gold rnafold.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test --noLP
testline "MFE prediction (RNAfold --noLP)"
RNAfold --noPS --noLP < ${DATADIR}/rnafold.small.seq > rnafold.fold
diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.small.noLP.mfe.gold rnafold.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test temperatures (25degC, 40degC)
for T in 25 40
do
  testline "MFE prediction (RNAfold -T ${T})"
  RNAfold --noPS -T ${T} < ${DATADIR}/rnafold.small.seq > rnafold.fold
  diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.small.T${T}.mfe.gold rnafold.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test maximum base pair span
testline "MFE prediction (RNAfold --maxBPspan=30)"
RNAfold --noPS --maxBPspan=30 -p0 < ${DATADIR}/rnafold.small.seq > rnafold.fold
diff=$(${DIFF} -I frequency ${RNAFOLD_RESULTSDIR}/rnafold.small.span30.gold rnafold.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test automatic ID generation
testline "MFE prediction and automatic sequence ID feature"
RNAfold --noPS --auto-id --id-prefix="blabla" --id-start=40 < ${DATADIR}/rnafold.small.seq > rnafold.fold
diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.small.id.mfe.gold rnafold.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test parameter file loading
for par in dna_mathews1999.par dna_mathews2004.par rna_andronescu2007.par rna_turner1999.par
do
  testline "MFE prediction (-P ${par})"
  RNAfold --noPS -P ${MISC_DIR}/${par} < ${DATADIR}/rnafold.small.seq > rnafold.fold
  diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.small.${par}.mfe.gold rnafold.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# clean up
rm rnafold.fold

exit ${RETURN}
