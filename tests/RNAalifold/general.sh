echo "Testing RNAalifold (MFE and general features):"

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
rnaalifold_version=$(RNAalifold --version)
rnaalifold_version=${rnaalifold_version/$CURRENT_VERSION/LATEST}
if [ "x${rnaalifold_version}" != "xRNAalifold LATEST" ] ; then failed; else passed; fi

# Test dangle models (d0, d2)
for dangles in 0 2
do
  testline "MFE prediction (RNAalifold -d${dangles})"
  RNAalifold -q --noPS -d ${dangles} ${DATADIR}/rfam_seed_many_short.stk > rnaalifold.out
  diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_many_short.d${dangles}.mfe.gold rnaalifold.out)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test --noLP
testline "MFE prediction (RNAalifold --noLP)"
RNAalifold -q --noPS --noLP ${DATADIR}/rfam_seed_selected.stk > rnaalifold.out
diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.noLP.mfe.gold rnaalifold.out)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test temperatures (25degC, 40degC)
for T in 25 40
do
  testline "MFE prediction (RNAalifold -T ${T})"
  RNAalifold -q --noPS -T ${T} ${DATADIR}/rfam_seed_selected.stk > rnaalifold.out
  diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.T${T}.mfe.gold rnaalifold.out)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test maximum base pair span
testline "MFE prediction (RNAalifold --maxBPspan=30)"
RNAalifold -q --noPS --maxBPspan=30 -p0 ${DATADIR}/rfam_seed_selected.stk > rnaalifold.out
diff=$(${DIFF} -I frequency ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.span30.gold rnaalifold.out)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test automatic ID generation
testline "MFE prediction and automatic alignment ID feature"
RNAalifold -q --noPS --auto-id --id-prefix="blabla" --id-start=40 ${DATADIR}/rfam_seed_selected.stk > rnaalifold.out
diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.id.mfe.gold rnaalifold.out)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test parameter file loading
for par in dna_mathews1999.par dna_mathews2004.par rna_andronescu2007.par rna_turner1999.par
do
  testline "MFE prediction (-P ${par})"
  RNAalifold -q --noPS -P ${MISC_DIR}/${par} ${DATADIR}/rfam_seed_selected.stk > rnaalifold.out
  diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.${par}.mfe.gold rnaalifold.out)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test ribosum scoring
testline "MFE prediction with ribosum scoring (RNAalifold -r)"
RNAalifold -q -r --noPS ${DATADIR}/rfam_seed_selected.stk > rnaalifold.out
diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_selected.ribosum.mfe.gold rnaalifold.out)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test different input formats (clustal, fasta, maf, stockholm)
for file in alignment_clustal.aln alignment_fasta.fa alignment_maf.maf alignment_stockholm.stk
do
  testline "MFE prediction with different file formats (${file})"
  RNAalifold -q --noPS -r ${DATADIR}/${file} > rnaalifold.out
  diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/${file}.mfe.gold rnaalifold.out)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# clean up
rm rnaalifold.out

exit ${RETURN}
