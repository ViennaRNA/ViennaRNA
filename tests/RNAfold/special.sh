echo "Testing RNAfold (special features):"

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

# Test G-Quadruplex feature
testline "G-Quadruplex feature (MFE, centroid, MEA)"
RNAfold --noPS -p --MEA -g --auto-id --id-prefix="rnafold_gquad_test" < ${DATADIR}/rnafold.gquad.fa > rnafold_special.fold
diff=$(${DIFF} -I frequency ${RNAFOLD_RESULTSDIR}/rnafold.gquad.gold rnafold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "G-Quadruplex feature (probabilities)"
for file in rnafold_gquad_test_00*dp.ps
do
  diff=$(${DIFF} -I CreationDate -I Creator ${RNAFOLD_RESULTSDIR}/${file} ${file})
  if [ "x${diff}" != "x" ] ; then break; fi
done
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test Ligand-motif feature
testline "Aptamer Motif - Theophylline (MFE, centroid, MEA)"
RNAfold --noPS -p --MEA --motif="GAUACCAG&CCCUUGGCAGC,(...((((&)...)))...),-9.22" --id-prefix="rnafold_theo_test" < ${DATADIR}/rnafold.theo.fa > rnafold_special.fold
diff=$(${DIFF} -I frequency ${RNAFOLD_RESULTSDIR}/rnafold.theo.gold rnafold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Aptamer Motif - Theophylline (probabilities)"
diff=$(${DIFF} -I CreationDate -I Creator ${RNAFOLD_RESULTSDIR}/rnafold_theo_test_0001_dp.ps rnafold_theo_test_0001_dp.ps)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test Constraints (Hard, Soft, Unstructured domains)
testline "Command file - Constraints and Ligand motifs (MFE, centroid, MEA)"
RNAfold --noPS --commands=${DATADIR}/rnafold.cmds -v -p --MEA --auto-id --id-prefix="rnafold_cmd_test" < ${DATADIR}/rnafold.small.seq > rnafold_special.fold
diff=$(${DIFF}  -I frequency -I MEA ${RNAFOLD_RESULTSDIR}/rnafold.small.cmds.gold rnafold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

testline "Command file - Constraints and Ligand motifs (probabilities)"
for file in rnafold_cmd_test_00*dp.ps
do
  diff=$(${DIFF} -I CreationDate -I Creator ${RNAFOLD_RESULTSDIR}/${file} ${file})
  if [ "x${diff}" != "x" ] ; then break; fi
done
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi


# Test SHAPE reactivity input (Deigan, Zarringhalam, and Washietl methods)
for f in TPP_riboswitch_E.coli Lysine_riboswitch_T._martima 5domain16S_rRNA_H.volcanii 5domain16S_rRNA_E.coli
do
  testline "SHAPE probing data ($f) Deigan et al. 2009 method"
  RNAfold --noPS --shape=${DATADIR}/${f}.shape_2rows --shapeMethod=D < ${DATADIR}/${f}.db > rnafold_special.fold
  diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.SHAPE.${f}.D.gold rnafold_special.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

  testline "SHAPE probing data ($f) Zarringhalam et al. 2012 method"
  RNAfold --noPS --shape=${DATADIR}/${f}.shape_2rows --shapeMethod=Z < ${DATADIR}/${f}.db > rnafold_special.fold
  diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.SHAPE.${f}.Z.gold rnafold_special.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

  testline "SHAPE probing data ($f) Washietl et al. 2012 method"
  RNAfold --noPS --shape=${DATADIR}/${f}.pvmin --shapeMethod=W < ${DATADIR}/${f}.db > rnafold_special.fold
  diff=$(${DIFF} ${RNAFOLD_RESULTSDIR}/rnafold.SHAPE.${f}.W.gold rnafold_special.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test circfold
testline "MFE prediction (RNAfold -p0 --circ)"
RNAfold --noPS -c -p0 < ${DATADIR}/rnafold.small.seq > rnafold_special.fold
diff=$(${DIFF} -I frequency ${RNAFOLD_RESULTSDIR}/rnafold.small.circ.gold rnafold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# clean up
rm rnafold_special.fold
rm rnafold_gquad_test_00*dp.ps rnafold_cmd_test_00*dp.ps rnafold_theo_test_0001_dp.ps

exit ${RETURN}
