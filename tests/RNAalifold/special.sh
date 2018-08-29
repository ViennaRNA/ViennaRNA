echo "Testing RNAalifold (special features):"

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
RNAalifold -q --noPS -g --auto-id --id-prefix="rnaalifold_gquad_test" -p --MEA \
  ${DATADIR}/adam10_5utr.aln \
  ${DATADIR}/bcl-2_5utr.aln \
  ${DATADIR}/crem_5utr.aln \
  ${DATADIR}/ebag9_5utr.aln \
  ${DATADIR}/ncam2_5utr.aln > rnaalifold_special.fold
diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rnaalifold.gquad.gold rnaalifold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test parallel processing support
testline "Parallel processing of input data"
RNAalifold -q --noPS -j5 ${DATADIR}/rfam_seed_many_short.stk > rnaalifold_special.fold
diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rfam_seed_many_short.d2.mfe.gold rnaalifold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi

# Test SHAPE reactivity input (Deigan method)
for f in 070313_ecoli_cdiff_16S_clustalw.aln
do
  testline "SHAPE probing data ($f) Deigan et al. 2009 method"
  RNAalifold -q --noPS --shape=2=${DATADIR}/cdiff_16S_1m7.shape,1=${DATADIR}/ecoli_16S_1m7.shape --shapeMethod=D ${DATADIR}/${f} > rnaalifold_special.fold
  diff=$(${DIFF} ${RNAALIFOLD_RESULTSDIR}/rnaalifold.SHAPE.gold rnaalifold_special.fold)
  if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi
done

# Test circfold
testline "MFE/PF prediction (RNAalifold -p0 --circ)"
RNAalifold -q --noPS -c -p0 ${DATADIR}/rfam_seed_many_short.stk > rnaalifold_special.fold
diff=$(${DIFF} -I frequency ${RNAALIFOLD_RESULTSDIR}/rfam_seed_many_short.circ.gold rnaalifold_special.fold)
if [ "x${diff}" != "x" ] ; then failed; echo -e "$diff"; else passed; fi


# clean up
rm rnaalifold_special.fold rnaalifold_gquad*_dp.ps rnaalifold_gquad*_ali.out

exit ${RETURN}
