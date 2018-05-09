#!/bin/bash
#

HELP2MAN=$1
PROGNAME=`basename $2 .ggo`
PROGPATH=`dirname $2`

${HELP2MAN} \
    -N \
    --help-option=--detailed-help \
    --opt-include=./include/ref_package.inc \
    --opt-include=./include/${PROGNAME}.inc \
    --opt-include=./include/ref_energy_par.inc \
    "./cmdlopt.sh ${PROGPATH}/${PROGNAME}"
