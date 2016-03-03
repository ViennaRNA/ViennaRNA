#!/bin/bash
#

PROGNAME=`basename $1 .ggo`
PROGPATH=`dirname $1`

help2man -N --help-option=--detailed-help \
  --opt-include=./include/ref_package.inc \
                --opt-include=./include/${PROGNAME}.inc \
                --opt-include=./include/ref_energy_par.inc \
                "./cmdlopt.sh ${PROGPATH}/${PROGNAME}"
