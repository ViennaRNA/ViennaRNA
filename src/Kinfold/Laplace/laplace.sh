#!/bin/sh

# strip leading path and tailing extension(s)
function to_basename {
    local basename=${1%/}
    basename=${basename##*/}
    basename=${basename%%.*}
    echo $basename
}

KINFOLD=../Kinfold
INFILE=$1

OUTFILE=$(to_basename $INFILE)

for phi in `seq 0.3 0.1 8.0`
do
 $KINFOLD \
 --silent --phi $phi --log $OUTFILE.$phi --num 2000 --time 1000000 < $INFILE
done

perl extract_data.pl *.log > $OUTFILE.data

R CMD BATCH to_boxplot.R

# End of file
