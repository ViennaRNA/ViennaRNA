$!
$!VAX-VMS cc make file for readseq
$!
$ echo := write sys$output
$ if p1.eqs."TEST" then goto tests
$
$ echo "compiling readseq..."
$ cc readseq, ureadseq
$!
$ echo "linking readseq..."
$ link readseq, ureadseq, sys$library:vaxcrtl/lib
$!
$tests:
$!
$ echo "defining readseq symbol:"
$ dd = f$environment("default")
$ readseq :== $ 'dd'readseq.exe
$ show symbol readseq
$!
$ echo ""
$ echo "test for general read/write of all chars:"
$ readseq -p alphabet.std -otest.alpha
$ diff test.alpha alphabet.std
$!
$ echo ""
$ echo "test for valid format conversions"
$!
$ readseq -v -p -f=ig   nucleic.std -otest.ig
$ readseq -v -p -f=gb   test.ig     -otest.gb
$ readseq -v -p -f=nbrf test.gb     -otest.nbrf
$ readseq -v -p -f=embl test.nbrf   -otest.embl
$ readseq -v -p -f=gcg  test.embl   -otest.gcg
$ readseq -v -p -f=strider test.gcg -otest.strider
$ readseq -v -p -f=fitch test.strider -otest.fitch
$ readseq -v -p -f=fasta test.fitch -otest.fasta
$ readseq -v -p -f=pir  test.fasta  -otest.pir
$ readseq -v -p -f=ig   test.pir    -otest.ig-b
$ diff test.ig test.ig-b
$!
$ echo ""
$ echo "Test for multiple-sequence format conversions:"
$ readseq -p -f=ig    multi.std   -otest.m-ig
$ readseq -p -f=gb    test.m-ig   -otest.m-gb
$ readseq -p -f=nbrf  test.m-gb   -otest.m-nbrf
$ readseq -p -f=embl  test.m-nbrf -otest.m-embl
$ readseq -p -f=fasta test.m-embl -otest.m-fasta
$ readseq -p -f=pir   test.m-fasta -otest.m-pir
$ readseq -p -f=msf   test.m-pir  -otest.m-msf
$ readseq -p -f=paup  test.m-msf  -otest.m-paup
$ readseq -p -f=ig    test.m-paup -otest.m-ig-b
$ diff test.m-ig test.m-ig-b
$ echo ""
$ echo "Expect differences in the header lines due to"
$ echo "different format headers.  If any sequence lines"
$ echo "differ, or if checksums differ, there is a problem."
$!
$! #cleanup
$! delete test.*;
$ echo "-----------"
$ echo ""
$ echo "To clean up test files, command me:
$ echo "  DELETE test.*;"
$!
