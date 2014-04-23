#! /bin/bash

# check the content of the output file, and print to stdout
# need sed, awk, sort etc

#
# Frank Liu, 02/2014
#

if [ $# -lt 1 ]; then
  echo "Usage: `basename $0` <sprnt_output_file>"
  exit 1
fi

if [ ! -e "$1" ]; then
  echo "bummer: unable to find the file $1 ..."
  exit -1
fi

grep -q "SPRNt Results" $1
if [ ! $? -eq 0 ]; then
  echo "bummer: input file \"$1\" does seem like a SPRNT output ..."
  exit -1
fi

echo "Output file $1 contains the following nodes:"
cat $1 | sed '/^\*\*\*/d' | awk '{print $1}' | sort -n -u |\
awk '{printf("%s ",$1); if ( NR % 13 == 0 ){printf("\n");}}END{printf("\n");}'
echo " "
echo "Output file $1 contains unsteady results at the following time points:"
cat $1 | sed '/^\*\*\*/d' | awk '{print $2}' | sort -n -u |\
awk '{printf("%s ",$1); if ( NR % 13 == 0 ){printf("\n");}}END{printf("\n");}'
echo " "
echo "It contains the following fields:"
head -n 2 $1 | sed 's/\*//g'
echo " "

# Local Variables:
# mode: shell-script
# End:
