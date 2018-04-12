#! /bin/bash

#
# shell script to extract unsteady results of ALL nodes at a particular time point
#  
#  see also: extract_node.sh 
#            chk_output.sh
#
#  if xgraph is installed, you can quickly pipe the results to xgrpah, e.g.:
#
#    ./extract_time.sh -f <file name> -t <time> -q | xgraph -bb -P
#
#  or save to a text file and use gnuplot
# 
#  needde:
#      awk, head, sed etc
# 

# 
# Frank Liu, 02/2014
#

find_loc_flow () 
{
   NUM=0
   for i in "$@"; do
       (( NUM += 1 ))
       if [[ "$i" = flow* ]]; then
	   return $NUM
       fi
   done
   return 0
}

find_loc_area () 
{
   NUM=0
   for i in "$@"; do
       (( NUM += 1 ))
       if [[ "$i" = wet_a* ]]; then
	   return $NUM
       fi
   done
   return 0
}

find_loc_depth () 
{
   NUM=0
   for i in "$@"; do
       (( NUM += 1 ))
       if [[ "$i" = depth* ]]; then
	   return $NUM
       fi
   done
   return 0
}

find_loc_surf () 
{
   NUM=0
   for i in "$@"; do
       (( NUM += 1 ))
       if [[ "$i" = surf* ]]; then
	   return $NUM
       fi
   done
   return 0
}

if [ $# -lt 2 ]; then
    echo "Usage: `basename $0` -f <output file> -t <time> [-q|-a|-d|-z]"
    echo "       `basename $0` --file <output file> --time <time> [-flow|-area|-depth|-surf]"
    echo " "
    echo "    extract results from sprnt output file at a particular time point"
    echo "    one of the following options can be used:"
    echo "    -a or --area  : plot wetted area"
    echo "    -q or --flow  : plot flow"
    echo "    -d or --depth : plot depth"
    echo "    -z or --surf  : plot surface elevation"
    echo "  if more than one options are specified, only the last one will be honored"
    echo "  "
    echo "  specifying no option will print all available fields"
    exit -1;
fi

PRT_A=0
PRT_Q=0
PRT_Z=0
PRT_D=0
while [[ $# > 0 ]]; do
    key="$1"
    shift

    case $key in
	-f|--file)
	FNAME="$1"
	shift
	;;
	-t|--time)
	TIME="$1"
	shift
	;;
	-a|--area)
	PRT_A=1
	;;
	-q|--flow)
	PRT_Q=1
	;;
	-d|--depth)
	PRT_D=1
	;;
	-z|--surf)
	PRT_Z=1
	;;
	*)

	;;
    esac
done

if [ ! -e "$FNAME" ]; then
    echo "Bummer: unable to find file named $FNAME ..."
    exit -2
fi

grep -q "SPRNt Results" $FNAME 
if [ ! $? -eq 0 ]; then
  echo "bummer: input file \"$FNAME\" does seem like a SPRNT output ..."
  exit -1
fi

if [ "X$TIME" = "X" ]; then
    echo "Need to specify the time ..."
    exit -2
fi


echo "###### File: $FNAME Time: $TIME"

(( TT = PRT_A + PRT_Q + PRT_D + PRT_Z ))
if [ $TT -eq 0 ]; then
    awk -v tt=$TIME -v pp=$ARG '
{ 
if (index($1, "*")==0 && $2==tt) {
  for (i=1;i<=NF;i++) {
    if (i!=2) {
      printf("%s ",$i);
    }
  }
  printf("\n");
}
}' $FNAME
    exit 1
fi

STR=`head -n 2 $FNAME | tail -n 1 | sed 's/\*//g'`
find_loc_flow $STR
HAS_Q=$?
find_loc_area $STR
HAS_A=$?
find_loc_depth $STR
HAS_D=$?
find_loc_surf $STR
HAS_Z=$?

if [ $PRT_Q -gt 0 ]; then
    let ARG=$HAS_Q
fi
if [ $PRT_A -gt 0 ]; then
    let ARG=$HAS_A
fi
if [ $PRT_D -gt 0 ]; then
    let ARG=$HAS_D
fi
if [ $PRT_Z -gt 0 ]; then
    let ARG=$HAS_Z
fi

if [ $ARG -eq 0 ]; then
    echo "Bummer: cannot find the data field specified ...."
    exit -2
fi

awk -v tt=$TIME -v pp=$ARG '
{ 
if (index($1, "*")==0 && $2==tt) {
   print $1,$pp
}
}' $FNAME

# Local Variables:
# mode: shell-script
# End:
