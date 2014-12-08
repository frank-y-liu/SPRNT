#! /bin/bash

# setup the library search path and run 
SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

if [ ! -r "$SCRIPTPATH/sprnt" ]; then
  echo "Need the SPRNT binary to continue ..."
  exit -1
fi
if [ ! -r "$SCRIPTPATH/../lib/libsolvers.so" ]; then
  echo "Need the solver library to proceed ..."
  exit -1 
fi
if [ $# -lt 1 ]; then
  $SCRIPTPATH/sprnt
  exit 1
fi

LD_LIBRARY_PATH=$SCRIPTPATH/../lib:$LD_LIBRARY_PATH $SCRIPTPATH/sprnt "$@"

## 
