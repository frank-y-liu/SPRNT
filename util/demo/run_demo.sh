#! /bin/bash

# run demo, but first add search path
SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
if [ ! -r "$SCRIPTPATH/../lib/libsolvers.so" ]; then
  echo "Need the solver library to proceed ..."
  exit -1 
fi

LD_LIBRARY_PATH=$SCRIPTPATH/../lib:$LD_LIBRARY_PATH $SCRIPTPATH/api_demo
