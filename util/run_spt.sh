#! /bin/bash

myrealpath() {
  OURPWD=$PWD
  cd "$(dirname "$1")"
  LINK=$(readlink "$(basename "$1")")
  while [ "$LINK" ]; do
    cd "$(dirname "$LINK")"
    LINK=$(readlink "$(basename "$1")")
  done
  REALPATH="$PWD/$(basename "$1")"
  cd "$OURPWD"
  echo "$REALPATH"
}

# setup the library search path and run 
SCRIPT=`myrealpath $0`
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
