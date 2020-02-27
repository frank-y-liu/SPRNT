#!/bin/sh
if [ $# -lt 1 ]; then
  echo "Usage: `basename $0` [options] <target.spt>"
  echo "    Available options are '-steadyonly' '-checkonly'"
  echo "   "
  exit -1
fi

IMAGE=frankliu1/sprnt:latest
exec docker run --rm -i --user="$(id -u):$(id -g)" --net=none -v "$PWD":/data "$IMAGE" run_spt.sh "$@"

## docker image is built from scratch
