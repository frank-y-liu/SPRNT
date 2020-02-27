#!/bin/sh
if [ $# -lt 1 ]; then
  echo "Usage: `basename $0` <COMMNAND> <target.spt>"
  echo "    Available COMMAND are 'run_sprnt.sh' 'chk_outputs.sh' 'extract_node.sh' 'extract_time.sh'"
  echo "     Example: `basename $0` run_sprnt.sh m1.spt"
  echo "   "
  exit -1
fi

IMAGE=frankliu1/sprnt:latest

exec docker run --rm -i --user="$(id -u):$(id -g)" --net=none -v "$PWD":/data "$IMAGE" "$@"

## docker image is built from scratch
