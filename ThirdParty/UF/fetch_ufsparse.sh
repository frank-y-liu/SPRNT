#! /bin/bash

# fetch SuiteSparse from Texas A&M CSE webpage
# 
# NOTE: the exact locations of the files could be changed without notice
#       modify this script if necessary!
# 
#  The repository at UF is gone

echo " "
echo "Fetching cached SuiteSparse .... "
echo " "

if [ ! -f SuiteSparse-3.0.0.tar.gz ]; then
## original location
##curl http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-3.0.0.tar.gz -o SuiteSparse-3.0.0.tar.gz
cp ../Repositories/UF/SuiteSparse-3.0.0.tar.gz ./.
fi

echo " "
echo "Expanding the source code ......"
echo " "
tar -xzvf SuiteSparse-3.0.0.tar.gz  > /dev/null 2>&1

# the original UFconfig.mk doesn't quite work for me
# use the hacked version
# 

echo " "
echo "Updating the config file ......"
echo " "
\mv -f SuiteSparse/UFconfig/UFconfig.mk SuiteSparse/UFconfig/UFconfig.mk.org
\cp -f ../../UFconfig.mk SuiteSparse/UFconfig/UFconfig.mk
\cp -f ../../uf_makefile.local ./makefile.local

echo "... Done!"
