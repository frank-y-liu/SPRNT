#! /bin/bash

# get packages from University of Florida CISE
# 
# NOTE: the exact locations of the files could be changed without notice
#       modify this script if necessary!

echo " "
echo "Fetching UMFPACK and KLU from U. Florida CISE webpage,"
echo " may take a few minutes depending on the speed of the internet connection ...."
echo " "
if [ ! -f UMFPACK-5.1.0.tar.gz ]; then
curl http://www.cise.ufl.edu/research/sparse/umfpack/UMFPACK-5.1.0.tar.gz -o UMFPACK-5.1.0.tar.gz
fi
if [ ! -f AMD-2.2.1.tar.gz ]; then
curl http://www.cise.ufl.edu/research/sparse/amd/AMD-2.2.1.tar.gz -o AMD-2.2.1.tar.gz
fi
if [ ! -f BTF-1.1.3.tar.gz ]; then
curl http://www.cise.ufl.edu/research/sparse/btf/BTF-1.1.3.tar.gz -o BTF-1.1.3.tar.gz
fi
if [ ! -f COLAMD-2.7.4.tar.gz ]; then
curl http://www.cise.ufl.edu/research/sparse/colamd/COLAMD-2.7.4.tar.gz -o COLAMD-2.7.4.tar.gz
fi
if [ ! -f KLU-1.1.4.tar.gz ]; then
curl http://www.cise.ufl.edu/research/sparse/klu/KLU-1.1.4.tar.gz -o KLU-1.1.4.tar.gz
fi
if [ ! -f UFconfig-3.7.1.tar.gz ]; then
curl http://www.cise.ufl.edu/research/sparse/UFconfig/UFconfig-3.7.1.tar.gz -o UFconfig-3.7.1.tar.gz
fi

echo " "
echo "Expanding the source code ......"
echo " "
tar -xzvf UMFPACK-5.1.0.tar.gz  > /dev/null 2>&1
tar -xzvf AMD-2.2.1.tar.gz > /dev/null 2>&1
tar -xzvf BTF-1.1.3.tar.gz > /dev/null 2>&1
tar -xzvf COLAMD-2.7.4.tar.gz  > /dev/null 2>&1 
tar -xzvf KLU-1.1.4.tar.gz  > /dev/null 2>&1
tar -xzvf UFconfig-3.7.1.tar.gz  > /dev/null 2>&1

# the original UFconfig.mk doesn't quite work for me
# use the hacked version
# 

echo " "
echo "Updating the config file ......"
echo " "
\mv -f UFconfig/UFconfig.mk UFconfig/UFconfig.mk.org
\cp -f UFconfig.mk.linux UFconfig/UFconfig.mk

echo "... Done!"
