#! /bin/bash

echo "Fetch cached modules ......"
echo " "
if [ ! -f src/src.tar.gz ]; then
## original location
##curl ftp://ftp.nist.gov/pub/cmlib/src.tar.Z -o src/src.tar.gz
cp ../Repositories/CMLIB/cmlib-3.0.tar.gz src/src.tar.gz
fi

echo " "
echo "Expanding the modules ......"
echo " "

cd src
tar -xzvf src.tar.gz > /dev/null 2>&1 
cd -

echo "... Done!"

#end
