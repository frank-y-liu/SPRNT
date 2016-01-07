#! /bin/bash

echo "Fetch modules from ftp.nist.gov ......"
echo " "
if [ ! -f src/src.tar.gz ]; then
curl ftp://ftp.nist.gov/pub/cmlib/src.tar.Z -o src/src.tar.gz
fi

echo " "
echo "Expanding the modules ......"
echo " "

cd src
tar -xzvf src.tar.gz > /dev/null 2>&1 
cd -

echo "... Done!"

#end
