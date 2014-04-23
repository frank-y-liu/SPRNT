#! /bin/bash

echo "Fetch modules from lib.stat.cmu.edu ......"
echo " "
if [ ! -f dbsplin_archive ]; then
curl http://lib.stat.cmu.edu/cmlib/src/dbsplin/archive -o dbsplin_archive
fi
if [ ! -f dbsplin_depend ]; then
curl http://lib.stat.cmu.edu/cmlib/src/dbsplin/depend -o dbsplin_depend
fi
if [ ! -f dtensbs_archive ]; then
curl http://lib.stat.cmu.edu/cmlib/src/dtensbs/archive -o dtensbs_archive
fi

chmod 700 dbsplin_archive
chmod 700 dbsplin_depend
chmod 700 dtensbs_archive

echo " "
echo "Expanding the modules ......"
echo " "
./dbsplin_archive > /dev/null 2>&1
./dbsplin_depend  > /dev/null 2>&1
./dtensbs_archive > /dev/null 2>&1

echo "... Done!"


#end
