#!/bin/sh
if [ ! $1 ] ; then
    echo "Please specify the code you want to compile by typing :"
    echo "./Make <Your-Code.C>"
    exit 1
fi
echo "================================================================"
#echo "====> Producing eventdict.cc and eventdict.h"
#rootcint -f eventdict.cc -c interface/myevent.h interface/LinkDef.h 
#echo "====> Compiling $1 linked with eventdict.cc"
filename=`echo $1 | awk -F"." '{print $1}'`
exefilename=${filename}.exe
rm -f $exefilename
rm *.o
#g++ -c HTT-utilities/LepEffInterface/src/ScaleFactor.cc -I./ `root-config --cflags --glibs`
#g++ -c -I../bTagReweighting/ ../bTagReweighting/bTagSF.cc `root-config --cflags --glibs`
#g++ $1 -o $exefilename `root-config --cflags --glibs` -I../bTagReweighting/ -I./ HTT-utilities/LepEffInterface/src/ScaleFactor.cc ../bTagReweighting/bTagSF.cc ../bTagReweighting/BTagCalibrationStandalone.cc #../bTagReweighting/bTagSF.o ScaleFactor.o 
g++ $1 -o $exefilename `root-config --cflags --glibs` -I./ HTT-utilities/LepEffInterface/src/ScaleFactor.cc #../bTagReweighting/bTagSF.o ScaleFactor.o 
echo ""
if [ -e $exefilename ]; then 
    echo "====> Created exe file : "
    ls -lrt $exefilename
    echo "====> Done."
else
    echo "====> Did not create the exe file!"
fi
echo "================================================================"
	
