#!/bin/bash  

rm ./fem_icpc 2> /dev/null
rm ./icpc_result.log 2> /dev/null
make -C ~/fem/ -f ~/fem/build/makefile.icpc clean
make -C ~/fem/ -f ~/fem/build/makefile.icpc build
cp ~/fem/fem_icpc ./fem_icpc

./fem_icpc --gtest_filter=cpu.china2015_test | tee -a icpc_result.log
zip -r $(date -d "today" +"%Y%m%d%H%M")"_icpc_log" icpc_result.log fem_icpc
rm ./fem_icpc
rm icpc_result.log