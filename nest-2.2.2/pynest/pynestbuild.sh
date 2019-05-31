#!/bin/sh

top_srcdir="/home/hp120306/k00825/Build/nest.git/"
flags1="-Xg -pthread -fPIC -fno-strict-aliasing -OPT:Olimit=0 -g -O2 -DNDEBUG -g -O3 -Wall -fPIC -DHAVE_NUMPY "
flags2="-O3 -Xg -fPIC -Kopenmp -fpermissive"
linkflags1="-Xg -pthread -fPIC -shared"
linkflags2="-lnest -lrandom -lsli -lfj90i -lfj90f -lgsl -lgslcblas -lm -lpthread -lpython2.7"
linkflags3="-Kopenmp"

mpiFCC ${flags1} -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/librandom -I/home/hp120306/k00825/Build/nest.git/sli -I/home/hp120306/k00825/Build/nest.git/nestkernel -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/nest -I/data/hp120306/k00825/nest.install/include -I${top_srcdir}/libltdl -I/opt/local/Python-2.7.3/lib/python2.7/site-packages/numpy/core/include -I/opt/local/Python-2.7.3/include/python2.7 -c /home/hp120306/k00825/Build/nest.git/pynest/pynestkernel.cpp -o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/pynestkernel.o ${flags2}

mpiFCC ${flags1} -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/librandom -I/home/hp120306/k00825/Build/nest.git/sli -I/home/hp120306/k00825/Build/nest.git/nestkernel -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/nest -I/data/hp120306/k00825/nest.install/include -I${top_srcdir}/libltdl -I/opt/local/Python-2.7.3/lib/python2.7/site-packages/numpy/core/include -I/opt/local/Python-2.7.3/include/python2.7 -c /home/hp120306/k00825/Build/nest.git/pynest/pydatum.cpp -o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/pydatum.o ${flags2} 

mpiFCC ${flags1} -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/librandom -I/home/hp120306/k00825/Build/nest.git/sli -I/home/hp120306/k00825/Build/nest.git/nestkernel -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/nest -I/data/hp120306/k00825/nest.install/include -I${top_srcdir}/libltdl -I/opt/local/Python-2.7.3/lib/python2.7/site-packages/numpy/core/include -I/opt/local/Python-2.7.3/include/python2.7 -c /home/hp120306/k00825/Build/nest.git/pynest/datumtopythonconverter.cpp -o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/datumtopythonconverter.o ${flags2} 

mpiFCC ${flags1} -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/librandom -I/home/hp120306/k00825/Build/nest.git/sli -I/home/hp120306/k00825/Build/nest.git/nestkernel -I/home/hp120306/k00825/Build/nest.git/libnestutil -I/home/hp120306/k00825/Build/nest.git/nest -I/data/hp120306/k00825/nest.install/include -I${top_srcdir}/libltdl -I/opt/local/Python-2.7.3/lib/python2.7/site-packages/numpy/core/include -I/opt/local/Python-2.7.3/include/python2.7 -c /home/hp120306/k00825/Build/nest.git/pynest/pynestpycsa.cpp -o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/pynestpycsa.o ${flags2} 

echo "*** linking ... ***"

mpiFCC ${linkflags1} /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/pynestkernel.o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/pydatum.o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/datumtopythonconverter.o /home/hp120306/k00825/Build/nest.git/pynest/build/temp.linux-s64fx-2.7/home/hp120306/k00825/Build/nest.git/pynest/pynestpycsa.o /home/hp120306/k00825/Build/nest.git/models/.libs/libmodelsmodule.a /home/hp120306/k00825/Build/nest.git/precise/.libs/libprecisemodule.a /home/hp120306/k00825/Build/nest.git/topology/.libs/libtopologymodule.a /home/hp120306/k00825/Build/nest.git/conngen/.libs/libconngenmodule.a -L/home/hp120306/k00825/Build/nest.git/nestkernel/.libs -L/home/hp120306/k00825/Build/nest.git/librandom/.libs -L/home/hp120306/k00825/Build/nest.git/sli/.libs -L/data/hp120306/k00825/nest.install/lib -L/opt/local/Python-2.7.3/lib ${linkflags2} -o /home/hp120306/k00825/Build/nest.git/pynest/build/lib.linux-s64fx-2.7/nest/pynestkernel.so ${linkflags3}

echo "*** copying ... ***"

cp /home/hp120306/k00825/Build/nest.git/pynest/build/lib.linux-s64fx-2.7/nest/pynestkernel.so /data/hp120306/k00825/nest.install/lib/python2.7/site-packages/nest/

