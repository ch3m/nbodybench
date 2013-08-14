# Name of the file(s) we're compiling

TARGET = nbody

# location of the Python header files

PYTHON = /usr/include/python2.6

# location of the Boost Python include files and library
#BOOST = /home/ch3m/Tools/boost/boost-install
BOOST_INC = $(BOOST)/include
BOOST_LIB = $(BOOST)/boost-install/lib

$(TARGET).so: $(TARGET).o
	icpc -openmp -O3 -shared -lmkl_rt -L$(BOOST_LIB) -lboost_python  \
	-L/usr/lib64/python2.6/config -lpython2.6 \
  $(TARGET).o -o NbodyBench.so

$(TARGET).o: $(TARGET).cc
	icpc -openmp -O3 -vec-report2 -fPIC -I$(PYTHON) -I$(BOOST_INC) -c $(TARGET).cc

clean:
	rm -f *.o *.dll *.a *~ core

