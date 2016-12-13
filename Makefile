INCLUDEDIRSGSL = -I//gwork/ilie/silky-kraken/gsl-19/include
LIBDIRS = -L//gwork/ilie/silky-kraken/gsl-19/lib
LIBS= -lgsl -lgslcblas -lm

MAKE=make
CXX=g++ -Wall


hitec: hitec.obj work.obj divsufsort.obj
	$(CXX) $(INCLUDEDIRSGSL) $(LIBDIRS) -O3 hitec.obj work.obj divsufsort.obj  -o $@ $(LIBS) 


work.obj : work.cc work.hh
	$(CXX) $(INCLUDEDIRSGSL) -O3 -c work.cc  -o work.obj

divsufsort.obj : divsufsort.c divsufsort.h
	$(CXX) -O3 -fomit-frame-pointer -Wall  -c divsufsort.c -o divsufsort.obj


hitec.obj : hitec.cc work.hh divsufsort.h
	$(CXX) $(INCLUDEDIRSGSL) -O3 -c hitec.cc  -o hitec.obj    

clean:
	rm -f hitec.obj
	rm -f work.obj
	rm -f hitec
	rm -f divsufsort.obj
