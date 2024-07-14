CXXFLAGS = -O3 -Wall -std=c++11

ALLPROGRAM = ../bin/bf

all: $(ALLPROGRAM)

clean:
	rm -f $(ALLPROGRAM) *.o* */*.o*

Obj = main.o PairHMM.o Sequence.o 

../bin/bf: $(Obj)
	$(CXX) $(CXXFLAGS) -o $@ $(Obj)

PairHMM.o: PairHMM.cc PairHMM.hh Sequence.hh PairHMMPar.hh Utils.hh
Sequence.o: Sequence.cc Sequence.hh
main.o: main.cc Sequence.hh PairHMM.hh PairHMMPar.hh Utils.hh
