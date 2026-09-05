CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O3 -march=native -flto -DNDEBUG -Wall -Wextra
LDFLAGS  ?= -pthread -flto

all: stab14 stab14_ext stab14_2stage

stab14: stab14.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

stab14_ext: stab14_ext.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

stab14_2stage: stab14_2stage.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run-2stage: stab14_2stage
	./stab14_2stage --mode hybrid

clean:
	rm -f stab14 stab14_ext stab14_2stage *.exe *.obj
