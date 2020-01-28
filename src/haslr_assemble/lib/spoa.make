all: clean-lib spoa/libspoa.a clean

INCS=-Ispoa/include
LIBS=
CXXFLAGS=-Wall -Wextra -pedantic -march=native -O3 -DNDEBUG -std=c++11 $(INCS)

objects=spoa/src/graph.o spoa/src/alignment_engine.o spoa/src/simd_alignment_engine.o spoa/src/sisd_alignment_engine.o 
# compile: $(objects)

spoa/libspoa.a: $(objects)
	$(AR) -qc $@ $(objects)
#	ranlib $@

# %.o: %.cpp
# 	$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

clean:
	@rm -f $(objects)
	
clean-lib:
	@rm -f spoa/libspoa.a
