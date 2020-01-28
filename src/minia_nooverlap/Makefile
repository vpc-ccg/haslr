all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
# build: clean-exe nooverlap clean
build: minia_nooverlap

LDFLAGS=
LIBS=-lz
CXXFLAGS=

objects=nooverlap.o

minia_nooverlap: nooverlap.o
	$(CXX) nooverlap.o -o $@ ${LDFLAGS} ${LIBS}

checkLinks: checkLinks.o
	$(CXX) checkLinks.o -o $@ ${LDFLAGS} ${LIBS}

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	@rm -f $(objects)
	
clean-exe:
	@rm -f nooverlap checkLinks

clean-all: clean-exe clean

DEBUG_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CXXFLAGS = $(CXXFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)
