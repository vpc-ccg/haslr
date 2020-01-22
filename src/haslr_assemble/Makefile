SPOA_V       := 1.1.3

CC           ?= gcc
CXX          ?= g++

SRCDIR       := src
LIBDIR       := lib
THRID_PARTY  := vendor
BINDIR       := bin
PROG         := haslr_assemble

INCS         := -I$(LIBDIR)/spoa/include/spoa
LIBS         := $(LIBDIR)/spoa/libspoa.a -lz -lpthread #-lm 
# passed for both C and C++
CPPFLAGS     := 
# passed for C
CFLAGS       := 
# passed for C++
CXXFLAGS     := -g -Wall -std=c++11 $(INCS) 

OBJS         := $(SRCDIR)/main.o \
                $(SRCDIR)/Common.o \
                $(SRCDIR)/Commandline.o \
                $(SRCDIR)/Compressed_sequence.o \
                $(SRCDIR)/Contig.o \
                $(SRCDIR)/Longread.o \
                $(SRCDIR)/Backbone_graph.o \
                $(SRCDIR)/Cleaning.o \
                $(SRCDIR)/Assemble.o \
                # $(SRCDIR)/Graph_repeat.o \
                # $(SRCDIR)/Align_LR2path.o \

all: OPTIMIZE_FLAGS debug
profile: PROFILE_FLAGS debug
debug: clean-exe $(PROG) clean
# debug: clean-exe $(PROG)

$(PROG): $(LIBDIR)/spoa/libspoa.a $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)

$(LIBDIR)/spoa/libspoa.a:
	@if [ ! -d $(LIBDIR)/spoa ]; then \
		mkdir -p $(LIBDIR); \
		git clone https://github.com/rvaser/spoa.git $(LIBDIR)/spoa; \
		cd $(LIBDIR)/spoa; git checkout -q $(SPOA_V); \
	fi
	@$(MAKE) -C $(LIBDIR) -f spoa.make

clean:
	@rm -f $(OBJS)

clean-exe:
	@rm -f $(PROG)

clean-all: clean-exe clean
	@$(MAKE) -C $(LIBDIR) -f spoa.make clean clean-lib

purge: clean-all
	@rm -rf $(LIBDIR)/spoa

OPTIMIZE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -O2)
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
	$(eval CXXFLAGS = $(CXXFLAGS) -pg)
	$(eval LIBS = $(LIBS) -pg)
