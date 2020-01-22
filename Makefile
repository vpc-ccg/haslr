MINIMAP2_V   := 2.17
MINIA_V      := 3.2.1

SRCDIR       := src
THRID_PARTY  := vendor
BINDIR       := bin

all: deps haslr nooverlap

haslr: nooverlap haslr_assemble

haslr_assemble:
	$(MAKE) -C $(SRCDIR)/haslr_assemble
	@cp $(SRCDIR)/haslr_assemble/haslr_assemble $(BINDIR)/haslr_assemble

nooverlap:
	$(MAKE) -C $(SRCDIR)/minia_nooverlap
	@cp $(SRCDIR)/minia_nooverlap/nooverlap $(BINDIR)/nooverlap

deps:
	@if [ ! -f $(THRID_PARTY)/minimap2-${MINIMAP2_V}_x64-linux/minimap2 ]; then \
		mkdir -p $(THRID_PARTY)/minimap2-${MINIMAP2_V}_x64-linux; \
		cd $(THRID_PARTY); wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_V}/minimap2-${MINIMAP2_V}_x64-linux.tar.bz2; \
		tar -xvf minimap2-${MINIMAP2_V}_x64-linux.tar.bz2; \
	fi
	@if [ ! -f $(BINDIR)/minimap2 ]; then \
		cp $(THRID_PARTY)/minimap2-${MINIMAP2_V}_x64-linux/minimap2 $(BINDIR)/minimap2; \
	fi
	@if [ ! -f $(THRID_PARTY)/minia-v${MINIA_V}-bin-Linux/bin/minia ]; then \
		mkdir -p $(THRID_PARTY)/minia-v${MINIA_V}-bin-Linux; \
		cd $(THRID_PARTY); wget https://github.com/GATB/minia/releases/download/v${MINIA_V}/minia-v${MINIA_V}-bin-Linux.tar.gz; \
		tar -zxvf minia-v${MINIA_V}-bin-Linux.tar.gz; \
	fi
	@if [ ! -f $(BINDIR)/minia ]; then \
		cp $(THRID_PARTY)/minia-v${MINIA_V}-bin-Linux/bin/minia $(BINDIR)/minia; \
	fi
	@if [ ! -f $(THRID_PARTY)/fastutils-master/fastutils ]; then \
		mkdir -p $(THRID_PARTY)/fastutils-master; \
		cd $(THRID_PARTY); wget https://github.com/haghshenas/fastutils/archive/master.tar.gz -O fastutils-master.tar.gz; \
		tar -zxvf fastutils-master.tar.gz; \
		$(MAKE) -C fastutils-master; \
	fi
	@if [ ! -f $(BINDIR)/fastutils ]; then \
		cp $(THRID_PARTY)/fastutils-master/fastutils $(BINDIR)/fastutils; \
	fi

clean_deps:
	@rm -rf $(BINDIR)/minimap2 $(BINDIR)/minia $(BINDIR)/fastutils $(THRID_PARTY)

# purge-all: clean-exe clean purge-deps
# 	@rm -rf $(LIBDIR)/edlib $(LIBDIR)/ksw2 $(LIBDIR)/minimap2 $(LIBDIR)/seq-align $(LIBDIR)/spoa 

# clean:
# 	@rm -f $(OBJS)

# clean-exe:
# 	@rm -f $(PROG)

clean:
	$(MAKE) -C $(SRCDIR)/minia_nooverlap clean-all
	@rm $(BINDIR)/nooverlap
	$(MAKE) -C $(SRCDIR)/haslr_assemble clean-all
	@rm $(BINDIR)/haslr_assemble

clean_all: clean_deps clean

# clean-all: clean-exe clean
# 	@if [ -d $(LIBDIR)/seq-align ]; then \
# 		$(MAKE) -C $(LIBDIR)/seq-align clean; \
# 		$(MAKE) -C $(LIBDIR)/seq-align/libs/string_buffer clean; \
# 	fi
# 	@if [ -d $(LIBDIR)/ksw2 ]; then \
# 		$(MAKE) -C $(LIBDIR)/ksw2 clean; \
# 	fi
# 	@if [ -d $(LIBDIR)/minimap2 ]; then \
# 		$(MAKE) -C $(LIBDIR)/minimap2 clean; \
# 	fi
# 	@$(MAKE) -C $(LIBDIR) -f spoa.make clean clean-lib
