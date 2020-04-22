MINIMAP2_V   := 2.17
MINIA_V      := 3.2.1
FASTUTILS_V  := 0.3

SRCDIR       := src
THRID_PARTY  := vendor
BINDIR       := bin

all: deps haslr

haslr: nooverlap haslr_assemble

haslr_assemble:
	$(MAKE) -C $(SRCDIR)/haslr_assemble
	@cp $(SRCDIR)/haslr_assemble/haslr_assemble $(BINDIR)/haslr_assemble

nooverlap:
	$(MAKE) -C $(SRCDIR)/minia_nooverlap
	@cp $(SRCDIR)/minia_nooverlap/minia_nooverlap $(BINDIR)/minia_nooverlap

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
	@if [ ! -f $(THRID_PARTY)/fastutils-${FASTUTILS_V}/fastutils ]; then \
		mkdir -p $(THRID_PARTY)/fastutils-${FASTUTILS_V}; \
		cd $(THRID_PARTY); wget https://github.com/haghshenas/fastutils/archive/v${FASTUTILS_V}.tar.gz -O fastutils-${FASTUTILS_V}.tar.gz; \
		tar -zxvf fastutils-${FASTUTILS_V}.tar.gz; \
		$(MAKE) -C fastutils-${FASTUTILS_V}; \
	fi
	@if [ ! -f $(BINDIR)/fastutils ]; then \
		cp $(THRID_PARTY)/fastutils-${FASTUTILS_V}/fastutils $(BINDIR)/fastutils; \
	fi

clean-deps:
	@rm -rf $(BINDIR)/minimap2 $(BINDIR)/minia $(BINDIR)/fastutils $(THRID_PARTY)

clean:
	$(MAKE) -C $(SRCDIR)/minia_nooverlap clean-all
	@rm -f $(BINDIR)/nooverlap
	$(MAKE) -C $(SRCDIR)/haslr_assemble clean-all
	@rm -f $(BINDIR)/haslr_assemble

clean-all: clean-deps
	$(MAKE) -C $(SRCDIR)/minia_nooverlap clean-all
	@rm -f $(BINDIR)/nooverlap
	$(MAKE) -C $(SRCDIR)/haslr_assemble purge
	@rm -f $(BINDIR)/haslr_assemble
