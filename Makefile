# Number of cores to use when invoking parallelism
#ifndef CORES
    CORES := 1 # temporarily set to 1 since I broke threading. :-(
#endif
ifndef PAUSE
    PAUSE := 100
endif
# Uncomment either of these to remove them (removing 7 implies removing 8)
EIGHT := 8
SEVEN := 7
ifdef NO8
    EIGHT := 
endif
ifdef NO7
    SEVEN :=
    EIGHT := # can't have 8 without 7
endif

ifdef DEBUG
    ifdef PROFILE
	SPEED=-O0 -ggdb -pg
	LIB_OPT=-g-pg
    else
	SPEED=-O0 -ggdb
	LIB_OPT=-g
    endif
else
    ifdef PROFILE
	SPEED=-O3 -pg
	LIB_OPT=-pg
    else
	SPEED=-O3
	LIB_OPT=
    endif
endif

# Waywe needs gcc-6 on MacOS:
GCC_VER=$(shell echo $(ARCH) $(HOME) | awk '/Darwin/&&/Users.wayne/{V="-6"}END{if(V)print V;else{printf "using default gcc: " > "/dev/null"; exit 1}}')
GCC=gcc$(GCC_VER)
CXX=g++

# Some architectures, eg CYGWIN 32-bit and MacOS("Darwin") need an 80MB stack.
export LIBWAYNE_HOME=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/libwayne
ARCH=$(shell uname -a | awk '{if(/CYGWIN/){V="CYGWIN"}else if(/Darwin/){V="Darwin"}else if(/Linux/){V="Linux"}}END{if(V){print V;exit}else{print "unknown OS" > "/dev/stderr"; exit 1}}')

# Darwin needs gcc-6 ever since a commit on 22 May 2022:
#GCC= $(shell $(CC) -v 2>&1 | awk '/gcc/{++gcc}{V=$$3}END{if(gcc && (V ~ /[0-9]\.[0-9]\.[0-9]*/))print "$(ARCH).gcc"V; else exit 1}')
STACKSIZE=$(shell ($(GCC) -v 2>&1; uname -a) | awk '/CYGWIN/{print "-Wl,--stack,83886080"}/gcc-/{actualGCC=1}/Darwin/&&actualGCC{print "-Wl,-stack_size -Wl,0x5000000"}')
CC=$(GCC) $(SPEED) -Wall -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wshadow $(PG) $(STACKSIZE)
LIBWAYNE_COMP=-I $(LIBWAYNE_HOME)/include $(STACKSIZE) $(SPEED)
LIBWAYNE_LINK=-L $(LIBWAYNE_HOME) -lwayne$(LIB_OPT) -lm $(STACKSIZE) $(SPEED)
LIBWAYNE_BOTH=$(LIBWAYNE_COMP) $(LIBWAYNE_LINK)

# Name of BLANT source directory
SRCDIR = src
# Put all c files in SRCDIR below.
BLANT_SRCS = blant.c \
			 blant-window.c \
			 blant-output.c \
			 blant-utils.c \
			 blant-sampling.c \
			 blant-predict.o \
			 blant-synth-graph.c \
			 importance.c \
			 odv.c

OBJDIR = _objs
OBJS = $(addprefix $(OBJDIR)/, $(BLANT_SRCS:.c=.o))

### Generated File Lists ###
K := 3 4 5 6 $(SEVEN) $(EIGHT)
canon_txt := canon_maps/canon_map canon_maps/canon_list canon_maps/canon-ordinal-to-signature canon_maps/orbit_map canon_maps/alpha_list_nbe canon_maps/alpha_list_mcmc
canon_bin := canon_maps/canon_map canon_maps/perm_map
canon_all := $(foreach k, $(K), $(addsuffix $(k).txt, $(canon_txt)) $(addsuffix $(k).bin, $(canon_bin)))
subcanon_txts := $(if $(EIGHT),canon_maps/subcanon_map8-7.txt) $(if $(SEVEN),canon_maps/subcanon_map7-6.txt) canon_maps/subcanon_map6-5.txt canon_maps/subcanon_map5-4.txt canon_maps/subcanon_map4-3.txt

#canon_map_bins := $(foreach k,$(K), canon_maps/canon_map$(k).bin)
#perm_map_bins := $(foreach k,$(K), canon_maps/perm_map$(k).bin)
#canon_map_txts := $(foreach k,$(K), canon_maps/canon_map$(k).txt)
#canon_list_txts := $(foreach k,$(K), canon_maps/canon_list$(k).txt)
#canon_ordinal_to_signature_txts := $(foreach k,$(K), canon_maps/canon-ordinal-to-signature$(k).txt)
#orbit_map_txts := $(foreach k,$(K), canon_maps/orbit_map$(k).txt)
#canon_map_files := $(canon_map_bins) $(perm_map_bins) $(canon_map_txts) $(canon_list_txts) $(canon_ordinal_to_signature_txts) $(orbit_map_txts)

# ehd takes up too much space and isn't used anywhere yet
#ehd_txts := $(foreach k,$(K), canon_maps/EdgeHammingDistance$(k).txt)
#alpha_nbe_txts := $(foreach k, $(K), canon_maps/alpha_list_nbe$(k).txt)
#alpha_mcmc_txts := $(foreach k, $(K), canon_maps/alpha_list_mcmc$(k).txt)
magic_table_txts := $(foreach k,$(K), orca_jesse_blant_table/UpperToLower$(k).txt)

#base: show-gcc-ver ./.notpristine libwayne $(alpha_nbe_txts) $(alpha_mcmc_txts) magic_table $(canon_map_files) blant
base: show-gcc-ver ./.notpristine libwayne $(canon_all) magic_table blant

show-gcc-ver:
	$(GCC) -v

./.notpristine:
	@echo '****************************************'
	@echo "If you haven't already, you should read the README at"
	@echo "	https://github.com/waynebhayes/BLANT#readme"
	@echo "BLANT can sample graphlets of up to k=8 nodes. The lookup table for k=8 can take"
	@echo "up to an hour to generate, but is needed if you want BLANT-seed to work, and so"
	@echo "we make it by default; set NO8=1 to turn it off."
	@echo "The fastest way to get started is to skip k=8 graphlets:"
	@echo "    PAUSE=0 NO8=1 make base"
	@echo "which will make everything needed to get started sampling up to k=7 graphlets".
	@echo "You will only see this message once on a 'pristine' repo. Pausing $(PAUSE) seconds..."
	@echo '****************************************'
	sleep $(PAUSE)
	@touch .notpristine

most: base Draw subcanon_maps

all: most $(ehd_txts)

gcc-version:
	@if ls $(SRCDIR)/*.gz | fgrep -q "$(GCC)"; then :; else echo "WARNING: gcc version not supported; prediction disabled" >&2; fi

canon_maps: base $(canon_all) subcanon_maps

.PHONY: all most pristine clean_canon_maps

### Executables ###

fast-canon-map: libwayne $(SRCDIR)/fast-canon-map.c | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CC) '-std=c99' -O3 -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/fast-canon-map.c $(LIBWAYNE_BOTH)

slow-canon-maps: libwayne $(SRCDIR)/slow-canon-maps.c | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CC) -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/slow-canon-maps.c $(LIBWAYNE_BOTH)

make-orbit-maps: libwayne $(SRCDIR)/make-orbit-maps.c | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CC) -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/make-orbit-maps.c $(LIBWAYNE_BOTH)

blant: gcc-version libwayne $(OBJS) $(OBJDIR)/convert.o $(OBJDIR)/libblant.o | $(LIBWAYNE_HOME)/MT19937/mt19937.o
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(OBJS) $(OBJDIR)/convert.o $(LIBWAYNE_HOME)/MT19937/mt19937.o $(LIBWAYNE_LINK)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(LIBWAYNE_COMP)

synthetic: libwayne $(SRCDIR)/synthetic.c $(SRCDIR)/syntheticDS.h $(SRCDIR)/syntheticDS.c | $(OBJDIR)/libblant.o
	$(CC) -c $(SRCDIR)/syntheticDS.c $(SRCDIR)/synthetic.c $(LIBWAYNE_COMP)
	$(CXX) -o $@ syntheticDS.o $(OBJDIR)/libblant.o synthetic.o $(LIBWAYNE_LINK)

CC: libwayne $(SRCDIR)/CC.c $(OBJDIR)/convert.o | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(SRCDIR)/CC.c $(OBJDIR)/convert.o $(LIBWAYNE_LINK)

makeEHD: $(OBJDIR)/makeEHD.o
	$(CXX) -o $@ $(OBJDIR)/libblant.o $(OBJDIR)/makeEHD.o $(LIBWAYNE_LINK)

compute-alphas-NBE: libwayne $(SRCDIR)/compute-alphas-NBE.c | $(OBJDIR)/libblant.o
	$(CC) -Wall -O3 -o $@ $(SRCDIR)/compute-alphas-NBE.c $(OBJDIR)/libblant.o $(LIBWAYNE_BOTH)

compute-alphas-MCMC: libwayne $(SRCDIR)/compute-alphas-MCMC.c | $(OBJDIR)/libblant.o
	$(CC) -Wall -O3 -o $@ $(SRCDIR)/compute-alphas-MCMC.c $(OBJDIR)/libblant.o $(LIBWAYNE_BOTH)

Draw: Draw/graphette2dot

Draw/graphette2dot: libwayne Draw/DrawGraphette.cpp Draw/Graphette.cpp Draw/Graphette.h Draw/graphette2dotutils.cpp Draw/graphette2dotutils.h  | $(SRCDIR)/blant.h $(OBJDIR)/libblant.o
	$(CXX) Draw/DrawGraphette.cpp Draw/graphette2dotutils.cpp Draw/Graphette.cpp $(OBJDIR)/libblant.o -o $@ -std=c++11 $(LIBWAYNE_BOTH)

make-subcanon-maps: libwayne $(SRCDIR)/make-subcanon-maps.c | $(OBJDIR)/libblant.o
	$(CC) -Wall -o $@ $(SRCDIR)/make-subcanon-maps.c $(OBJDIR)/libblant.o $(LIBWAYNE_BOTH)

make-orca-jesse-blant-table: libwayne $(SRCDIR)/magictable.cpp | $(OBJDIR)/libblant.o
	$(CXX) -Wall -o $@ $(SRCDIR)/magictable.cpp $(OBJDIR)/libblant.o -std=c++11 $(LIBWAYNE_BOTH)

$(OBJDIR)/blant-predict.o:
	if [ -f $(SRCDIR)/EdgePredict/blant-predict.c ]; then (cd $(SRCDIR)/EdgePredict && ../../libwayne/bin/wgcc -c -o ../blant-predict.o blant-predict.c; cd ..; cp -p blant-predict.o ../_objs); elif [ -f $(SRCDIR)/blant-predict.o ]; then cat $(SRCDIR)/blant-predict.o; elif [ "$(ARCH)" = Darwin ]; then gunzip < $(SRCDIR)/blant-predict.Darwin.o.gz; elif [ -f $(SRCDIR)/blant-predict.$(GCC).o.gz ]; then gunzip < $(SRCDIR)/blant-predict.$(GCC).o.gz; else $(CC) -c -o $(SRCDIR)/blant-predict.o $(SRCDIR)/blant-predict-stub.c $(LIBWAYNE_BOTH); cat $(SRCDIR)/blant-predict.o; fi > $@

### Object Files/Prereqs ###

$(OBJDIR)/convert.o: $(SRCDIR)/convert.cpp
	@mkdir -p $(dir $@)
	$(CXX) -c $(SRCDIR)/convert.cpp -o $@ -std=c++11

$(LIBWAYNE_HOME)/MT19937/mt19937.o: libwayne
	cd $(LIBWAYNE_HOME)/MT19937 && $(MAKE)

$(OBJDIR)/libblant.o: libwayne $(SRCDIR)/libblant.c
	@mkdir -p $(dir $@)
	$(CC) -c $(SRCDIR)/libblant.c -o $@ $(LIBWAYNE_COMP)

$(OBJDIR)/makeEHD.o: libwayne $(SRCDIR)/makeEHD.c | $(OBJDIR)/libblant.o
	@mkdir -p $(dir $@)
	$(CC) -c $(SRCDIR)/makeEHD.c -o $@ $(LIBWAYNE_COMP)

libwayne: $(LIBWAYNE_HOME)/Makefile $(LIBWAYNE_HOME)/made

$(LIBWAYNE_HOME)/made:
	cd $(LIBWAYNE_HOME) && $(MAKE) all

### Generated File Recipes ###

canon_maps/canon_map%.bin canon_maps/perm_map%.bin canon_maps/orbit_map%.txt canon_maps/alpha_list_mcmc%.txt: libwayne $(SRCDIR)/create-bin-data.c | $(OBJDIR)/libblant.o $(SRCDIR)/blant.h canon_maps/canon_list%.txt canon_maps/canon_map%.txt make-orbit-maps compute-alphas-MCMC
	$(CC) '-std=c99' "-Dkk=$*" "-DkString=\"$*\"" -o create-bin-data$* $(SRCDIR)/libblant.c $(SRCDIR)/create-bin-data.c $(LIBWAYNE_BOTH)
	[ -f canon_maps/canon_map$*.bin -a -f canon_maps/perm_map$*.bin ] || ./create-bin-data$*
	./make-orbit-maps $* > canon_maps/orbit_map$*.txt
	@if [ -f canon_maps.correct/alpha_list_mcmc$*.txt ]; then echo "computing MCMC alphas for k=$* takes days, so just copy it"; cp canon_maps.correct/alpha_list_mcmc$*.txt canon_maps/ && touch $@; else ./compute-alphas-MCMC $* > canon_maps/alpha_list_mcmc$*.txt; fi

canon_maps/canon_map%.txt canon_maps/canon_list%.txt canon_maps/canon-ordinal-to-signature%.txt: fast-canon-map
	mkdir -p canon_maps
	# It's cheap to make all but k=8 canon maps, so make all but skip 8 if it already exists
	[ $* -eq 8 -a '(' -f canon_maps/canon_map$*.txt -o -f canon_maps/canon_map$*.txt.gz ')' ] || ./fast-canon-map $* | tee canon_maps/canon_map$*.txt | awk -F '	' 'BEGIN{n=0}!seen[$$1]{seen[$$1]=$$0;map[n++]=$$1}END{print n;for(i=0;i<n;i++)print seen[map[i]]}' | cut -f1,3- | tee canon_maps/canon_list$*.txt | awk 'NR>1{print NR-2, $$1}' > canon_maps/canon-ordinal-to-signature$*.txt
	if [ $* -eq 8 -a -f canon_maps/canon_map$*.txt -a ! -f canon_maps/canon_map$*.txt.gz ]; then gzip canon_maps/canon_map$*.txt & fi

canon_maps/EdgeHammingDistance%.txt: makeEHD | canon_maps/canon_list%.txt canon_maps/canon_map%.bin
	@if [ ! -f canon_maps.correct/EdgeHammingDistance$*.txt.xz ]; then ./makeEHD $* > $@; cmp canon_maps.correct/EdgeHammingDistance$*.txt $@; else echo "EdgeHammingDistance8.txt takes weeks to generate; uncompressing instead"; unxz < canon_maps.correct/EdgeHammingDistance$*.txt.xz > $@ && touch $@; fi
	#(cd canon_maps.correct && ls EdgeHammingDistance$*.txt*) | awk '{printf "cmp canon_maps.correct/%s canon_maps/%s\n",$$1,$$1}' | sh

canon_maps/alpha_list_nbe%.txt: compute-alphas-NBE canon_maps/canon_list%.txt
	./compute-alphas-NBE $* > $@

.INTERMEDIATE: .created-subcanon-maps
subcanon_maps: $(subcanon_txts) ;
$(subcanon_txts): .created-subcanon-maps
.created-subcanon-maps: make-subcanon-maps | $(canon_all) #$(canon_list_txts) $(canon_map_bins)
	# only do it for k > 3 since it's 4-3, 5-4, etc.
	for k in $(K); do if [ $$k -gt 3 ]; then ./make-subcanon-maps $$k > canon_maps/subcanon_map$$k-$$(($$k-1)).txt; fi; done

magic_table: $(magic_table_txts) ;
$(magic_table_txts): make-orca-jesse-blant-table | $(canon_all) #$(canon_list_txts) $(canon_map_bins)
	./make-orca-jesse-blant-table $(if $(EIGHT),8,$(if $(SEVEN),7,6))

### Cleaning ###

clean:
	@/bin/rm -f *.[oa] blant canon-sift fast-canon-map make-orbit-maps compute-alphas-MCMC compute-alphas-NBE makeEHD make-orca-jesse-blant-table Draw/graphette2dot blant-sanity make-subcanon-maps create-bin-data*
	@/bin/rm -rf $(OBJDIR)/*

pristine: clean clean_canon_maps
	@cd $(LIBWAYNE_HOME); $(MAKE) clean
	@/bin/rm -f canon_maps/* .notpristine .firsttime # .firsttime is the old name but remove it anyway

clean_canon_maps:
	@/bin/rm -f canon_maps/*[3-7].* # don't remove 8 since it takes too long to create
	@/bin/rm -f orca_jesse_blant_table/UpperToLower*.txt
