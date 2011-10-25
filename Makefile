# required input from environment:
# $(BIOINFOCONFDIR) and $(BICOSN)

# ---------------- setup symbols needed in the rest of the file --------
CC = gcc
CFLAGS = -m64 -Wall -Wno-parentheses -Wno-sign-compare -Wno-unknown-pragmas -g
CFLAGSO = -m64 -Wall -Wno-parentheses -Wno-sign-compare -Wno-unknown-pragmas -Wno-uninitialized -g -O2
CCEXPAND = $(CC) -x c -E


##include $(BIOINFOCONFDIR)/biosdefs.make
.SUFFIXES:
SHELL = /bin/bash


# ----------------------- entry points --------------

PROGRAMS=readBgr sequenceFilter annotateTARs mergeTARs tarintron2interval
#convert2interval

MODULES=

all: allprogs 

allprogs: $(MODULES) $(PROGRAMS) 

clean: 
	/bin/rm -f $(PROGRAMS) $(MODULES)


sequenceFilter: sequenceFilter.c 
	-@/bin/rm -f sequenceFilter
	$(CC) $(CFLAGS) $(CPPFLAGS) sequenceFilter.c -o sequenceFilter $(LDFLAGS) -lbios

convert2interval: convert2interval.c 
	-@/bin/rm -f convert2interval
	$(CC) $(CFLAGS) $(CPPFLAGS) convert2interval.c -o convert2interval $(LDFLAGS) -lm

readBgr: readBgr.c
	-@/bin/rm -f readBgr
	$(CC) $(CFLAGS) $(CPPFLAGS) readBgr.c -o readBgr $(LDFLAGS) -lbios -lm 

mergeTARs: mergeTARs.c
	-@/bin/rm -f mergeTARs
	$(CC) $(CFLAGS) $(CPPFLAGS) mergeTARs.c -o mergeTARs $(LDFLAGS) -lbios -lm -lgsl -lgslcblas

tarintron2interval: tarintron2interval.c
	-@/bin/rm -f tarintron2interval
	$(CC) $(CFLAGS) $(CPPFLAGS) tarintron2interval.c -o tarintron2interval $(LDFLAGS) -lbios -lm

annotateTARs: annotateTARs.c
	-@/bin/rm -f annotateTARs
	$(CC) $(CFLAGS) $(CPPFLAGS) annotateTARs.c -o annotateTARs $(LDFLAGS) -lbios -lm

