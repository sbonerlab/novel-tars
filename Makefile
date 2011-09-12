# required input from environment:
# $(BIOINFOCONFDIR) and $(BICOSN)

# ---------------- setup symbols needed in the rest of the file --------
CC = gcc
CFLAGS = -m64 -Wall -Wno-parentheses -Wno-sign-compare -Wno-unknown-pragmas -g
CFLAGSO = -m64 -Wall -Wno-parentheses -Wno-sign-compare -Wno-unknown-pragmas -Wno-uninitialized -g -O2
CCEXPAND = $(CC) -x c -E
BIOSINC = /home/as898/include
BIOSLNK = /home/as898/lib
BIOSLIB = $(BIOSLNK)/libbios.a

##include $(BIOINFOCONFDIR)/biosdefs.make
.SUFFIXES:
SHELL = /bin/bash


# ----------------------- entry points --------------

PROGRAMS=readBgr sequenceFilter annotateTARs
#convert2interval

MODULES=

all: allprogs 

allprogs: $(MODULES) $(PROGRAMS) 

clean: 
	/bin/rm -f $(PROGRAMS) $(MODULES)


sequenceFilter: sequenceFilter.c $(BIOSLIB)
	-@/bin/rm -f sequenceFilter
	$(CC) $(CFLAGS) -I$(BIOSINC) sequenceFilter.c -o sequenceFilter -L$(BIOSLNK) -lbios

convert2interval: convert2interval.c $(BIOSLIB)
	-@/bin/rm -f convert2interval
	$(CC) $(CFLAGS) -I/$(BIOSINC) convert2interval.c -o convert2interval -L$(BIOSLNK) -lm

readBgr: readBgr.c
	-@/bin/rm -f readBgr
	$(CC) $(CFLAGS) -I$(BIOSINC) readBgr.c -o readBgr -L$(BIOSLNK) -lbios -ltabix -lz

annotateTARs: annotateTARs.c
	-@/bin/rm -f annotateTARs
	$(CC) $(CFLAGS) -I$(BIOSINC) annotateTARs.c -o annotateTARs -L$(BIOSLNK) -lbios -lm
