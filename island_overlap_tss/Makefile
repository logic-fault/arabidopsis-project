# genometools project 2012
# Brock Anderson 2012


# set prefix for libgenometools
prefix ?= /usr/local

CC:=gcc
GT_CFLAGS:=-Wall -Wunused-parameter -g -Os -pipe \
           -I$(prefix)/include/genometools \
           -I/usr/include/cairo \
           -I/usr/local/include/cairo \
           -I/sw/include/cairo \
           -I/opt/local/include/cairo
LD:=$(CC)
GT_LDFLAGS:=-L$(prefix)/lib \
            -L/usr/X11R6/lib \
            -L/usr/local/lib \
            -L/opt/local/lib

SOURCES=island_overlap_tss.c CpGIOverlap_stream/CpGIOverlap_stream.c
OBJECTS=$(SOURCES:.c=.o)

all: $(SOURCES) island_overlap_tss

island_overlap_tss: $(OBJECTS)
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $(OBJECTS) -lm -lgenometools -lcairo -o $@


# generic compilation rule which creates dependency file on the fly
.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read dependencies
-include $(wildcard *.d)

.PHONY: clean
clean:
	rm -f *.[od] island_overlap_tss
