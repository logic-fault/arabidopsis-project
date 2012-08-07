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

TSS_SOURCES=island_overlap_tss.c CpGIOverlap_stream/CpGIOverlap_stream.c 
SCORE_SOURCES=island_score.c CpGI_score_stream/CpGI_score_stream.c
NUC_SOURCES=nuc_score.c island_nuc_score_stream/island_nuc_score_stream.c
EXPRESSION_SOURCES=gene_expression_score.c gene_expression_score_stream/gene_expression_score_stream.c
TSS_OBJECTS=$(TSS_SOURCES:.c=.o)
SCORE_OBJECTS=$(SCORE_SOURCES:.c=.o)
EXPRESSION_OBJECTS=$(EXPRESSION_SOURCES:.c=.o)
NUC_OBJECTS=$(NUC_SOURCES:.c=.o)

all: $(TSS_SOURCES) $(SCORE_SOURCES) island_overlap_tss island_score expression_score nuc_score

island_overlap_tss: $(TSS_OBJECTS)
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $(TSS_OBJECTS) -lm -lgenometools -lcairo -o $@

island_score: $(SCORE_OBJECTS)
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $(SCORE_OBJECTS) -lm -lgenometools -lcairo -o $@

expression_score: $(EXPRESSION_OBJECTS)
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $(EXPRESSION_OBJECTS) -lm -lgenometools -lcairo -o $@

nuc_score: $(NUC_OBJECTS)
	$(LD) $(LDFLAGS) $(GT_LDFLAGS) $(NUC_OBJECTS) -lm -lgenometools -lcairo -o $@

# generic compilation rule which creates dependency file on the fly
.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) $(GT_CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read dependencies
-include $(wildcard *.d)

.PHONY: clean
clean:
	rm -f *.[od] island_overlap_tss island_score expression_score
