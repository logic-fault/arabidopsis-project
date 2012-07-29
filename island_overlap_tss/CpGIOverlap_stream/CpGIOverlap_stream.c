/*
* @file
* @author Brock Anderson <brock.wright.anderson@gmail.com>
* @section LICENSE
* Released to public domain without restriction
*
* @section DESCRIPTION
* Find CpGI that overlap TSS, mark as such
*
*
*
* Scenarios:
*
* CpGI              <---------------------------->
* genes     
*  (+) contained               <-----|     
*  (-) contained                       |----->
*  (+) periphery                             <--------|
*  (-) periphery |-------->
*
*
*  When you read a CpGI:
*     latest_CpGINode  <- node_pointer
*     cpgi_start       <- node_pointer->start
*     cpgi_end         <- node_pointe->end
*
*     // see if (-) periphery case satisfied
*     if (in_range(latest_tss, cpgi_start, cpgi_end))
*         mark(latest_CpGINode)  
*         mark(latest_gene)  // tell which CpGINode also
*       
*
*  When read gene: 
*                  latest_tss        <- gene's tss
*                  latest_gene       <- gene reference
*
*                  // below case covers contained statements and
*                  // (+) periphery
*                  if (in_range(latest_tss,cpgi_start,cpgi_end)
*                      mark(latest_CpGINode)
*                      mark(latest_gene) // tell which CpGINode also
                   

*************************************************/
#include <genometools.h>
#include "CpGIOverlap_stream.h"


typedef struct
{
    unsigned long start;
    unsigned long end;
} island_t;

struct CpGIOverlap_stream {
    const GtNodeStream parent_instance;
    GtNodeStream * in_stream;
    unsigned long index;
    GtGenomeNode * latest_CpGI_node;
    GtGenomeNode * latest_gene_node;
    island_t island;
    unsigned long latest_tss;
};

#define CpGIOverlap_stream_cast(GS) gt_node_stream_cast(CpGIOverlap_stream_class(), GS);


static int CpGIOverlap_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    return 0;
}

static void CpGIOverlap_stream_free(GtNodeStream * ns)
{

    return;
}

const GtNodeStreamClass * CpGIOverlap_stream_class(void)
{
    static const GtNodeStreamClass * c = NULL;
    if (!c)
    {
        c = gt_node_stream_class_new( sizeof(CpGIOverlap_stream),
                                      CpGIOverlap_stream_free,
                                      CpGIOverlap_stream_next
                                    );
    }
    
    return c;
}

GtNodeStream * CpGIOverlap_stream_new(GtNodeStream * in_stream)
{
    GtNodeStream * ns = gt_node_stream_create(CpGIOverlap_stream_class(), 
                                              true); // must be sorted
    CpGIOverlap_stream * overlap_stream = CpGIOverlap_stream_cast(ns);
    gt_assert(in_stream);
    overlap_stream->in_stream = gt_node_stream_ref(in_stream);

    overlap_stream->index = 0;

    // we haven't found an island yet
    overlap_stream->island.start = 0;
    overlap_stream->island.end   = 0;
    overlap_stream->latest_tss   = 0;

    // haven't found a cpGI or gene yet
    overlap_stream->latest_CpGI_node = NULL;
    overlap_stream->latest_gene_node = NULL;

    return ns;
 
}
