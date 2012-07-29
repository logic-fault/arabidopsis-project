/* @file
 * @author brock a <brock.wright.anderson@gmail.cm>
 *
 *
 *
 *
 *
 *
 */

#include <genometools.h>
#include "CpGIOverlap_stream.h"


typedef struct
{
    unsigned long start;
    unsigned long end;
} island_t

struct CpGIOverlap_stream {
    const GtNodeStream parent_instance;
    GtNodeStream * in_stream;
    unsigned long index;
    gtGenomeNode * latest_CpGI_node;
    gtGenomeNode * latest_gene_node;
    island_t island;
    unsigned long latest_tss;
};

#define CpGIOverlap_stream_cast(GS) CpGIOverlap_stream_cast(CpGI_overlap_stream_class(), GS);


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


 
}

static int CpGIOverlap_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    return 0;
}
