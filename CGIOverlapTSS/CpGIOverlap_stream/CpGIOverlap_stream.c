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

#define CpGIOverlap_stream_cast(GS)\
        CpGIOverlap_stream_cast(CpGI_overlap_stream_class(), GS);


static int CpGIOverlap_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
}
