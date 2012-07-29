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
#include "stdio.h"

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
    GtArray      * node_buffer;
    island_t island;
    unsigned long latest_tss;
    int emptying_buffer;
    int building_buffer;
};

typedef enum { FT_CPGI, FT_GENE, FT_OTHER } overlap_feature_type_t;

const char * feature_type_CpGI = "CpGI";
const char * feature_type_gene = "gene";


const GtNodeStreamClass * CpGIOverlap_stream_class(void);

#define CpGIOverlap_stream_cast(GS) gt_node_stream_cast(CpGIOverlap_stream_class(), GS);

overlap_feature_type_t determine_feature_type(const char * t)
{
    if (strcmp(t, feature_type_CpGI) == 0)
        return FT_CPGI;
    else if (strcmp(t, feature_type_gene) == 0)
        return FT_GENE;
    return FT_OTHER;
}

static int CpGIOverlap_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    /*   Buffer upon encountering a gene entry until we 
     *   determine whether a CpGI overlapts it's TSS
     *   That is, buffer it until we encounter the next
     *   CpGI or until we marked it as being overlapped
     *   with the last discovered CpGI or until end of stream
     *
     */
    overlap_feature_type_t feature_type;
    CpGIOverlap_stream * overlap_stream;
    GtGenomeNode * cur_node;
    int err_num = 0;
    *gn = NULL;


    overlap_stream = CpGIOverlap_stream_cast(ns);

    // if we're currently getting rid of buffer backlog dump it out
    //  until the buffer is empty
    if (overlap_stream->emptying_buffer && 
        gt_array_size(overlap_stream->node_buffer) > 0
       )
    {
        //ensure array has been reverse before doing this
        *gn = gt_array_pop(overlap_stream->node_buffer);
        return;
    }
    else
    {
        overlap_stream->emptying_buffer = 0;
    }

    // find genes. If not a gene send on its way, possibly saving CpGI
    //  location.

    while (*gn == NULL)
    {
        if(!gt_node_stream_next(overlap_stream->in_stream,
                               &cur_node,
                               err
                              ) && cur_node != NULL
          )
        {
            // try casting as a feature node so we can test type
            if(!gt_genome_node_try_cast(gt_feature_node_class(), cur_node))
            {
                if (overlap_stream->building_buffer)
                   gt_array_add(overlap_stream->node_buffer, cur_node);
                else
                  *gn = cur_node;
            }
            else // we found a feature node
            {
                const char * szType = gt_feature_node_get_type(cur_node);
                if (szType = NULL)
                {
                    if (overlap_stream->building_buffer)
                       gt_array_add(overlap_stream->node_buffer, cur_node);
                    else
                       *gn = cur_node;
                    break;
                }
                else // the feature node has an associated type
                {
                    switch(determine_feature_type(szType))
                    {
                    case FT_GENE:
                        overlap_stream->latest_gene_node = cur_node;
                        overlap_stream->latest_tss = (gt_feature_node_get_strand(cur_node) == GT_STRAND_FORWARD) ? gt_genome_node_get_start(cur_node) : 
                        gt_genome_node_get_end(cur_node);
                        break;
                    case FT_CPGI:
                        overlap_stream->island.start = gt_genome_node_get_start(cur_node);
                        overlap_stream->island.end   = gt_genome_node_get_end(cur_node);
                        overlap_stream->latest_CpGI_node = cur_node;
                        break;
                    default:
                        if (overlap_stream->building_buffer)
                            gt_array_add(overlap_stream->node_buffer, cur_node);
                        else
                            *gn = cur_node;
                        break;
                    }
                }
            }
        }
        else // we couldn't get anything else out of the in stream
        {
            return 1;
        }
    }

    return err_num;
}

static void CpGIOverlap_stream_free(GtNodeStream * ns)
{
    CpGIOverlap_stream * overlap_stream;
    
    overlap_stream = CpGIOverlap_stream_cast(ns);
    gt_array_delete(overlap_stream->node_buffer);
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

    overlap_stream->node_buffer = gt_array_new(sizeof(GtGenomeNode *));
    overlap_stream->index = 0;

    // we haven't found an island yet
    overlap_stream->island.start = 0;
    overlap_stream->island.end   = 0;
    overlap_stream->latest_tss   = 0;
 
    // stream eof hasn't been found yet
    overlap_stream->emptying_buffer = 0;

    // we aren't buffering for status of gene overlap
    overlap_stream->building_buffer = 0;

    // haven't found a cpGI or gene yet
    overlap_stream->latest_CpGI_node = NULL;
    overlap_stream->latest_gene_node = NULL;

    return ns;
 
}
