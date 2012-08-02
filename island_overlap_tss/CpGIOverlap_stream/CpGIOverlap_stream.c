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
#include <stdio.h>
#include <stdlib.h>
#include "CpGIOverlap_stream.h"

#define MAX_CPGI_DIGITS 16
#define DEBUG_OVERLAP 0

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
    GtGenomeNode * latest_pseudo_node;
    GtArray      * node_buffer;
    island_t island;
    unsigned long latest_tss;
    int emptying_buffer;
    int building_buffer;
    int eof_found;
    unsigned long CpGI_number;
    char latest_cpgi_name[MAX_CPGI_DIGITS + 6]; // leave space for CpGI_
};

typedef enum { FT_CPGI, FT_GENE, FT_OTHER } overlap_feature_type_t;

const char * feature_type_CpGI = "CpGI";
const char * feature_type_gene = "gene";


const GtNodeStreamClass * CpGIOverlap_stream_class(void);

#define CpGIOverlap_stream_cast(GS) gt_node_stream_cast(CpGIOverlap_stream_class(), GS);

static inline int in_island(unsigned long tss, island_t island)
{
    return (tss >= island.start && tss <= island.end) ? 1 : 0;
}

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
    char cpgi_name_buf[MAX_CPGI_DIGITS + 6];
    overlap_feature_type_t feature_type;
    CpGIOverlap_stream * overlap_stream;
    GtGenomeNode * cur_node, * next_node, * pseudo_node;
    GtFeatureNode * rep;
    GtFeatureNodeIterator * iter;
    int err_num = 0;
    *gn = NULL;
    const char * szType;


    overlap_stream = CpGIOverlap_stream_cast(ns);

    #if DEBUG_OVERLAP
    printf("Before buffer emptying checks\n");
    #endif

    // if we're currently getting rid of buffer backlog dump it out
    //  until the buffer is empty
    if (overlap_stream->emptying_buffer && 
        gt_array_size(overlap_stream->node_buffer) > 0
       )
    {
        //ensure array has been reverse before doing this
        *gn = *(GtGenomeNode **)gt_array_pop(overlap_stream->node_buffer);
        return 0;
    }
    else if (overlap_stream->emptying_buffer && !overlap_stream->eof_found)
    {
        // if we were emptying a buffer take care of getting ready here
        overlap_stream->emptying_buffer = 0;

        // now put the latest gene into the buffer array if needed
        if (overlap_stream->building_buffer)
            gt_array_add(overlap_stream->node_buffer, overlap_stream->latest_pseudo_node);
    }
    else if (overlap_stream->eof_found)
    {
        // we've emptied buffer and nothing is left;
        *gn = NULL;
        return 0;
    }

    // find genes. If not a gene send on its way, possibly saving CpGI
    //  location.

    while (*gn == NULL)
    {

        #if DEBUG_OVERLAP
        printf("Loop state: buffer_size=%d\n", gt_array_size(overlap_stream->node_buffer));
        printf("eof_found = %d\n emptying_buffer = %d\n building_buffer=%d\n", overlap_stream->eof_found, overlap_stream->emptying_buffer, overlap_stream->building_buffer);
        #endif


        if(!gt_node_stream_next(overlap_stream->in_stream,
                               &cur_node,
                               err
                              ) && cur_node != NULL
          )
        {
            pseudo_node = cur_node;

            // try casting as a feature node so we can test type
            if(!gt_genome_node_try_cast(gt_feature_node_class(), cur_node))
            {
                #if DEBUG_OVERLAP
                printf("Not feature\n");
                #endif

                if (overlap_stream->building_buffer)
                   gt_array_add(overlap_stream->node_buffer, cur_node);
                else
                  *gn = cur_node;
            }
            else // we found a feature node
            {
                #if DEBUG_OVERLAP
                printf("node_start=%d\n", gt_genome_node_get_start(cur_node));
                printf("Feature\n");
                #endif

                //first thing we need to do is determine if we're dealing with a pseudo-node
                // if it's a pseudo node we save the gene as cur_node but want to 
                // use the pseudo node for purposes of buffering / sending on to the next
                // link in the stream

                if(gt_feature_node_is_pseudo(cur_node))
                {
                    iter = gt_feature_node_iterator_new(cur_node);
                    while (next_node = gt_feature_node_iterator_next(iter))
                    {
                        if ( (szType = gt_feature_node_get_type(next_node)) && determine_feature_type(szType) == FT_GENE)
                        {
                            #if DEBUG_OVERLAP
                            printf("Got gene within pseudo node");
                            #endif
                            cur_node = next_node;
                            break;
                        }
                    }
                    gt_feature_node_iterator_delete(iter);
                }

                szType = gt_feature_node_get_type(cur_node);
                if (szType == NULL)
                {
                    if (overlap_stream->building_buffer)
                       gt_array_add(overlap_stream->node_buffer, cur_node);
                    else
                       *gn = cur_node;
                    continue;
                }
                else // the feature node has an associated type
                {
                    switch(determine_feature_type(szType))
                    {
                    case FT_GENE:
                        #if DEBUG_OVERLAP
                        printf("In a GENE\n");
                        printf("Gene Id=%s\n", gt_feature_node_get_attribute(cur_node, "ID"));
                        #endif
                    

                        overlap_stream->latest_pseudo_node = pseudo_node;
                        overlap_stream->latest_gene_node   = cur_node;

                        overlap_stream->latest_tss = (gt_feature_node_get_strand(cur_node) == GT_STRAND_FORWARD) ? gt_genome_node_get_start(cur_node) : 
                        gt_genome_node_get_end(cur_node);
                        
                        #if DEBUG_OVERLAP
                        printf ("GENE COORD TSS=%d start=%d end=%d\n", overlap_stream->latest_tss, overlap_stream->island.start, overlap_stream->island.end);
                        #endif

                        // now we see if we can mark for overlap
                        if (in_island(overlap_stream->latest_tss, overlap_stream->island))
                        {
                           #if DEBUG_OVERLAP
                           printf("GENE found and marked\n");
                           #endif
                           gt_feature_node_set_attribute(overlap_stream->latest_gene_node,
                                                         "cpgi_at_tss",
                                                        gt_feature_node_get_attribute(overlap_stream->latest_CpGI_node, "Name")
                                                        );
                           // TODO: mark current node
                           // mark it, we're on our way
                           overlap_stream->building_buffer = 0;
                           overlap_stream->emptying_buffer = 1;
                           gt_array_add(overlap_stream->node_buffer, pseudo_node);
                           gt_array_reverse(overlap_stream->node_buffer);
                           *gn = *(GtGenomeNode **)gt_array_pop(overlap_stream->node_buffer);
                           return 0;
                        }
                        else // possible in next island
                        {
                            #if DEBUG_OVERLAP
                            printf("GENE found and not marked\n");
                            #endif

                            overlap_stream->building_buffer = 1;
                            if (gt_array_size(overlap_stream->node_buffer) > 0)
                            {
                                gt_array_reverse(overlap_stream->node_buffer);
                                *gn = *(GtGenomeNode **)gt_array_pop(overlap_stream->node_buffer);
                                overlap_stream->emptying_buffer = 1;
                            }
                            // if there wasn't a previous buffer start building this one
                            else
                            {
                                gt_array_add(overlap_stream->node_buffer, pseudo_node);
                                continue;
                            }
                        }
                        break;
                    case FT_CPGI:
                        #if DEBUG_OVERLAP
                        printf("CPGI found\n");
                        #endif
                        overlap_stream->latest_CpGI_node = cur_node;
                        (overlap_stream->CpGI_number)++;
                        sprintf(cpgi_name_buf,"CpGI_%lu", overlap_stream->CpGI_number);
                        
                        gt_feature_node_set_attribute(cur_node, 
                                                      "Name", 
                                                      cpgi_name_buf
                                                      );
 
                        // save island
                        overlap_stream->island.start = gt_genome_node_get_start(cur_node);
                        overlap_stream->island.end   = gt_genome_node_get_end(cur_node);

                        // TODO: we can only clear the buffer here if this island's end is greater than last gene's TSS
                        if (overlap_stream->building_buffer && overlap_stream->island.end < overlap_stream->latest_tss)
                        {
                           // this one didn't fit but we keep the buffer because another island may show up
                           // keep the record of the last gene
                           gt_array_add(overlap_stream->node_buffer, cur_node);
                           break;
                        }

                        // test previous gene, then clear the buffer
                        if (overlap_stream->building_buffer && in_island(overlap_stream->latest_tss, overlap_stream->island))
                        {
                           gt_feature_node_set_attribute(overlap_stream->latest_gene_node,
                                                         "cpgi_at_tss",
                                                        gt_feature_node_get_attribute(overlap_stream->latest_CpGI_node, "Name")
                                                        );
                           
                        }
                       
                        overlap_stream->building_buffer = 0;
                        overlap_stream->emptying_buffer = 1;
                        gt_array_add(overlap_stream->node_buffer, cur_node);
                        #if DEBUG_OVERLAP
                        printf("Before reverse\n");
                        #endif
                        gt_array_reverse(overlap_stream->node_buffer);
                        
                        #if DEBUG_OVERLAP
                        printf("buffer size before pop=%d\n", gt_array_size(overlap_stream->node_buffer));
                        #endif
                        *gn = *(GtGenomeNode **)gt_array_pop(overlap_stream->node_buffer);
                        break;
                    case FT_OTHER:
                    default:
                        #if DEBUG_OVERLAP
                        printf("FEATURE NOT CPGI OR GENE\n");
                        #endif

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
            overlap_stream->building_buffer = 0;
            overlap_stream->emptying_buffer = 1;
            overlap_stream->eof_found       = 1;
            if (gt_array_size(overlap_stream->node_buffer) > 0)
            {
               gt_array_reverse(overlap_stream->node_buffer);
               *gn = *(GtGenomeNode **)gt_array_pop(overlap_stream->node_buffer);
            }
            else
               *gn = NULL;
            return 0;
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
    overlap_stream->eof_found = 0;

    // we aren't buffering for status of gene overlap
    overlap_stream->building_buffer = 0;
    overlap_stream->emptying_buffer = 0;

    // haven't found a cpGI or gene yet
    overlap_stream->latest_pseudo_node = NULL;
    overlap_stream->latest_CpGI_node   = NULL;
    overlap_stream->latest_gene_node   = NULL;
    overlap_stream->CpGI_number        = 0;

    return ns;
 
}
