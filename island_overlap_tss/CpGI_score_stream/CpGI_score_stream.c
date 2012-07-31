/*
* @file
* @author Brock Anderson <brock.wright.anderson@gmail.com>
* @section LICENSE
* Released to public domain without restriction
*
* @section DESCRIPTION
* Find CpGI and score them based upon methylome db
*
*
*************************************************/
#include <genometools.h>
#include <stdio.h>
#include <stdlib.h>
#include "CpGI_score_stream.h"


typedef struct
{
    unsigned long start;
    unsigned long end;
} island_t;

struct CpGI_score_stream {
    const GtNodeStream parent_instance;
    GtNodeStream * in_stream;
    FILE * methylome_file;
  
    // we store these in case we fscanf'd an entry too far
    unsigned long previous_methylome_position;
    float         previous_methylome_fraction;
    int           previous_methylome_chromosome;
};

static const char * feature_type_CpGI = "CpGI";


const GtNodeStreamClass * CpGI_score_stream_class(void);

#define CpGI_score_stream_cast(GS) gt_node_stream_cast(CpGI_score_stream_class(), GS);

static float CpGI_score_stream_score_island(CpGI_score_stream * context,
                                            int island_chromosome_num,
                                            unsigned long island_start,
                                            unsigned long island_end
                                         )
{
    // iterate through the methylome db to find all entries
    // score is sum(entries in island range) / (island_end - island_start + 1)
    //
    // assume both streams are sorted

    unsigned long position = 0;
    int chromosome_num = island_chromosome_num;
    float methylation;
    float score = 0.0f;

    // first see if the previous data matches this CpGI
    if (context->previous_methylome_chromosome == island_chromosome_num &&
        island_start <= context->previous_methylome_position &&
        island_end   >= context->previous_methylome_position
       )
          score += context->previous_methylome_fraction; 

    while ((position < island_end && chromosome_num == island_chromosome_num) || chromosome_num < island_chromosome_num)
    {
       if (3 != fscanf(context->methylome_file, "%d %lu %f", &chromosome_num, &position, &methylation))
           break;
       if (position >= island_start && position <= island_end && island_chromosome_num == chromosome_num)
           score += methylation;
    }

    context->previous_methylome_position = position;
    context->previous_methylome_fraction = methylation;
    context->previous_methylome_chromosome = chromosome_num;

    score = score / (float)(island_end - island_start + 1);

    return score;
}

static int CpGI_score_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    GtGenomeNode * cur_node;
    int err_num = 0;
    *gn = NULL;
    CpGI_score_stream * score_stream;
    unsigned long island_start;
    unsigned long island_end;
    float island_score;
    int chromosome_num;
    GtStr * seqID_gtstr;
    char *  seqID_str;

    score_stream = CpGI_score_stream_cast(ns);

    // find the CpGI's, process methylome score
     if(!gt_node_stream_next(score_stream->in_stream,
                            &cur_node,
                            err
                           ) && cur_node != NULL
       )
     {
         *gn = cur_node;

         // try casting as a feature node so we can test type
         if(!gt_genome_node_try_cast(gt_feature_node_class(), cur_node))
         {
               return 0;
         }
         else // we found a feature node
         {
              if(!gt_feature_node_has_type(cur_node, feature_type_CpGI))
                  return 0;

              #if DEBUG_SCORE
              printf("found CpGI\n");
              #endif 
 
              island_start = gt_genome_node_get_start(cur_node);
              island_end   = gt_genome_node_get_end(cur_node);

              seqID_gtstr = gt_genome_node_get_seqid(cur_node);
              seqID_str   = gt_str_get(seqID_gtstr);

              sscanf(seqID_str, "Chr%d", &chromosome_num);

              // now figure out the score
              island_score = CpGI_score_stream_score_island(score_stream ,
                                                            chromosome_num,
                                                            island_start,
                                                            island_end);
//              gt_str_delete(seqID_gtstr);

              // save the score into the node
              gt_feature_node_set_score(cur_node, island_score);
              
              return 0;

         }
     }

    return err_num;
}

static void CpGI_score_stream_free(GtNodeStream * ns)
{
    CpGI_score_stream * score_stream;
    
    score_stream = CpGI_score_stream_cast(ns);
    fclose(score_stream->methylome_file);
    return;
}

const GtNodeStreamClass * CpGI_score_stream_class(void)
{
    static const GtNodeStreamClass * c = NULL;

    if (!c)
    {	
        c = gt_node_stream_class_new( sizeof(CpGI_score_stream),
                                      CpGI_score_stream_free,
                                      CpGI_score_stream_next
                                    );
    }
    
    return c;
}

GtNodeStream * CpGI_score_stream_new(GtNodeStream * in_stream, const char * methylome_db)
{
    GtNodeStream * ns = gt_node_stream_create(CpGI_score_stream_class(), 
                                              true); // must be sorted
    CpGI_score_stream * score_stream = CpGI_score_stream_cast(ns);
    gt_assert(in_stream);
    score_stream->in_stream = gt_node_stream_ref(in_stream);
    score_stream->previous_methylome_position = 0;
    score_stream->previous_methylome_fraction = 0.0f;
    score_stream->previous_methylome_chromosome = 0;

    if ((score_stream->methylome_file = fopen(methylome_db, "r")) == NULL)
    {
       gt_node_stream_delete(ns);
       fprintf(stderr, "Failed to open methylome db file %s\n", methylome_db);
       return NULL;
    }

    return ns;
 
}
