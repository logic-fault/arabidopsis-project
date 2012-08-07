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
#include "island_nuc_score_stream.h"


typedef struct
{
    unsigned long start;
    unsigned long end;
} island_t;

struct island_nuc_score_stream {
    const GtNodeStream parent_instance;
    GtNodeStream * in_stream;
    FILE * nucleosome_file;
  
    // we store these in case we fscanf'd an entry too far
    unsigned long previous_nucleosome_position;
    float         previous_nucleosome_reads;
    int           previous_nucleosome_chromosome;
};

static const char * feature_type_CpGI = "CpGI";


const GtNodeStreamClass * island_nuc_score_stream_class(void);

#define island_nuc_score_stream_cast(GS) gt_node_stream_cast(island_nuc_score_stream_class(), GS);

static float island_nuc_score_stream_score_island(island_nuc_score_stream * context,
                                            int island_chromosome_num,
                                            unsigned long island_start,
                                            unsigned long island_end
                                         )
{
    // iterate through the methylome db to find all entries
    // score is sum(entries in island range) / (num_cg)
    //
    // assume both streams are sorted

    unsigned long position = 0;
    int chromosome_num = island_chromosome_num;
    float reads;
    float score = 0.0f;


    // first see if the previous data matches this CpGI
    if (context->previous_nucleosome_chromosome == island_chromosome_num &&
        island_start <= context->previous_nucleosome_position &&
        island_end   >= context->previous_nucleosome_position
       )
          score += context->previous_nucleosome_reads; 

    while ((position < island_end && chromosome_num == island_chromosome_num) || chromosome_num < island_chromosome_num)
    {
       if (3 != fscanf(context->nucleosome_file, "%d %lu %f", &chromosome_num, &position, &reads))
           break;
       if (position >= island_start && position <= island_end && island_chromosome_num == chromosome_num)
           score += reads;
    }

    context->previous_nucleosome_position = position;
    context->previous_nucleosome_reads = reads;
    context->previous_nucleosome_chromosome = chromosome_num;

    score = score / (float)(island_end - island_start + 1);

    return score;
}

static int island_nuc_score_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    GtGenomeNode * cur_node;
    int err_num = 0;
    *gn = NULL;
    island_nuc_score_stream * score_stream;
    unsigned long island_start;
    unsigned long island_end;
    float island_score;
    int chromosome_num;
    GtStr * seqID_gtstr;
    char *  seqID_str;
    char *  num_cg_str;
    char score_str[255];

    score_stream = island_nuc_score_stream_cast(ns);

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

              num_cg_str = gt_feature_node_get_attribute(cur_node, "sumcg");
              
              // now figure out the score
              island_score = island_nuc_score_stream_score_island(score_stream ,
                                                            chromosome_num,
                                                            island_start,
                                                            island_end);

              sprintf(score_str, "%f", island_score);

              // save the score into the node
              gt_feature_node_set_attribute(cur_node, "nuc_density",score_str); 
              return 0;

         }
     }

    return err_num;
}

static void island_nuc_score_stream_free(GtNodeStream * ns)
{
    island_nuc_score_stream * score_stream;
    
    score_stream = island_nuc_score_stream_cast(ns);
    fclose(score_stream->methylome_file);
    return;
}

const GtNodeStreamClass * island_nuc_score_stream_class(void)
{
    static const GtNodeStreamClass * c = NULL;

    if (!c)
    {	
        c = gt_node_stream_class_new( sizeof(island_nuc_score_stream),
                                      island_nuc_score_stream_free,
                                      island_nuc_score_stream_next
                                    );
    }
    
    return c;
}

GtNodeStream * island_nuc_score_stream_new(GtNodeStream * in_stream, const char * nucleosome_db)
{
    GtNodeStream * ns = gt_node_stream_create(island_nuc_score_stream_class(), 
                                              true); // must be sorted
    island_nuc_score_stream * score_stream = island_nuc_score_stream_cast(ns);
    gt_assert(in_stream);
    score_stream->in_stream = gt_node_stream_ref(in_stream);
    score_stream->previous_nucleosome_position = 0;
    score_stream->previous_nucleosome_reads = 0.0f;
    score_stream->previous_nucleosome_chromosome = 0;

    if ((score_stream->nucleosome_file = fopen(nucleosome_db, "r")) == NULL)
    {
       gt_node_stream_delete(ns);
       fprintf(stderr, "Failed to open nucleosome db file %s\n", methylome_db);
       return NULL;
    }

    return ns;
 
}
