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
#include "gene_expression_score_stream.h"


struct gene_expression_score_stream {
    const GtNodeStream parent_instance;
    GtNodeStream * in_stream;
    FILE * rnaseq_file;
  
    // we store these in case we fscanf'd an entry too far
};

static const char * feature_type_gene = "Gene";


const GtNodeStreamClass * gene_expression_score_stream_class(void);

#define gene_expression_score_stream_cast(GS) gt_node_stream_cast(CpGI_score_stream_class(), GS);

static float gene_expression_score_stream_score_gene(gene_expression_score_stream * context,
                                                     const char * gene_name;
                                                    )
{
    float score = 0.0f;

    return score;
}

static int gene_expression_score_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    GtGenomeNode * cur_node;
    int err_num = 0;
    *gn = NULL;
    gene_expression_score_stream * score_stream;
    const char * gene_name = NULL;


    context = gene_expression_score_stream_cast(ns);

    // find the genes, determine expression level
     if(!gt_node_stream_next(context->in_stream,
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

 
              island_start = gt_genome_node_get_start(cur_node);
              island_end   = gt_genome_node_get_end(cur_node);

              seqID_gtstr = gt_genome_node_get_seqid(cur_node);
              seqID_str   = gt_str_get(seqID_gtstr);

              sscanf(seqID_str, "Chr%d", &chromosome_num);

              // now figure out the score
              island_score = gene_expression_score_stream_score_island(score_stream ,
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

static void gene_expression_score_stream_free(GtNodeStream * ns)
{
    gene_expression_score_stream * score_stream;
    
    score_stream = gene_expression_score_stream_cast(ns);
    fclose(score_stream->rnaseq_file);
    return;
}

const GtNodeStreamClass * gene_expression_score_stream_class(void)
{
    static const GtNodeStreamClass * c = NULL;

    if (!c)
    {	
        c = gt_node_stream_class_new( sizeof(gene_expression_score_stream),
                                      gene_expression_score_stream_free,
                                      gene_expression_score_stream_next
                                    );
    }
    
    return c;
}

GtNodeStream * gene_expression_score_stream_new(GtNodeStream * in_stream, const char * rnaseq_db)
{
    GtNodeStream * ns = gt_node_stream_create(gene_expression_score_stream_class(), 
                                              true); // must be sorted
    gene_expression_score_stream * context = gene_expression_score_stream_cast(ns);
    gt_assert(in_stream);
    context->in_stream = gt_node_stream_ref(in_stream);

    if ((context->rnaseq_file = fopen(rnaseq_db, "r")) == NULL)
    {
       gt_node_stream_delete(ns);
       fprintf(stderr, "Failed to open RNA seq db file %s\n", rnaseq_db);
       return NULL;
    }

    return ns;
 
}
