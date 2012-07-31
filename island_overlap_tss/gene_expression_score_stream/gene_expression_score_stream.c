/*
* @file
* @author Brock Anderson <brock.wright.anderson@gmail.com>
* @section LICENSE
* Released to public domain without restriction
*
* @section DESCRIPTION
*   Score genes based upon input rna-seq db
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

};

static const char * feature_type_gene = "gene";


const GtNodeStreamClass * gene_expression_score_stream_class(void);

#define gene_expression_score_stream_cast(GS) gt_node_stream_cast(gene_expression_score_stream_class(), GS);

static float gene_expression_score_stream_score_gene( gene_expression_score_stream * context,
                                                     const char * gene_name
                                                    )
{
    char found_name[255];
    char trash_buffer[255];
    float found_expression;
    float score = 0.0f;

    // this is slow but we rewind the rna-seq each time because it's out of sequence
    rewind(context->rnaseq_file);

    while (3 == fscanf(context->rnaseq_file, "%s %s %f", found_name, trash_buffer, &found_expression) )
    {
        if (strcmp(found_name, gene_name) == 0)
            score += found_expression;
    }
    return score;
}

static int gene_expression_score_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    GtGenomeNode * cur_node, * next_node;
    GtFeatureNodeIterator * iter;
    int err_num = 0;
    *gn = NULL;
    gene_expression_score_stream * context;
    const char * gene_name = NULL;

    float gene_expression_score;


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
              // first check if it is a pseudo node, if so find the gene in it if available
              if (gt_feature_node_is_pseudo(cur_node))
              {
                  iter = gt_feature_node_iterator_new(cur_node);
                  if (iter == NULL)
                      return;
                  while ((next_node = gt_feature_node_iterator_next(iter)) && !gt_feature_node_has_type(next_node, feature_type_gene));
                  gt_feature_node_iterator_delete(iter);
                  if (NULL == (cur_node = next_node))
                     return 0;
              }


              if(!gt_feature_node_has_type(cur_node, feature_type_gene))
                  return 0;

              // find name of gene
              gene_name = gt_feature_node_get_attribute(cur_node, "Name");

              if (gene_name == NULL)
                  return;

              // now figure out the score
              gene_expression_score = gene_expression_score_stream_score_gene(context, gene_name);

              // save the score into the node
              gt_feature_node_set_score(cur_node, gene_expression_score);
              
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
