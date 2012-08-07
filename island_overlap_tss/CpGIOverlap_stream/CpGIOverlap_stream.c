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
#include "CpGIOverlap_stream.h"


struct CpGIOverlap_stream {
    const GtNodeStream parent_instance;
    GtNodeStream * in_stream;
    FILE * cpgi_file;

};

static const char * feature_type_gene = "gene";


const GtNodeStreamClass * CpGIOverlap_stream_class(void);

#define CpGIOverlap_stream_cast(GS) gt_node_stream_cast(CpGIOverlap_stream_class(), GS);


static inline int in_range(unsigned long n, unsigned long s, unsigned long e)
{
    return (n >= s && n <= e) ? 1 : 0;
}

// not thread safe
const char * CpGIOverlap_stream_find_gene_overlap( CpGIOverlap_stream * context,
                                                   unsigned long        TSS,
                                                   int                  chromosome
                                                 )
{
    static char found_island_name[255];
    unsigned long island_start = 0;
    unsigned long island_end   = 0;
    int           island_chromosome = 0;
    char          buf;
    int           err;

    int passed_TSS = false;
    
    // file positioning
    long int file_start_search_pos;
    file_start_search_pos = ftell(context->cpgi_file);

    // to find the associated CpGI, we check for valid islands moving forward in file until we get beyond the TSS.  If
    // none is found, we move backwards in the file until we get beyond the TSS and then declare an overlap not found

    while (island_chromosome < chromosome || (island_chromosome == chromosome && island_end <= TSS))
    {
        if (4 != fscanf(context->cpgi_file, "%s %d %lu %lu", found_island_name, &island_chromosome, &island_start, &island_end) )
            break; // reached EOF
        if (island_chromosome != chromosome)
            continue; // keep scanning forward until we get to our chromosome
        if (in_range(TSS, island_start, island_end))
            return found_island_name;
    }

    // if we didn't find an overlapping cpgi moving forward, rewind back to where we were and search file backwards
    fseek(context->cpgi_file, file_start_search_pos + 1, SEEK_SET);
    

    // fake this info so we can search 
    island_chromosome = chromosome;
    island_start = TSS + 1;
    island_end   = TSS + 1;
 
    while ( island_chromosome > chromosome || (island_chromosome == chromosome && island_start >= TSS))
    {
        // we have to search backwards from the current position for a newline marker
        while (! ( err = fseek(context->cpgi_file, -2, SEEK_CUR)))
        {
           fread(&buf, 1, 1, context->cpgi_file);
           if (buf == '\n') // look for newline
               break;
        }
       
        // if error we were at beginning of file and it's time to return without finding anything
        if (err)
            return (const char *)0;

        file_start_search_pos = ftell(context->cpgi_file);

        // now fscanf, check for a match, if no match rewind to the beginning of the search line
        if (4 != fscanf(context->cpgi_file, "%s %d %lu %lu", found_island_name, &island_chromosome, &island_start, &island_end) )
            break; // something went wrong

        if (chromosome == island_chromosome && in_range(TSS, island_start, island_end))
        {
             printf("Found special");
             return found_island_name;
        }

        fseek(context->cpgi_file, file_start_search_pos, SEEK_SET);
    }

    return (const char *)0;
}

static int CpGIOverlap_stream_next(GtNodeStream * ns,
                                   GtGenomeNode ** gn,
                                   GtError * err)
{
    GtGenomeNode * cur_node, * next_node;
    GtFeatureNodeIterator * iter;
    int err_num = 0;
    *gn = NULL;
    CpGIOverlap_stream * context;
    const char * gene_name = NULL;
    const char * overlap_name = NULL;
    char  chr_str[255];
    int  chr_num;
    unsigned int TSS;

    float CpGIOverlap;


    context = CpGIOverlap_stream_cast(ns);

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

              if ( 1 != sscanf(gt_str_get(gt_genome_node_get_seqid(cur_node)), "Chr%d", &chr_num))
                  return 0;

              TSS = (gt_feature_node_get_strand(cur_node) == GT_STRAND_FORWARD) ? gt_genome_node_get_start(cur_node) : gt_genome_node_get_end(cur_node);

              // now figure out the overlapping gene 
              if (! (overlap_name = CpGIOverlap_stream_find_gene_overlap( context, TSS, chr_num)))
                 return 0;

              // save the score into the node
              gt_feature_node_set_attribute(cur_node, "cpgi_at_tss", overlap_name);
              
              return 0;

         }
     }

    return err_num;
}

static void CpGIOverlap_stream_free(GtNodeStream * ns)
{
    CpGIOverlap_stream * score_stream;
    
    score_stream = CpGIOverlap_stream_cast(ns);
    fclose(score_stream->cpgi_file);
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

GtNodeStream * CpGIOverlap_stream_new(GtNodeStream * in_stream, const char * cpgi_db)
{
    GtNodeStream * ns = gt_node_stream_create(CpGIOverlap_stream_class(), 
                                              true); // must be sorted
    CpGIOverlap_stream * context = CpGIOverlap_stream_cast(ns);
    gt_assert(in_stream);
    context->in_stream = gt_node_stream_ref(in_stream);

    if ((context->cpgi_file = fopen(cpgi_db, "r")) == NULL)
    {
       gt_node_stream_delete(ns);
       fprintf(stderr, "Failed to open CpG Island db file %s\n", cpgi_db);
       return NULL;
    }

    return ns;
 
}
