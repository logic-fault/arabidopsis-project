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
#include "genometools.h"	
#include <stdio.h>


void usage(const char * name)
{
   printf("Usage: %s <fileName>\n", name);
}

inline int in_range(unsigned long num, unsigned long min, unsigned long max)
{
   return (num >= min && num <= max) ? 1 : 0;
}

int main(int argc, char ** argv)
{
    int i;
    char * type;
    char * seqid;
    GtError *err;
    GtNodeStream * in_stream;
    GtGenomeNode * node, *latest_gene, *latest_cpgi;
    GtFeatureNodeIterator* feat_iter;
    GtFeatureNode * cur_feat;
    GtFeatureIndex *feature_index;
    GtArray * feature_array;

   
    unsigned long latest_tss        = 0; // tss for the last entry gene we read
    unsigned long latest_cpgi_start = 0;
    unsigned long latest_cpgi_end   = 0;

    int err_num;

    if (argc != 2)
    {
       usage(argv[0]);
       exit(1);
    }

    // initilaize genometools
    gt_lib_init();
    err = gt_error_new();

    feature_index = gt_feature_index_memory_new();

    if (gt_feature_index_add_gff3file(feature_index, argv[1], err))
    {
       printf("Falied gt_feature_index_add_gff3file %s\n", argv[1]);
    }



    seqid = gt_feature_index_get_first_seqid(feature_index);
    feature_array = gt_feature_index_get_features_for_seqid(feature_index, seqid);
    printf("Seqid = %s\n size =%d", seqid, gt_array_size(feature_array));

    for (i = 0; i < gt_array_size(feature_array) && i < 2; i++)
    {
       cur_feat = gt_array_get(feature_array,i);
       printf("%s\n", gt_feature_node_get_type(cur_feat));
    }

    // close genome tools
    gt_feature_index_delete(feature_index);
    gt_error_delete(err);
    gt_lib_clean();
    return 0;
}
