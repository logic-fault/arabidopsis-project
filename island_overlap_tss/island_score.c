/*
* @file
* @author Brock Anderson <brock.wright.anderson@gmail.com>
* @section LICENSE
* Released to public domain without restriction
*
* @section DESCRIPTION
*  score the CpGI based upon methylome
*
*************************************************/
#include "genometools.h"	
#include "CpGI_score_stream/CpGI_score_stream_api.h"
#include <stdio.h>


void usage(const char * name)
{
   printf("Usage: %s <in fileName> <out fileName> <methylome db>\n", name);
}


int main(int argc, char ** argv)
{
    GtNodeStream * in, * score, * out;
    GtFile * out_file;
    GtError * err;

    if (argc != 4)
    {
       usage(argv[0]);
       exit(1);
    }

    // initilaize genometools
    gt_lib_init();
    err = gt_error_new();

    if (!(in = gt_gff3_in_stream_new_sorted(argv[1])))
    {
        fprintf(stderr, "Failed to open input stream with arg %s\n", argv[1]);
        exit(1);
    }

    if (!(out_file = gt_file_new(argv[2], "w+", err)))
    {
        gt_node_stream_delete(in);
        fprintf(stderr, "Failed to create output file %s\n", argv[2]);
        exit(1);
    }

    if (!(score = CpGI_score_stream_new(in, argv[3])))
    {

        gt_file_delete(out_file);
        gt_node_stream_delete(in);
        fprintf(stderr, "Failed to create CpGI score stream\n");
        exit(1);
    }
    out = gt_gff3_out_stream_new(in, out_file);
    
    if (!(out = gt_gff3_out_stream_new(score, out_file)))
    {
        gt_node_stream_delete(score);
        gt_file_delete(out_file);
        gt_node_stream_delete(in);
        fprintf(stderr, "Failed to create output stream\n");
        exit(1);
    }

    if (gt_node_stream_pull(out, err))
    {
        fprintf(stderr, "Failed to pull through out stream\n");
    }

    // close genome tools
    gt_node_stream_delete(out);
    gt_node_stream_delete(score);
    gt_file_delete(out_file);
    gt_node_stream_delete(in);
    gt_error_delete(err);
    gt_lib_clean();
    return 0;
}
