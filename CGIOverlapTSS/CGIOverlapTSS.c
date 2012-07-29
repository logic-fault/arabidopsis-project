/*
* @file
* @author Brock Anderson <brock.wright.anderson@gmail.com>
* @section LICENSE
* Released to public domain without restriction
*
* @section DESCRIPTION
* Find CpGI that overlap TSS, mark as such
*
*************************************************/
#include "genometools.h"	
#include "CpGIOverlap_stream/CpGIOverlap_stream_api.h"
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
    GtNodeStream * in, * overlap, * out;
    GtFile * out_file;
    GtError * err;

    if (argc != 2)
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

    if (!(overlap = CpGIOverlap_stream_new(in)))
    {
        fprintf(stderr, "Failed to create CpGI overlap stream\n");
        exit(1);
    }

    if (!(out = gt_gff3_out_stream_new(overlap, out_file)))
    {
        fprintf(stderr, "Failed to create output stream\n");
        exit(1);
    }


    // close genome tools
    gt_error_delete(err);
    gt_lib_clean();
    return 0;
}
