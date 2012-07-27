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
#include <stdio.h>

void usage(const char * name)
{
   printf("Usage: %s <fileName>\n", name);
}

int main(int argc, char ** argv)
{
    GtError *err;
    GtNodeStream * in_stream;
    GtGenomeNode * node;
    int err_num;

    if (argc != 2)
    {
       usage(argv[0]);
       exit(1);
    }

    // initilaize genometools
    gt_lib_init();
    err = gt_error_new();

    // open the gff3 file we are evaluating
    in_stream = gt_gff3_in_stream_new_sorted(argv[1]);  
    gt_gff3_in_stream_show_progress_bar((GtGFF3InStream *)in_stream);

    // close genome tools
    gt_error_delete(err);
    gt_lib_clean();
    return 0;
}
