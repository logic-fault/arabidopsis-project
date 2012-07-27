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
*                  second_latest_tss <- latest_tss
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

int main(int argc, char ** argv)
{
    GtError *err;
    GtNodeStream * in_stream;
    GtGenomeNode * node;
   
    unsigned long latest_tss; // tss for the last entry gene we read
    unsigned long second_latest_tss; // tss for the second to last entry we read

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
      
    // parse the file
    
    // close genome tools
    gt_error_delete(err);
    gt_lib_clean();
    return 0;
}
