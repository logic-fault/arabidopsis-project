/*      
 *    @file 
 *    @author Brock Anderson <brock.wright.anderson@gmail.com>
 *    @section LICENSE
 *       Released to public domain without restriction
 *
 *    @section DESCRIPTION
 *       convert CPGISLE format to GFF3
 *       TODO: insert into current GFF3 document
 *
 *************************************************/

#include <boost/tokenizer.hpp>
#include <iostream>
#include <string>


typedef enum { ID_FIELD, FT_FIELD, UNKNOWN_FIELD } field_t;

int main(int argc, char ** argv)
{
    std::string seqid;
    std::string line;
    std::string start_base, end_base;

    std::cout << "##gff-version 3" << std::endl;

    bool waiting_for_sum = false;

    while (std::cin)
    {
        field_t field = UNKNOWN_FIELD;
        int col = 1;

        bool has_cpg_coords  = false;
        bool cpg_coords_found = false;
        bool sum_found = false;

        size_t dot_loc;
        std::getline(std::cin, line);
        boost::tokenizer<> tok(line);
        boost::tokenizer<>::iterator cur = tok.begin();
 
   

        if (tok.begin() != tok.end())
        {
           if (cur->compare("ID") == 0)
              field = ID_FIELD;
           else if (cur->compare("FT") == 0)
              field = FT_FIELD;
           else
              continue; // bad field
        }
        else
          continue;
 
        cur++;

        for (; cur!= tok.end(); cur++)
        { 
            col++;
            switch(field)
            {
            case FT_FIELD:
                if (col == 2 && !cur->compare("CpG"))
                {
                   has_cpg_coords = true;
                   break;
                }
                if (col == 2 && !cur->compare("Sum"))
                {
                   sum_found = true;
                }
                if (col == 4 && has_cpg_coords)
                {
                   start_base = *cur;
                }
                if (col == 5 && has_cpg_coords)
                {
                   end_base = *cur;
                   cpg_coords_found = true;
                   waiting_for_sum  = true;
                }
                if (col == 5 && waiting_for_sum && sum_found)
                {
                   std::cout << "sumcg=" << *cur << std::endl;
                   waiting_for_sum = false;
                }
                break;
            case ID_FIELD:
                if (col == 2)
                   seqid = *cur;
                break;
            }
        }
       if (cpg_coords_found)
           std::cout << seqid << "\t"
                     << "."   << "\t" 
                     << "CpGI" << "\t"
                     << start_base << "\t" 
                     << end_base <<"\t"
                     << "."      << "\t"
                     << "."      << "\t"
                     << "."      << "\t";
                     //wait for annotations until sumcg known
    }
   return 0;
}
