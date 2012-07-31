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
#include <fstream>
#include <iostream>
#include <string>

std::string locate_island_score(std::ifstream * island_stream, std::string name)
{
 
    // this is really inneficient and needs to be fixed  
 
    island_stream->seekg(0, std::ios::beg);
    while (*island_stream)
    {
        std::string cpgi_name;
        std::string line;
   
        std::getline(*island_stream, line);
        boost::char_separator<char> sep(", \t", "", boost::drop_empty_tokens);
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator cur = tok.begin();

        std::string value;

        // don't process empty line
        if (cur == tok.end())
           continue;

        value = *cur;

        // see if the island matched
        if (value != name)
            continue;

        cur++; // go to next field
        
        if (cur == tok.end())
           continue;

        return *cur;  

    }
    
    return "NA";
}


int main(int argc, char ** argv)
{

    if (argc < 3)
    {
       std::cout << "Usage: " << argv[0] << " <gene file> <island file>" << std::endl;
       return 0;
    }

    std::ifstream gene_file(argv[1]);

    if (!gene_file.is_open())
    {
        std::cout << "Error opening gene file" << std::endl;
        return 0;
    } 

    std::ifstream island_file(argv[2]);

    if (!island_file.is_open())
    {
        std::cout << "Error opening island file" << std::endl;
        return 0;
    } 

    while (gene_file)
    {
        int col = 0;
        bool methylation_found  = false;
        std::string island_name;     
        std::string expression_level;
        std::string line;
        std::string  methylation_score;


        std::getline(gene_file, line);
        boost::char_separator<char> sep(", \t", "", boost::drop_empty_tokens);
        boost::tokenizer<boost::char_separator<char> > tok(line, sep);
        boost::tokenizer<boost::char_separator<char> >::iterator cur = tok.begin();

        for (; cur!= tok.end(); cur++)
        { 
            col++;
            
            // first column is expression level
            if (col == 1)
            {
                expression_level = *cur;
            }

            // we y want to find the CpGI
            if (col == 3)
            {
                island_name = *cur;
//                std::cout << island_name << std::endl;
                methylation_score = locate_island_score(&island_file, island_name);
                if (methylation_score != "NA")
                    methylation_found = true;
                break;
            }
        }

       if (methylation_found)
           std::cout << methylation_score << "\t"
                     << expression_level  <<  std::endl;
    }
   return 0;
}
