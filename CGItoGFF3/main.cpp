/*      
 *    @file 
 *    @author Brock Anderson <brock.wright.anderson@gmail.com>
 *    @section LICENSE
 *       Released to public domain without restriction
 *
 *    @section DESCRIPTION
 *       convert CPGISLE format to GFF3
 *       TODO: insert into current GFF3 document
 *       This code adapted from sample BOOST library
 *       code located at http://www.boost.org/doc/libs/1_50_0/more/getting_started/unix-variants.html#id33
 *
 *************************************************/


#include <boost/regex.hpp>
#include <iostream>
#include <string>

const std::string regex_island = "IS";

int main(int argc, char ** argv)
{
    std::string line;
    boost::regex pat(regex_island);

    while (std::cin)
    {
        std::getline(std::cin, line);
        boost::smatch matches;
        if (boost::regex_match(line, matches, pat))
            std::cout << matches[2] << std::endl;
    }
   return 0;
}
