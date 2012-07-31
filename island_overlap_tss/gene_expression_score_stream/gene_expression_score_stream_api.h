/* 
 * @author brock a<brock.wright.anderson@gmail.cm>
 *
 *   
 *
 *
 */

#ifndef  GENE_EXPRESSION_SCORE_STREAM_API_H
#define  GENE_EXPRESSION_SCORE_STREAM_API_H

typedef struct gene_expression_score_stream CpGI_score_stream;

GtNodeStream* gene_expression_score_stream_new(GtNodeStream * in_stream, const char * methylome_db);

#endif
