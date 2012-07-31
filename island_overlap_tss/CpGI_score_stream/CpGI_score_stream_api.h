/* 
 * @author brock a<brock.wright.anderson@gmail.cm>
 *
 *   
 *
 *
 */

#ifndef  CPGI_OVERLAP_STREAM_API_H
#define  CPGI_OVERLAP_STREAM_API_H

typedef struct CpGI_score_stream CpGI_score_stream;

GtNodeStream* CpGI_score_stream_new(GtNodeStream * in_stream, const char * methylome_db);

#endif
