/* 
 * @author brock a<brock.wright.anderson@gmail.cm>
 *
 *   
 *
 *
 */

#ifndef  ISLAND_NUC_SCORE_STREAM_API_H
#define  ISLAND_NUC_SCORE_STREAM_API_H

typedef struct island_nuc_score_stream island_nuc_score_stream;

GtNodeStream* island_nuc_score_stream_new(GtNodeStream * in_stream, const char * methylome_db);

#endif
