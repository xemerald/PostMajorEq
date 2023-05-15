/*
 *
 */
#pragma once
/* */
#include <sachead.h>

#define SAC_FILE_NAME_FORMAT  "%s/%s.%s.%s.%s"

/* */
int sac_proc_station_data_extract(
	const char *, const char *, const char *, const char *[], const char *, const float [], float *[], int *, float *, double *
);
int sac_proc_sac_load( const char *, struct SAChead *, float ** );
struct SAChead *sac_proc_scnl_modify( struct SAChead *, const char *, const char *, const char *, const char * );
const char *sac_proc_scnl_print( struct SAChead * );
double sac_proc_reftime_fetch( struct SAChead * );
void sac_proc_data_preprocess( float *, const int, const float );
