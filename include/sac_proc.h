/*
 *
 */
#pragma once
/* */
#include <sachead.h>

/* */
int sac_proc_station_data_extract( const char *, const char *, const float [], float *[], int *, float *, double * );
int sac_proc_sac_load( const char *, struct SAChead *, float ** );
struct SAChead *sac_proc_scnl_modify( struct SAChead *, const char *, const char *, const char *, const char * );
const char *sac_proc_scnl_print( struct SAChead * );
