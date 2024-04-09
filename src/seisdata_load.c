/**
 * @file seisdata_load.c
 * @author Benjamin Yang @ National Taiwan University (b98204032@gmail.com)
 * @brief
 * @version 1.0.1
 * @date 2024-04-07
 *
 * @copyright Copyright (c) 2024
 *
 */
/* */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
/* */
#include <postmajor.h>
#include <sac.h>

/* */
#define MAX_FILE_NAME          512
#define SAC_FILE_NAME_FORMAT  "%s/%s.%s.%s.%s"

/**
 * @brief
 *
 * @param snl_info
 * @param path
 * @return int
 */
int seisdata_load_sac( SNL_INFO *snl_info, const char *path )
{
	char   filename[MAX_FILE_NAME] = { 0 };
	struct SAChead sh;
	float *_seis = NULL;

/* Just a initialization */
	snl_info->npts      = -1;
	snl_info->delta     = -1.0;
	snl_info->starttime = -1.0;
	for ( int i = 0; i < NUM_CHANNEL_SNL; i++ )
		snl_info->seis[i] = NULL;

/* Open all the three channels' SAC files */
	for ( int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
	/* Opening the SAC files */
		sprintf(filename, SAC_FILE_NAME_FORMAT, path, snl_info->sta, snl_info->chan[i], snl_info->net, snl_info->loc);
		if ( sac_file_load( filename, &sh, &_seis ) < 0 )
			return -1;

	/* Check the consistency of npts */
		if ( snl_info->npts < 0 ) {
			snl_info->npts = (int)sh.npts;
		}
		else if ( snl_info->npts != (int)sh.npts ) {
			fprintf(
				stderr, "WARNING! There is a different npts within the SAC files of %s.%s.%s.%s.\n",
				snl_info->sta, snl_info->chan[i], snl_info->net, snl_info->loc
			);
			if ( (int)sh.npts < snl_info->npts )
				snl_info->npts = (int)sh.npts;
		}
	/* Check the consistency of delta */
		if ( snl_info->delta < 0.0 ) {
			snl_info->delta = sh.delta;
		}
		else if ( fabs(snl_info->delta - sh.delta) > FLT_EPSILON ) {
			fprintf(
				stderr, "ERROR! There is a different delta within the SAC files of %s.%s.%s.%s.\n",
				snl_info->sta, snl_info->chan[i], snl_info->net, snl_info->loc
			);
			free(_seis);
			return -2;
		}
	/* Check the consistency of start time of waveform */
		if ( snl_info->starttime < 0.0 ) {
			snl_info->starttime = sac_reftime_fetch( &sh );
		}
		else if ( fabs(snl_info->starttime - sac_reftime_fetch( &sh )) > FLT_EPSILON ) {
			fprintf(
				stderr, "ERROR! There is a different start time within the SAC files of %s.%s.%s.%s.\n",
				snl_info->sta, snl_info->chan[i], snl_info->net, snl_info->loc
			);
			free(_seis);
			return -2;
		}
	/* */
		sac_data_preprocess( &sh, _seis, snl_info->gain[i] );
	/* */
		snl_info->seis[i] = _seis;
	}

	return 0;
}

/**
 * @brief
 *
 * @param snl_info
 * @param path
 * @return int
 */
int seisdata_load_ms( SNL_INFO *snl_info, const char *path )
{
	return -1;
}

/**
 * @brief
 *
 * @param snl_info
 * @param path
 * @return int
 */
int seisdata_load_tank( SNL_INFO *snl_info, const char *path )
{
	return -1;
}