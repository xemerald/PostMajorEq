/**
 * @file seisdata_load.c
 * @author Benjamin Yang @ National Taiwan University (b98204032@gmail.com)
 * @brief
 * @version 1.0.0
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
#include <fcntl.h>
/* */
#include <sys/mman.h>
#include <sys/stat.h>
/* */
#include <postmajor.h>
#include <sac.h>
#include <libmseed.h>

/* */
#define MAX_FILE_NAME          512
#define SAC_FILE_NAME_FORMAT  "%s/%s.%s.%s.%s"
/* */
static float *subs_gap2nan( float [], const int, const float );
static float *apply_gain2data( float [], const int, const float );
static float *dmean_data( float [], const int, const double, int * );
/* */
static void *LoadContextSAC  = NULL;
static void *LoadContextMS   = NULL;
static void *LoadContextTANK = NULL;

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
	int    gap   = 0;

/* Just a initialization */
	snl_info->npts      = -1;
	snl_info->delta     = -1.0;
	snl_info->starttime = -1.0;
	for ( register int i = 0; i < NUM_CHANNEL_SNL; i++ )
		snl_info->seis[i] = NULL;

/* Open all the three channels' SAC files */
	for ( register int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
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
		subs_gap2nan( _seis, snl_info->npts, SACUNDEF );
		apply_gain2data( _seis, snl_info->npts, snl_info->gain[i] );
		dmean_data( _seis, snl_info->npts, 1.0 / snl_info->delta, &gap );
		if ( gap ) {
			fprintf(
				stderr, "Found %d gaps within total %d samples in %s, filled with mean value!\n",
				gap, snl_info->npts, sac_scnl_print( &sh )
			);
		}
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
	MS3TraceID *tid[NUM_CHANNEL_SNL] = { NULL };
	int         offset   = 0;
	int         npts     = 0;
	int         seis_idx = 0;
	float      *_seis    = NULL;
	nstime_t    earliest = 0;
	nstime_t    latest   = 0;
	nstime_t    lastend  = 0;
	double      samprate = -1.0;
	uint8_t     samplesize;
	char        sampletype;
	char        sid[LM_SIDLEN] = { 0 };

/* Only mapping the trace list at the first time */
	if ( !LoadContextMS ) {
	/* */
		LoadContextMS = mstl3_init(NULL);
	/* Read all miniSEED from the path, accumulate in MS3TraceList */
		fprintf(stderr, "Mapping the miniSEED file %s into memory...\n", path);
		if ( ms3_readtracelist((MS3TraceList **)&LoadContextMS, path, NULL, 0, MSF_VALIDATECRC | MSF_RECORDLIST, 0) != MS_NOERROR ) {
			fprintf(
				stderr, "ERROR! Cannot read miniSEED from file: %s\n", path
			);
			return -2;
		}
	}

/* */
	for ( register int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
	/* */
		snl_info->seis[i] = NULL;
	/* */
		ms_nslc2sid(sid, LM_SIDLEN, 0, snl_info->net, snl_info->sta, strcmp("--", snl_info->loc) ? snl_info->loc : NULL, snl_info->chan[i]);
		if ( !(tid[i] = mstl3_findID((MS3TraceList *)LoadContextMS, sid, 0, NULL)) ) {
			fprintf(stderr, "ERROR! Cannot find the SID: %s in the miniSEED file: %s\n", sid, path);
			return -1;
		}
	/* */
		if (
			(earliest && tid[i]->latest < earliest) ||
			(latest && tid[i]->earliest > latest)
		) {
			fprintf(stderr, "ERROR! There is an out of time range trace within the miniSEED files of SID: %s\n", sid);
			return -2;
		}
	/* */
		if ( !earliest || earliest > tid[i]->earliest )
			earliest = tid[i]->earliest;
	/* */
		if ( !latest || latest < tid[i]->latest )
			latest = tid[i]->latest;
	/* Check the consistency of delta */
		if ( samprate < 0.0 ) {
			samprate = tid[i]->first->samprate;
		}
		else if ( fabs(samprate - tid[i]->first->samprate) > FLT_EPSILON ) {
			fprintf(stderr, "ERROR! There is a different sampleing rate within the miniSEED files of SID: %s\n", sid);
			return -2;
		}
	}
/* Just derive the maximum number of samples we need here */
	npts = (latest - earliest) * 1.0e-9 * samprate + 1;
/* Again, go thru all the traces */
	for ( register int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
	/* Create the buffer space for storaging the seismic data, and then fill with the gap value, NAN */
		if ( !(_seis = (float *)calloc(npts, sizeof(float))) ) {
			fprintf(stderr, "ERROR! Out of memory for %d float samples\n", npts);
			return -2;
		}
		for ( register int j = 0; j < npts; j++ )
			_seis[j] = NAN;
	/* Start to read in the trace segments and copy the data into buffer */
		lastend = earliest;
		seis_idx = 0;
		for ( MS3TraceSeg *seg = tid[i]->first; seg; seg = seg->next ) {
		/* Check the data sample size & type and then unpack it */
			ms_encoding_sizetype(seg->recordlist->first->msr->encoding, &samplesize, &sampletype);
			mstl3_unpack_recordlist(tid[i], seg, NULL, 0, 0);
		/* Calculate the gap and offset the index */
			if ( (offset = (seg->starttime - lastend) * 1.0e-9 * seg->samprate - 1) > 0 )
				seis_idx += offset;
		/* */
			switch ( sampletype ) {
			case 'i':
				for ( register int j = 0; j < seg->numsamples; j++ )
					_seis[seis_idx++] = *((int32_t *)seg->datasamples + j);
				break;
			case 'f':
				for ( register int j = 0; j < seg->numsamples; j++ )
					_seis[seis_idx++] = *((float *)seg->datasamples + j);
				break;
			case 'd':
				for ( register int j = 0; j < seg->numsamples; j++ )
					_seis[seis_idx++] = *((double *)seg->datasamples + j);
				break;
			default:
				break;
			}
		/* Save the last end time for next loop usage */
			lastend = seg->endtime;
		}
	/* Preprocess the seismic data */
		apply_gain2data( _seis, npts, snl_info->gain[i] );
		dmean_data( _seis, npts, samprate, &seis_idx );
		if ( seis_idx ) {
			fprintf(
				stderr, "Found %d gaps within total %d samples in SID: %s, filled with mean value!\n",
				seis_idx, npts, tid[i]->sid
			);
		}
	/* Keep the buffer pointer */
		snl_info->seis[i] = _seis;
	}
/* */
	snl_info->npts      = npts;
	snl_info->delta     = 1.0 / samprate;
	snl_info->starttime = earliest * 1.0e-9;

	return 0;
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

/**
 * @brief
 *
 */
void seisdata_release_sac( void )
{
	if ( LoadContextSAC )
		free(LoadContextSAC);

	return;
}

/**
 * @brief
 *
 */
void seisdata_release_ms( void )
{
	if ( LoadContextMS )
		mstl3_free((MS3TraceList **)&LoadContextMS, 0);

	return;
}

/**
 * @brief
 *
 */
void seisdata_release_tank( void )
{
	if ( LoadContextTANK )
		free(LoadContextTANK);

	return;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param gap_value
 * @return float*
 */
static float *subs_gap2nan( float input[], const int npts, const float gap_value )
{
	for ( register int i = 0; i < npts; i++ )
		if ( fabs(input[i] - gap_value) < FLT_EPSILON )
			input[i] = NAN;

	return input;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param gain
 * @return float*
 */
static float *apply_gain2data( float input[], const int npts, const float gain )
{
/* */
	for ( register int i = 0; i < npts; i++ ) {
		if ( !isnan(input[i]) )
			input[i] *= gain;
	}

	return input;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param samprate
 * @param gap_count
 * @return float*
 */
static float *dmean_data( float input[], const int npts, const double samprate, int *gap_count )
{
	register int   mean_count = 0;
	register int   _gap_count = 0;
	register float mean       = 0.0;
	register int   i_head;

/* First, calculate the mean value from the head part of data */
	i_head = (int)(npts * 0.1);
	i_head = i_head >= (int)samprate ? i_head : npts;
	for ( register int i = 0; i < i_head; i++ ) {
		if ( !isnan(input[i]) ) {
			mean += input[i];
			mean_count++;
		}
	}
	mean /= mean_count;
/* */
	for ( register int i = 0; i < npts; i++ ) {
		if ( !isnan(input[i]) ) {
			input[i] -= mean;
		}
		else {
			input[i] = mean;
			_gap_count++;
		}
	}
/* */
	if ( gap_count )
		*gap_count = _gap_count;

	return input;
}
