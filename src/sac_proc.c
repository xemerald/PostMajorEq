/*
 *
 */

/* */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <float.h>
/* */
#include <sachead.h>
#include <sac_proc.h>

/* Internal Function Prototypes */
static int    read_sac_header( FILE *, struct SAChead * );
static void   swap_order_4byte( void * );
static void   preprocess_sac_data( float *, const int, const float );
static double fetch_sac_time( const struct SAChead * );
static char  *trim_sac_string( char *, const int );

/*
 *
 */
int sac_proc_station_data_extract(
	const char *sta, const char *net, const char *loc, const char *chan[], const char *path,
	const float gain[], float *seis[], int *npts, float *delta, double *starttime
) {
/* */
	char   filename[512] = { 0 };
	struct SAChead sh;

	int    i;
	float *_seis = NULL;

/* Just a initialization */
	*npts      = -1;
	*delta     = -1.0;
	*starttime = -1.0;
	for ( i = 0; i < 3; i++ )
		seis[i] = NULL;

/* Open all the three channels' SAC files */
	for ( i = 0; i < 3; i++ ) {
	/* Opening the SAC files */
		sprintf(filename, SAC_FILE_NAME_FORMAT, path, sta, chan[i], net, loc);
		if ( sac_proc_sac_load( filename, &sh, &_seis ) < 0 )
			return -1;

	/* Check the consistency of npts */
		if ( *npts < 0 ) {
			*npts = (int)sh.npts;
		}
		else if ( *npts != (int)sh.npts ) {
			fprintf(stderr, "WARNING! There is a different npts within the SAC files of %s.%s.%s.%s.\n", sta, chan[i], net, loc);
			if ( (int)sh.npts < *npts )
				*npts = (int)sh.npts;
		}
	/* */
		if ( *delta < 0.0 ) {
			*delta = sh.delta;
		}
		else if ( fabs(*delta - sh.delta) > FLT_EPSILON ) {
			fprintf(stderr, "ERROR! There is a different delta within the SAC files of %s.%s.%s.%s.\n", sta, chan[i], net, loc);
			free(_seis);
			return -2;
		}
	/* */
		if ( *starttime < 0.0 ) {
			*starttime = fetch_sac_time( &sh );
		}
		else if ( fabs(*starttime - fetch_sac_time( &sh )) > FLT_EPSILON ) {
			fprintf(stderr, "ERROR! There is a different start time within the SAC files of %s.%s.%s.%s.\n", sta, chan[i], net, loc);
			free(_seis);
			return -2;
		}
	/* */
		preprocess_sac_data( _seis, *npts, gain[i] );
	/* */
		seis[i] = _seis;
	}

	return 0;
}

/*
 *
 */
int sac_proc_sac_load( const char *filename, struct SAChead *sh, float **seis )
{
	FILE  *fd;
	float *_seis = NULL;
	int    i, swap_flag, tmp;
	int    result;

/* */
	if ( (fd = fopen(filename, "rb")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening %s\n", filename);
		return -1;
	}
	if ( (swap_flag = read_sac_header(fd, sh)) < 0 ) {
		fclose(fd);
		return -1;
	}
/* Read the sac data into a buffer */
	result = tmp = sh->npts * sizeof(float);
	if ( (_seis = (float *)malloc((size_t)tmp)) == (float *)NULL ) {
		fprintf(stderr, "Out of memory for %d float samples\n", sh->npts);
		fclose(fd);
		return -3;
	}
	if ( (tmp = (int)fread(_seis, sizeof(float), sh->npts, fd)) != sh->npts ) {
		fprintf(stderr, "Error reading SAC data: %s\n", strerror(errno));
		fclose(fd);
		free(_seis);
		return -1;
	}
/* */
	if ( swap_flag == 1 )
		for ( i = 0; i < sh->npts; i++ )
			swap_order_4byte( &(_seis[i]) );

/* */
	*seis   = _seis;
	result += sizeof(struct SAChead);
	fclose(fd);

	return result;
}

/*
 *
 */
struct SAChead *sac_proc_scnl_modify(
	struct SAChead *sh, const char *n_sta, const char *n_chan, const char *n_net, const char *n_loc
) {
	int i;

/* Station name */
	if ( n_sta != NULL ) {
		strncpy(sh->kstnm, n_sta, strlen(n_sta));
		for ( i = (int)strlen(n_sta); i < K_LEN; i++ )
			sh->kstnm[i] = ' ';
	}
/* Channel code */
	if ( n_chan != NULL ) {
		strncpy(sh->kcmpnm, n_chan, strlen(n_chan));
		for ( i = (int)strlen(n_chan); i < K_LEN; i++ )
			sh->kcmpnm[i] = ' ';
	}
/* Network code */
	if ( n_net != NULL ) {
		strncpy(sh->knetwk, n_net, strlen(n_net));
		for ( i = (int)strlen(n_net); i < K_LEN; i++ )
			sh->knetwk[i] = ' ';
	}
/* Location code */
	if ( n_loc != NULL ) {
		strncpy(sh->khole, n_loc, strlen(n_loc));
		for ( i = (int)strlen(n_loc); i < K_LEN; i++ )
			sh->khole[i] = ' ';
	}

	return sh;
}

/*
 *
 */
const char *sac_proc_scnl_print( struct SAChead *sh )
{
	static char result[64] = { 0 };

	char sta[K_LEN]  = { 0 };
	char chan[K_LEN] = { 0 };
	char net[K_LEN]  = { 0 };
	char loc[K_LEN]  = { 0 };

/* */
	strncpy(sta, sh->kstnm, K_LEN);
	trim_sac_string( sta, K_LEN );
/* */
	strncpy(chan, sh->kcmpnm, K_LEN);
	trim_sac_string( chan, K_LEN );
/* */
	strncpy(net, sh->knetwk, K_LEN);
	trim_sac_string( net, K_LEN );
/* */
	strncpy(loc, sh->khole, K_LEN);
	trim_sac_string( loc, K_LEN );

	sprintf(result, "%s.%s.%s.%s", sta, chan, net, loc);

	return result;
}

/*
 *
 */
double sac_proc_reftime_fetch( struct SAChead *sh )
{
	return fetch_sac_time( sh );
}

/*
 *
 */
void sac_proc_data_preprocess( float *seis, const int npts, const float gain_fac )
{
	preprocess_sac_data( seis, npts, gain_fac );
	return;
}

/*
 * read_sac_header: read the header portion of a SAC file into memory.
 *  arguments: file pointer: pointer to an open file from which to read
 *             filename: pathname of open file, for logging.
 * returns: 0 on success
 *          1 on success and if byte swapping is needed
 *         -1 on error reading file
 *     The file is left open in all cases.
 */
static int read_sac_header( FILE *fp, struct SAChead *psh )
{
	int i;
	int filesize;
	int result = 0;
	struct SAChead2 *psh2;

/* obtain file size */
	fseek(fp, 0, SEEK_END);
	filesize = ftell(fp);
	rewind(fp);

/* */
	psh2 = (struct SAChead2 *)psh;
	if ( fread(psh, sizeof(struct SAChead2), 1, fp) != 1 ) {
		fprintf(stderr, "Error reading SAC file: %s!\n", strerror(errno));
		return -1;
	}

/* */
	if ( filesize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
		result = 1;
		fprintf(stderr, "WARNING: Swapping is needed! (filesize %d, psh.npts %d)\n", filesize, psh->npts);
		for ( i = 0; i < NUM_FLOAT; i++ )
			swap_order_4byte( &(psh2->SACfloat[i]) );
		for ( i = 0; i < MAXINT; i++ )
			swap_order_4byte( &(psh2->SACint[i]) );
		if ( filesize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
			fprintf(stderr, "ERROR: Swapping is needed again! (filesize %d, psh.npts %d)\n", filesize, psh->npts);
			result = -1;
		}
	}

	return result;
}

/*
 * Do byte swapping on the given 4-byte integer or float.
 */
static void swap_order_4byte( void *data )
{
	uint8_t temp;

	union {
		uint8_t c[4];
	} _data;

	memcpy( &_data, data, sizeof(uint32_t) );
	temp       = _data.c[0];
	_data.c[0] = _data.c[3];
	_data.c[3] = temp;
	temp       = _data.c[1];
	_data.c[1] = _data.c[2];
	_data.c[2] = temp;
	memcpy( data, &_data, sizeof(uint32_t) );

	return;
}

/*
 *
 */
static double fetch_sac_time( const struct SAChead *sh )
{
	struct tm tms;
    double result;

    tms.tm_year  = sh->nzyear - 1900;
    tms.tm_mon   = 0;           /* Force the month to January */
    tms.tm_mday  = sh->nzjday;  /* tm_mday is 1 - 31; nzjday is 1 - 366 */
    tms.tm_hour  = sh->nzhour;
    tms.tm_min   = sh->nzmin;
    tms.tm_sec   = sh->nzsec;
    tms.tm_isdst = 0;
    result       = (double)timegm(&tms);
	result      += sh->nzmsec / 1000.0;

    return result;
}

/*
 *
 */
static char *trim_sac_string( char *input, const int size )
{
	char *p = input + size - 1;

/* */
	if ( input != NULL && *input != '\0' ) {
		while ( p >= input && isspace(*p) ) {
			*p-- = '\0';
		}
	}

	return input;
}

/*
 *
 */
static void preprocess_sac_data( float *input, const int npts, const float gain )
{
	int   i, i_head;
	int   gap_count = 0;
	float mean = 0.0;

/* */
	i_head = (int)(npts * 0.1);
	i_head = i_head ? i_head : npts;
	for ( i = 0; i < npts; i++ ) {
		if ( input[i] == SACUNDEF ) {
			input[i] = 0.0;
			if ( i < i_head )
				gap_count++;
			continue;
		}
		input[i] *= gain;
		if ( i < i_head )
			mean += input[i];
	}
	mean /= (float)(i_head - gap_count);
/* */
	for ( i = 0; i < npts; i++ )
		if ( input[i] != 0.0 )
			input[i] -= mean;
		else
			input[i] = mean;

	return;
}
