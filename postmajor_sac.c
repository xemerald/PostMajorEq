/*
 * Standalone program to read SAC data files and compute
 * the peak values include acceleration, velocity & displacement.
 *
 * It can also derive the warning time for each station.
 *
 * Benjamin Yang; Feb 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <errno.h>
/* */
#include <sachead.h>
#include <iirfilter.h>
#include <picker_wu.h>


#define EV_LON 120.54f
#define EV_LAT 22.92f
#define EV_DEP 14.6f

#define PI  3.141592653589793238462643383279f
#define PI2 6.283185307179586476925286766559f

/* Internal Function Prototypes */
static char *trim_sac_string( char *, const int );
static float fetch_sac_time( const struct SAChead * );
static int read_sac_header( FILE *, struct SAChead * );
static void swap_order_4byte( void * );
static float *integral_waveform( float *, const int, const double );
static float *highpass_filter( float *, const int, const double, const int );
static float cal_tau_c( float *, const int, const float, const int );
static double coor2distf( const double, const double, const double, const double );

/* */
typedef struct {
/* Fixed information */
	char sta[8];
	char net[8];
	char loc[8];
	float latitude;
	float longitude;

/* Computed information */
	float pga;
	float pgv;
	float pgd;
	float pd;
	float tc;
	float pga_leadtime;
	float pgv_leadtime;
	float p_arrival;
	float s_arrival;
	float p_weight;
	float s_weight;
} STAINFO;

/*
 *
 */
int main( int argc, char **argv )
{
	struct SAChead sh;

	FILE  *fp, *fp_list;
	char   sta_c[K_LEN + 1] = { 0 };
	char   filename[3][128] = { 0 };
	float *seis[3];       /* input trace buffer */

	int   i, j;
	int   _flag = 0;
	int   npts, nparrival;
	int   samprate = 0;
	float delta = -1.0;
	float float_tmp = 0.0;
	float lat, lon;
	float gain_factor;

/* Final output data */
	float epic_dist;
	float pga, pgv, pgd;
	float pa3, pv3, pd3, tc;
	float pga_time, pgv_time, pgd_time;
	float pd35_time, pga80_time, pga4_time, pga2_time;
	float pga_leadtime, pgv_leadtime;

	double snr;

/* Check command line arguments */
	if ( argc != 2 ) {
		fprintf(stderr, "Usage: %s <Station List>\n", argv[0]);
		//fprintf(stderr, "Usage: %s <Z Component File> <N Component File> <E Component File>\n", argv[0]);
		exit(0);
	}

	if ( (fp_list = fopen(argv[1], "r")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening %s\n", argv[1]);
		exit(1);
	}


	while ( 1 ) {
		if ( fscanf(fp_list, "%s %f %f %*f %*f %*f %*f\n", sta_c, &lat, &lon) != 3 ) {
			return 0;
		}

		sprintf(filename[0], "%s.HLZ.TW.--", sta_c);
		sprintf(filename[1], "%s.HLN.TW.--", sta_c);
		sprintf(filename[2], "%s.HLE.TW.--", sta_c);
		if ( access(filename[0], F_OK) != 0 )
			continue;
	/* Open all the three SAC files */
		//printf("Opening %s\n", filename[0]);
		for ( i = 0; i < 3; i++ ) {
		/* Opening the SAC files */
			if ( (fp = fopen(filename[i], "rb")) == (FILE *)NULL ) {
				fprintf(stderr, "Error opening %s\n", argv[i+1]);
				continue;
				exit(1);
			}
			if ( (_flag = read_sac_header(fp, &sh)) < 0 ) {
				fclose(fp);
				exit(1);
			}
		/* */
			//strncpy(sta_c, sh.kstnm, K_LEN);
			//trim_sac_string( sta_c, K_LEN );
		/* */
			npts = (int)sh.npts;
			j    = npts * sizeof(float);
			if ( (seis[i] = (float *)malloc((size_t)j)) == (float *)NULL ) {
				fprintf(stderr, "Out of memory for %d float samples\n", npts);
				fclose(fp);
				exit(1);
			}
		/* Read the sac data into a buffer */
			if ( (j = (int)fread(seis[i], sizeof(float), sh.npts, fp)) != npts ) {
				fprintf(stderr, "Error reading SAC data: %s\n", strerror(errno));
				fclose(fp);
				exit(1);
			}
			fclose(fp);
		/* */
			if ( _flag == 1 )
				for ( j = 0; j < npts; j++ )
					swap_order_4byte( &(seis[i][j]) );

		/* */
			if ( delta < 0.0 ) {
				delta      = sh.delta;
				samprate   = (int)(1.0 / delta);
			}
			else if ( fabs(delta - sh.delta) > FLT_EPSILON ) {
				fprintf(stderr, "There is a different delta within the SAC files, exiting!\n");
				exit(1);
			}
		/* */
			switch ( sh.kcmpnm[2] ) {
			case 'Z':
				gain_factor = 0.059814;
				break;
			case 'N':
				gain_factor = 0.059814;
				break;
			case 'E':
				gain_factor = 0.059814;
				break;
			default:
				//fprintf(stderr, "The SAC header do not define the component or not clear at all, we can't define the direction!\n");
				gain_factor = 0.059814;
				break;
			}
		/* */
			_flag     = 0;
			float_tmp = 0.0;
			for ( j = 0; j < npts; j++ ) {
				if ( seis[i][j] == SACUNDEF ) {
					seis[i][j] = 0.0;
					_flag++;
					continue;
				}
				seis[i][j] *= gain_factor;
				float_tmp  += seis[i][j];
			}
			float_tmp /= (float)(npts - _flag);
		/* */
			for ( j = 0; j < npts; j++ )
				if ( seis[i][j] != 0.0 )
					seis[i][j] -= float_tmp;
		}

	/* */
		_flag = 0;
		if ( (nparrival = pickwu_p_arrival_pick( seis[0], npts, delta, 2, 0 )) ) {
			if ( pickwu_p_trigger_check( seis[0], npts, delta, nparrival ) ) {
				if ( pickwu_p_arrival_quality_calc( seis[0], npts, delta, nparrival, &snr ) <= 3 ) {
					_flag = 1;
				}
			}
		}

		if ( _flag == 0 ) {
			fprintf(stderr, "Can't find the P arrival time, just skip this station.\n");
			continue;
		}


	/*
	 *
	 */
		pga        = 0.0;
		pga4_time  = npts * delta;
		pga80_time = npts * delta;
		for ( i = 0; i < 3; i++ ) {
			for ( j = nparrival; j < npts; j++ ) {
			/* */
				if ( fabs(seis[i][j]) > 4.0 ) {
				/* */
					if ( (j * delta) < pga4_time )
						pga4_time = j * delta;
				/* */
					if ( fabs(seis[i][j]) > 80.0 ) {
						if ( (j * delta) < pga80_time )
							pga80_time = j * delta;
					}
				}
			/* */
				if ( fabs(seis[i][j]) > pga ) {
					pga      = fabs(seis[i][j]);
					pga_time = j * delta;
				}
			}
		}
	/* */
		pa3 = 0.0;
		for ( i = nparrival; i < nparrival + samprate * 3; i++ ) {
		/* */
			if ( fabs(seis[0][i]) > pa3 )
				pa3 = fabs(seis[0][i]);
		}

	/* Transform the acceleration sample to velocity sample */
		for ( i = 0; i < 3; i++ )
			integral_waveform( seis[i], npts, delta );
	/* */
		pgv = 0.0;
		for ( i = 0; i < 3; i++ ) {
			for ( j = nparrival; j < npts; j++ ) {
			/* */
				if ( fabs(seis[i][j]) > pgv ) {
					pgv      = fabs(seis[i][j]);
					pgv_time = j * delta;
				}
			}
		}
	/* */
		pv3 = 0.0;
		for ( i = nparrival; i < nparrival + samprate * 3; i++ ) {
		/* */
			if ( fabs(seis[0][i]) > pv3 )
				pv3 = fabs(seis[0][i]);
		}

	/* Transform the velocity sample to displacement sample */
		for ( i = 0; i < 3; i++ )
			integral_waveform( seis[i], npts, delta );
	/* */
		pgd       = 0.0;
		pd35_time = npts * delta;
		for ( i = 0; i < 3; i++ ) {
			for ( j = nparrival; j < npts; j++ ) {
			/* */
				if ( fabs(seis[i][j]) > 0.35 ) {
					if ( (j * delta) < pd35_time )
						pd35_time = j * delta;
				}
			/* */
				if ( fabs(seis[i][j]) > pgd ) {
					pgd      = fabs(seis[i][j]);
					pgd_time = j * delta;
				}
			}
		}

	/* Computation of Pd & Tau-c at 3 seconds */
		pd3 = 0.0;
		for ( i = nparrival; i < nparrival + samprate * 3; i++ ) {
		/* */
			if ( fabs(seis[0][i]) > pd3 )
				pd3 = fabs(seis[0][i]);
		}
		tc = cal_tau_c( seis[0] + nparrival, npts, delta, 3 );

	/*
	 *
	 */
	/* */
		if ( pd35_time <= 0.0 && pga80_time <= 0.0 ) {
			pga_leadtime = pgv_leadtime = 0.0;
		}
		else {
		/* */
			if ( (pga_time - pd35_time) <= (pga_time - pga80_time) )
				pga_leadtime = pga_time - pga80_time;
			else
				pga_leadtime = pga_time - pd35_time;
		/* */
			if ( (pgv_time - pd35_time) <= (pgv_time - pga80_time) )
				pgv_leadtime = pgv_time - pga80_time;
			else
				pgv_leadtime = pgv_time - pd35_time;
		/* */
			if ( pga_leadtime < 0.0 )
				pga_leadtime = 0.0;
			if ( pgv_leadtime < 0.0 )
				pgv_leadtime = 0.0;
		}

	/* */
		for ( i = 0; i < 3; i++ )
			free(seis[i]);

		epic_dist = coor2distf( EV_LON, EV_LAT, lon, lat );

		fprintf(
			stdout, "%s %lf %lf %lf %lf %lf %lf %lf %lf\n",
			sta_c, pga, pgv, pgd, pd3, tc, pga_leadtime, pgv_leadtime, epic_dist
		);
	}

	fclose(fp_list);

	return 0;
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
static int read_sac_header(FILE *fp, struct SAChead *psh)
{
	int i;
	struct SAChead2 *psh2;
	int fileSize;
	int ret = 0;

/* obtain file size */
	fseek(fp, 0, SEEK_END);
	fileSize = ftell(fp);
	rewind(fp);

	psh2 = (struct SAChead2 *)psh;

	if ( fread(psh, sizeof(struct SAChead2), 1, fp) != 1 ) {
		fprintf(stderr, "read_sac_header: error reading SAC file: %s\n", strerror(errno));
		return -1;
	}

/* mtheo 2007/10/19
* Guessing if byte swapping is needed
* fileSize should be equal to sizeof(struct SAChead) + (psh->npts * sizeof(float)) */
	if( fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
		ret = 1;
		fprintf(stderr, "WARNING: Swapping is needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
		for ( i=0; i<NUM_FLOAT; i++ )
			swap_order_4byte( &(psh2->SACfloat[i]) );
		for ( i=0; i<MAXINT; i++ )
			swap_order_4byte( &(psh2->SACint[i]) );
		if( fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
			fprintf(stderr, "ERROR: Swapping is needed again! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
			ret = -1;
		}
	} else {
		fprintf(stderr, "Swapping is not needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
	}

	return ret;
}

/* Do byte swapping on the given 4-byte integer or float. */
static void swap_order_4byte( void *data )
{
   uint8_t temp;

   union {
      uint8_t c[4];
   } dat;

   memcpy( &dat, data, sizeof(uint32_t) );
   temp     = dat.c[0];
   dat.c[0] = dat.c[3];
   dat.c[3] = temp;
   temp     = dat.c[1];
   dat.c[1] = dat.c[2];
   dat.c[2] = temp;
   memcpy( data, &dat, sizeof(uint32_t) );
   return;
}

/*
 *
 */
static float fetch_sac_time( const struct SAChead *sh )
{
	struct tm tms;
    double sec;

    tms.tm_year  = sh->nzyear - 1900;
    tms.tm_mon   = 0;    /* Force the month to January */
    tms.tm_mday  = sh->nzjday;  /* tm_mday is 1 - 31; nzjday is 1 - 366 */
    tms.tm_hour  = sh->nzhour;
    tms.tm_min   = sh->nzmin;
    tms.tm_sec   = sh->nzsec;
    tms.tm_isdst = 0;
    sec          = (double)timegm(&tms);

    return (sec + (sh->nzmsec / 1000.0));
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
static float *integral_waveform( float *input, const int npts, const double delta )
{
	int   i;
	float last_seis  = 0.0;
	float this_pseis = 0.0;

	const float half_delta = delta * 0.5;

/* */
	input[0] = 0.0;
	for ( i = 1; i < npts; i++ ) {
		this_pseis = (input[i] + last_seis) * half_delta + input[i-1];
		last_seis  = input[i];
		input[i]   = this_pseis;
	}
/* */
	highpass_filter( input, npts, delta, 0 );

	return input;
}

/*
 *
 */
static float *highpass_filter( float *input, const int npts, const double delta, const int zero_phase )
{
	int        i;
	IIR_FILTER filter;
	IIR_STAGE *stage;

	filter = designfilter( 2, IIR_HIGHPASS_FILTER, IIR_BUTTERWORTH, 0.075, 0.0, delta );
	stage  = (IIR_STAGE *)calloc(filter.nsects, sizeof(IIR_STAGE));
/* First time, forward filtering */
	memset(stage, 0, sizeof(IIR_STAGE) * filter.nsects);
	for ( i = 0; i < npts; i++ )
		input[i] = applyfilter( input[i], &filter, stage );

/* Second time, backward filtering */
	if ( zero_phase ) {
		memset(stage, 0, sizeof(IIR_STAGE) * filter.nsects);
		for ( i = npts - 1; i >= 0; i-- )
			input[i] = applyfilter( input[i], &filter, stage );
	}

	free(stage);

	return input;
}

/*
 *
 */
static float cal_tau_c( float *input, const int npts, const float delta, const int sec )
{
	int   i;
	float _tmp    = 0.0;
	float sum_dis = 0.0;
	float sum_vel = 0.0;

	const int i_end = (int)(sec / delta);

/* */
	if ( i_end > npts )
		return 0.0;

	for ( i = 1; i < i_end; i++ ) {
		_tmp     = (input[i] - input[i-1]) / delta;
		sum_dis += input[i] * input[i];
		sum_vel += _tmp * _tmp;
	}

	_tmp = PI2 * sqrt(sum_dis / sum_vel);

	return _tmp > 10.0 ? 10.0 : _tmp;
}

/*
 * coor2distf() - Transforms the coordinate(latitude & longitude) into distance(unit: km)
 */
static double coor2distf( const double elon, const double elat, const double slon, const double slat )
{
	const double avlat = (elat + slat)*0.5;

	double a = 1.840708 + avlat*(.0015269 + avlat*(-.00034 + avlat*(1.02337e-6)));
	double b = 1.843404 + avlat*(-6.93799e-5 + avlat*(8.79993e-6 + avlat*(-6.47527e-8)));

	a *= (slon - elon) * 60.0;
	b *= (slat - elat) * 60.0;

	return sqrt(a * a + b * b);
}
