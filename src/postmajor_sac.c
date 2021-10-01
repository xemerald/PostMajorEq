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
#include <time.h>
#include <math.h>
#include <errno.h>
/* */
#include <sachead.h>
#include <iirfilter.h>
#include <picker_wu.h>
#include <sac_proc.h>

/* */
#define  MAX_STR_SIZE    512

#define  PI  3.141592653589793238462643383279f
#define  PI2 6.283185307179586476925286766559f

/* Internal Function Prototypes */
static int    parse_stalist_line( const char *, char *, float *, float *, float *, float [] );
static int    parse_eqinfo_file( const char *, float *, float *, float *, float * );
static float *integral_waveform( float *, const int, const double );
static float *highpass_filter( float *, const int, const double, const int );
static float  calc_tau_c( const float *, const int, const float, const int );
static float  calc_peak_value( const float *, const int, const float, const int );
static double coor2distf( const double, const double, const double, const double );

/*
 *
 */
int main( int argc, char **argv )
{
	int    i, j;
	FILE  *fd;
	char   list_line[MAX_STR_SIZE] = { 0 };
	char   sta_c[8] = { 0 };
	float *seis[3];       /* input trace buffer */
/* */
	int   eq_flag = 0;
	float otime, elat, elon, edep;
/* */
	int   flag = 0;
	int   npts;
	float delta = -1.0;
	float gain[3];
	float starttime;
	float lat, lon, elev;

/* Derived from waveform */
	int    nparrival;
	float  epic_dist;
	float  pga, pgv, pgd;
	float  pa3, pv3, pd3, tc;
	float  pga_time, pgv_time, pgd_time;
	float  pd35_time, pga80_time, pga4_time;
	float  pga_leadtime, pgv_leadtime;
	double snr;

/* Check command line arguments */
	if ( argc != 4 ) {
		fprintf(stderr, "Usage: %s <Eq. Info> <Station List> <SAC Files Path>\n", argv[0]);
		exit(0);
	}
/* */
	if ( !parse_eqinfo_file( argv[1], &otime, &elat, &elon, &edep ) ) {
		eq_flag = 1;
	}
/* */
	if ( (fd = fopen(argv[2], "r")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening station list %s\n", argv[2]);
		exit(1);
	}
/* */
	fprintf(stdout, "# Station  PGA  PGV  PGD  PA3  PV3  PD3  TauC3  PGA_Leading  PGV_Leading  Epc_Dist  S/N_Ratio\n");

/* */
	while ( fgets(list_line, sizeof(list_line) - 1, fd) != NULL ) {
	/* */
		if ( parse_stalist_line( list_line, sta_c, &lat, &lon, &elev, gain ) )
			continue;
	/* */
		if ( sac_proc_station_data_extract( sta_c, argv[3], gain, seis, &npts, &delta, &starttime ) ) {
			for ( i = 0; i < 3; i++ ) {
				if ( seis[i] != NULL ) {
					free(seis[i]);
					seis[i] = NULL;
				}
			}
			continue;
		}

		fprintf(
			stderr, "Processing data of %s (start at %lf, npts %d, delta %.2lf)... \n",
			sta_c, starttime, npts, delta
		);
	/* */
		flag = 0;
		if ( (nparrival = pickwu_p_arrival_pick( seis[0], npts, delta, 2, 0 )) ) {
			if ( pickwu_p_trigger_check( seis[0], npts, delta, nparrival ) ) {
				if ( pickwu_p_arrival_quality_calc( seis[0], npts, delta, nparrival, &snr ) < 4 ) {
					flag = 1;
				}
			}
		}
		if ( flag == 0 ) {
			fprintf(
				stderr, "Can't find valid P arrival time (SNR: %lf), just skip the station %s.\n",
				snr, sta_c
			);
			continue;
		}

	/*
	 *
	 */
		pga        = 0.0;
		pga_time   = 0.0;
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
		pa3 = calc_peak_value( seis[0] + nparrival, npts, delta, 3 );

	/* Transform the acceleration sample to velocity sample */
		for ( i = 0; i < 3; i++ )
			integral_waveform( seis[i], npts, delta );
	/* */
		pgv      = 0.0;
		pgv_time = 0.0;
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
		pv3 = calc_peak_value( seis[0] + nparrival, npts, delta, 3 );

	/* Transform the velocity sample to displacement sample */
		for ( i = 0; i < 3; i++ )
			integral_waveform( seis[i], npts, delta );
	/* */
		pgd       = 0.0;
		pgd_time  = 0.0;
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
		pd3 = calc_peak_value( seis[0] + nparrival, npts, delta, 3 );
		tc  = calc_tau_c( seis[0] + nparrival, npts, delta, 3 );

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
	/* */
		if ( eq_flag )
			epic_dist = coor2distf( elon, elat, lon, lat );
		else
			epic_dist = 0.0;

		fprintf(
			stdout, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
			sta_c, pga, pgv, pgd, pa3, pv3, pd3, tc, pga_leadtime, pgv_leadtime, epic_dist, snr
		);
	}

	fclose(fd);

	return 0;
}

/*
 *
 */
static int parse_stalist_line( const char *line, char *sta, float *lat, float *lon, float *elev, float gain[] )
{
	int i;

/* */
	if ( strlen(line) ) {
	/* */
		for ( i = 0; i < MAX_STR_SIZE; i++ ) {
			if ( line[i] == '#' || line[i] == '\n' ) {
				break;
			}
			else if ( line[i] == '\t' || line[i] == ' ' ) {
				continue;
			}
			else {
				if ( sscanf(line, "%s %f %f %f %f %f %f\n", sta, lat, lon, elev, &gain[0], &gain[1], &gain[2]) == 7 )
					return 0;
			}
		/* */
			break;
		}
	}

	return -1;
}

/*
 *
 */
static int parse_eqinfo_file( const char *path, float *otime, float *epc_lat, float *epc_lon, float *dep )
{
	int   i;
	FILE *fd;
	char  line[MAX_STR_SIZE] = { 0 };
	float year, mon, day, hour, min, sec;
	struct tm _otime;

/* */
	if ( (fd = fopen(path, "r")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening Eq. information %s\n", path);
		return -2;
	}
/* */
	while ( fgets(line, sizeof(line) - 1, fd) != NULL ) {
		if ( strlen(line) ) {
		/* */
			for ( i = 0; i < MAX_STR_SIZE; i++ ) {
				if ( line[i] == '#' || line[i] == '\n' ) {
					break;
				}
				else if ( line[i] == '\t' || line[i] == ' ' ) {
					continue;
				}
				else {
					if (
						sscanf(
							line, "%f %f %f %f %f %f %f %f %f\n",
							&year, &mon, &day, &hour, &min, &sec, epc_lat, epc_lon, dep
						) == 9
					) {
					/* */
						_otime.tm_year  = (int)year - 1900;
						_otime.tm_mon   = (int)mon - 1;
						_otime.tm_mday  = (int)day;
						_otime.tm_hour  = (int)hour;
						_otime.tm_min   = (int)min;
						_otime.tm_sec   = 0;
						_otime.tm_isdst = 0;

						*dep   = -(*dep);
						*otime = (double)timegm(&_otime) + sec;

						fclose(fd);
						return 0;
					}
				}
			/* */
				break;
			}
		}
	}
/* */
	fclose(fd);
	return -1;
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
static float calc_tau_c( const float *input, const int npts, const float delta, const int sec )
{
	int   i;
	float result  = 0.0;
	float sum_dis = 0.0;
	float sum_vel = 0.0;

	const int i_end = (int)(sec / delta);

/* */
	if ( i_end > npts )
		return 0.0;

	for ( i = 1; i < i_end; i++ ) {
		result   = (input[i] - input[i-1]) / delta;
		sum_dis += input[i] * input[i];
		sum_vel += result * result;
	}

	result = PI2 * sqrt(sum_dis / sum_vel);

	return result > 10.0 ? 10.0 : result;
}

/*
 *
 */
static float calc_peak_value( const float *input, const int npts, const float delta, const int sec )
{
	int   i;
	float result = 0.0;

	const int i_end = (int)(sec / delta);
/* */
	if ( i_end > npts )
		return 0.0;

	for ( i = 1; i < i_end; i++ ) {
	/* */
		if ( fabs(input[i]) > result )
			result = fabs(input[i]);
	}

	return result;
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
