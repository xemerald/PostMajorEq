/**
 * @file postmajor.c
 * @author your name (you@domain.com)
 * @brief Standalone program to read seismic data files and compute the peak values include acceleration, velocity & displacement.
 *        It can also derive the warning time for each station.
 * @version 2.0.0
 * @date 2024-04-05
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <errno.h>
/* */
#include <postmajor.h>
#include <seisdata_load.h>
#include <iirfilter.h>
#include <picker_wu.h>

/* Internal Function Prototypes */
static int    proc_argv( int, char * [] );
static void   usage( void );
static void   init_snl_info_params( SNL_INFO * );
static int    parse_stalist_line( SNL_INFO *, const char * );
static int    parse_stalist( SNL_INFO **, const char * );
static int    parse_eqinfo_file( const char *, float *, float *, float *, double * );
static int    pick_pwave_arrival( SNL_INFO *, const double );
static void   proc_acc( SNL_INFO *, const int );
static void   proc_vel( SNL_INFO *, const int );
static void   proc_disp( SNL_INFO *, const int );
static void   proc_leadtime( SNL_INFO * );
static void   integral_waveforms( SNL_INFO *, const _Bool );
static void   differential_waveform( SNL_INFO * );
static float *_integral_waveform( float *, const int, const double, const _Bool );
static float *_differential_waveform( float *, const int, const double );
static float *highpass_filter( float *, const int, const double, const _Bool );
static float  calc_tau_c( const float *, const float *, const int, const float, const int );
static float  calc_peak_value( const float *, const int, const float, const int );
static double coor2distf( const double, const double, const double, const double );
/* */
static _Bool HeaderSwitch      = true;
static _Bool CoordinateSwitch  = false;
static _Bool IgnStaWithoutData = false;
static _Bool IgnStaWithoutPick = false;
static _Bool TwoStageIntegral  = false;
static _Bool VecSumSwitch      = false;
static char *EqInfoFile        = NULL;
static char *StaListFile       = NULL;
static char *SeisDataFile      = NULL;
static int (*LoadSeisdataFunc)( SNL_INFO *, const char * ) = seisdata_load_sac;

/**
 * @brief
 *
 * @param argc
 * @param argv
 * @return int
 */
int main( int argc, char **argv )
{
/* */
	float  elat;
	float  elon;
	float  edep;
	double otime = 0.0;
/* */
	SNL_INFO *snl_infos    = NULL;
	int       totalsnl     = 0;
	int       end_pos      = 0;

/* Check command line arguments */
	if ( proc_argv( argc, argv ) ) {
		usage();
		return -1;
	}

/* */
	if ( parse_eqinfo_file( EqInfoFile, &elat, &elon, &edep, &otime ) < 0 )
		return -1;
/* */
	if ( (totalsnl = parse_stalist( &snl_infos, StaListFile )) <= 0 )
		return -1;

/* */
	for ( int i = 0; i < totalsnl; i++ ) {
	/* */
		snl_infos[i].epic_dist = coor2distf( elon, elat, snl_infos[i].longitude, snl_infos[i].latitude );

	/* */
		if ( LoadSeisdataFunc( &snl_infos[i], SeisDataFile ) < 0 ) {
			init_snl_info_params( &snl_infos[i] );
			continue;
		}

	/* */
		fprintf(
			stderr, "Processing data of %s.%s.%s (start at %lf, npts %d, delta %.2lf)... \n",
			snl_infos[i].sta, snl_infos[i].net, snl_infos[i].loc, snl_infos[i].starttime, snl_infos[i].npts, snl_infos[i].delta
		);

	/* Set the time before origin time 1 sec. as the start point for scaning */
		if ( !( snl_infos[i].pick_flag = pick_pwave_arrival( &snl_infos[i], otime )) ) {
			fprintf(
				stderr, "Can't find valid P arrival (Np: %d, SNR: %lf), skip those time related parameters for SNL %s.%s.%s.\n",
				snl_infos[i].parrival_pos, snl_infos[i].snr, snl_infos[i].sta, snl_infos[i].net, snl_infos[i].loc
			);
		}
	/* */
		if ( (end_pos = snl_infos[i].parrival_pos + (int)(EV_DURATION / snl_infos[i].delta) + 1) > snl_infos[i].npts )
			end_pos = snl_infos[i].npts;
	/* */
		snl_infos[i].sum_vel = calloc(snl_infos[i].npts, sizeof(float));
		snl_infos[i].sum_dis = calloc(snl_infos[i].npts, sizeof(float));
	/* First of all, process the raw acceleration sample */
		proc_acc( &snl_infos[i], end_pos );
	/* Fork the process depends on the two stage integral switch */
		if ( TwoStageIntegral ) {
		/* Transform the acceleration sample to velocity sample */
			integral_waveforms( &snl_infos[i], true );
		/* */
			proc_vel( &snl_infos[i], end_pos );
		/* Transform the velocity sample to displacement sample */
			integral_waveforms( &snl_infos[i], true );
		/* */
			proc_disp( &snl_infos[i], end_pos );
		}
		else {
		/* Transform the acceleration sample directly to displacement sample */
			integral_waveforms( &snl_infos[i], false );
			integral_waveforms( &snl_infos[i], true );
		/* */
			proc_disp( &snl_infos[i], end_pos );
		/* Then transform the displacement sample back to velocity sample */
			differential_waveform( &snl_infos[i] );
		/* */
			proc_vel( &snl_infos[i], end_pos );
		}
	/* Computation of Tau-c at 3 seconds */
		snl_infos[i].tc = calc_tau_c(
			&snl_infos[i].sum_dis[snl_infos[i].parrival_pos],
			&snl_infos[i].sum_vel[snl_infos[i].parrival_pos],
			end_pos, snl_infos[i].delta, 3
		);
	/* Finally, derive the lead time information */
		proc_leadtime( &snl_infos[i] );

	/* End of seismic data processing */
		fprintf(
			stderr, "Finished the processing data of %s.%s.%s (start at %lf, npts %d, delta %.2lf)!\n",
			snl_infos[i].sta, snl_infos[i].net, snl_infos[i].loc, snl_infos[i].starttime, snl_infos[i].npts, snl_infos[i].delta
		);

	/* After the processing, free the seismic data memory space */
		for ( int j = 0; j < NUM_CHANNEL_SNL; j++ ) {
			if ( snl_infos[i].seis[j] ) {
				free(snl_infos[i].seis[j]);
				snl_infos[i].seis[j] = NULL;
			}
		}
	/* */
		if ( snl_infos[i].sum_vel )
			free(snl_infos[i].sum_vel);
		if ( snl_infos[i].sum_dis )
			free(snl_infos[i].sum_dis);
	}

/* Output the result, first the header... */
	if ( HeaderSwitch ) {
		fprintf(stdout, OUTPUT_FILE_HEADER);
		if ( CoordinateSwitch )
			fprintf(stdout, OUTPUT_FILE_COOR_HEADER);
		fprintf(stdout, "\n");
	}
/* Then, all the stations' result */
	for ( int i = 0; i < totalsnl; i++ ) {
	/* */
		if ( (IgnStaWithoutData && snl_infos[i].npts < 0) || (IgnStaWithoutPick && !snl_infos[i].pick_flag) )
			continue;
	/* */
		fprintf(
			stdout, OUTPUT_DATA_FORMAT,
			snl_infos[i].sta, snl_infos[i].net, snl_infos[i].loc,
			snl_infos[i].pga, snl_infos[i].pgv, snl_infos[i].pgd,
			snl_infos[i].pa3, snl_infos[i].pv3, snl_infos[i].pd3, snl_infos[i].tc,
			snl_infos[i].pga_leadtime, snl_infos[i].pgv_leadtime, snl_infos[i].epic_dist, snl_infos[i].snr
		);
		if ( CoordinateSwitch )
			fprintf(
				stdout, OUTPUT_DATA_COOR_FORMAT,
				snl_infos[i].latitude, snl_infos[i].longitude, snl_infos[i].elevation
			);
		fprintf(stdout, "\n");
	}

/* */
	free(snl_infos);

	return 0;
}

/**
 * @brief
 *
 * @param argc
 * @param argv
 * @return int
 */
static int proc_argv( int argc, char *argv[] )
{
	char informat[16] = { 0 };

	for ( int i = 1; i < argc; i++ ) {
		if ( !strcmp(argv[i], "-v") ) {
			fprintf(stdout, "%s\n", PROG_NAME);
			fprintf(stdout, "Version: %s\n", VERSION);
			fprintf(stdout, "Author:  %s\n", AUTHOR);
			fprintf(stdout, "Compiled at %s %s\n", __DATE__, __TIME__);
			exit(0);
		}
		else if ( !strcmp(argv[i], "-h") ) {
			usage();
			exit(0);
		}
		else if ( !strcmp(argv[i], "-c") ) {
			CoordinateSwitch = true;
		}
		else if ( !strcmp(argv[i], "-n") ) {
			HeaderSwitch = false;
		}
		else if ( !strcmp(argv[i], "-t") ) {
			TwoStageIntegral = true;
		}
		else if ( !strcmp(argv[i], "-s") ) {
			VecSumSwitch = true;
		}
		else if ( !strcmp(argv[i], "-i") ) {
			IgnStaWithoutData = true;
		}
		else if ( !strcmp(argv[i], "-ip") ) {
			IgnStaWithoutPick = true;
		}
		else if ( !strcmp(argv[i], "-f") ) {
			strncpy(informat, argv[++i], sizeof(informat) - 1);
			for ( char *c = informat; *c; c++ )
				*c = toupper(*c);
		}
		else if ( i == argc - 1 ) {
			SeisDataFile = argv[i];
		}
		else if ( i == argc - 2 ) {
			StaListFile = argv[i];
		}
		else if ( i == argc - 3 ) {
			EqInfoFile = argv[i];
		}
		else {
			fprintf(stderr, "Unknown option: %s\n\n", argv[i]);
			return -1;
		}
	}
/* */
	if ( !SeisDataFile || !StaListFile || !EqInfoFile ) {
		fprintf(stderr, "No input file was specified; ");
		fprintf(stderr, "exiting with error!\n\n");
		return -1;
	}
/* */
	if ( !strlen(informat) ) {
		strcpy(informat, "SAC");
	}
/* */
	if ( !strcmp(informat, "SAC") ) {
		LoadSeisdataFunc = seisdata_load_sac;
	}
	else if ( !strcmp(informat, "MSEED") || !strcmp(informat, "MSEED3") ) {
		LoadSeisdataFunc = seisdata_load_ms;
	}
	else if ( !strcmp(informat, "TANK") ) {
		LoadSeisdataFunc = seisdata_load_tank;
	}
	else {
		fprintf(stderr, "Unknown format: %s\n", informat);
		return -1;
	}

	return 0;
}

/**
 * @brief
 *
 */
static void usage( void )
{
	fprintf(stdout, "\n%s\n", PROG_NAME);
	fprintf(stdout, "Version: %s\n", VERSION);
	fprintf(stdout, "Author:  %s\n", AUTHOR);
	fprintf(stdout, "Compiled at %s %s\n", __DATE__, __TIME__);
	fprintf(stdout, "***************************\n");
	fprintf(stdout, "Usage: %s [options] <input eq. info> <input station list> <input seismic data>\n", PROG_NAME);
	fprintf(stdout, "       or %s [options] <input eq. info> <input station list> <input seismic data> > <output path>\n\n", PROG_NAME);
	fprintf(stdout,
		"*** Options ***\n"
		" -v              Report program version\n"
		" -h              Show this usage message\n"
		" -c              Append the station coordinate in output, default is off\n"
		" -n              Turn off the output header, default is on\n"
		" -t              Turn on the two stage integral process, default is only one stage\n"
		" -s              Turn on the vector summation process, default is off\n"
		" -i              Ignore the station without input seismic data, default is on\n"
		" -ip             Ignore the station without valid picking, default is on\n"
		" -f format       Specify input format, there are SAC, MSEED|MSEED3 & TANK, default is SAC\n"
		//" -o output_file  Specify output file name, it will turn off the standard output & create a new output file\n"
		"\n"
		"This program will program to read SAC data files and compute\n"
		"the peak values include acceleration, velocity & displacement.\n"
		"It can also derive the warning lead time for each station.\n"
		"\n"
	);

	return;
}

/**
 * @brief
 *
 * @param snl_info
 */
static void init_snl_info_params( SNL_INFO *snl_info )
{
/* */
	snl_info->npts         = -1;
	snl_info->delta        = -1.0;
	snl_info->starttime    = -1.0;
	snl_info->pick_flag    = 0;
	snl_info->parrival_pos = -1;
	snl_info->sarrival_pos = -1;
	snl_info->snr          = 0.0;
/* Derived from waveform */
	snl_info->pga = 0.0;
	snl_info->pgv = 0.0;
	snl_info->pgd = 0.0;
	snl_info->pa3 = 0.0;
	snl_info->pv3 = 0.0;
	snl_info->pd3 = 0.0;
	snl_info->tc  = 0.0;
/* */
	snl_info->sum_vel = NULL;
	snl_info->sum_dis = NULL;
/* */
	snl_info->pga_pos   = -1;
	snl_info->pgv_pos   = -1;
	snl_info->pgd_pos   = -1;
	snl_info->pd35_pos  = -1;
	snl_info->pga80_pos = -1;
	snl_info->pga4_pos  = -1;
/* */
	snl_info->pga_leadtime = -1.0;
	snl_info->pgv_leadtime = -1.0;

	return;
}

/**
 * @brief
 *
 * @param snl_info
 * @param line
 * @return int
 */
static int parse_stalist_line( SNL_INFO *snl_info, const char *line )
{
/* */
	if ( strlen(line) ) {
		for ( int i = 0; i < MAX_STR_SIZE; i++ ) {
			if ( line[i] == '#' || line[i] == '\n' ) {
				break;
			}
			else if ( line[i] == '\t' || line[i] == ' ' ) {
				continue;
			}
			else if (
				sscanf(
					line, "%s %s %s %f %f %f %s %e %s %e %s %e\n",
					snl_info->sta, snl_info->net, snl_info->loc, &snl_info->latitude, &snl_info->longitude, &snl_info->elevation,
					snl_info->chan[0], &snl_info->gain[0], snl_info->chan[1], &snl_info->gain[1], snl_info->chan[2], &snl_info->gain[2]
				) == 12
			) {
					return 1;
			}
		/* */
			break;
		}
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
static int parse_stalist( SNL_INFO **snl_info, const char *path )
{
	FILE *fd = NULL;
	char  line[MAX_STR_SIZE] = { 0 };
	int   totalline = 0;

/* */
	if ( (fd = fopen(path, "r")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening station list %s\n", path);
		return -1;
	}

/* */
	if ( !(*snl_info = (SNL_INFO *)calloc(1, sizeof(SNL_INFO))) ) {
		fprintf(stderr, "Error allocating memory space for SNLs!\n");
		totalline = -1;
		goto close_list;
	}
/* */
	while ( fgets(line, sizeof(line) - 1, fd) != NULL )
		if ( parse_stalist_line( *snl_info, line ) )
			totalline++;
/* */
	rewind(fd);

/* */
	if ( totalline > 1 ) {
		free(*snl_info);
		if ( !(*snl_info = (SNL_INFO *)calloc(totalline + 1, sizeof(SNL_INFO))) ) {
			fprintf(stderr, "Error allocating memory space for SNLs!\n");
			totalline = -1;
			goto close_list;
		}
	}
/* */
	totalline = 0;
	while ( fgets(line, sizeof(line) - 1, fd) != NULL ) {
		memset(*snl_info + totalline, 0 , sizeof(SNL_INFO));
		init_snl_info_params( *snl_info + totalline );
		if ( parse_stalist_line( *snl_info + totalline, line ) )
			totalline++;
	}

close_list:
	fclose(fd);

	return totalline;
}

/*
 *
 */
static int parse_eqinfo_file( const char *path, float *epc_lat, float *epc_lon, float *dep, double *otime )
{
	FILE *fd;
	char  line[MAX_STR_SIZE] = { 0 };
/* */
	float year;
	float mon;
	float day;
	float hour;
	float min;
	float sec;
	struct tm _otime;

/* */
	if ( (fd = fopen(path, "r")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening Eq. information %s\n", path);
		return -2;
	}
/* */
	while ( fgets(line, sizeof(line) - 1, fd) != NULL ) {
	/* */
		if ( strlen(line) ) {
			for ( int i = 0; i < MAX_STR_SIZE; i++ ) {
				if ( line[i] == '#' || line[i] == '\n' ) {
					break;
				}
				else if ( line[i] == '\t' || line[i] == ' ' ) {
					continue;
				}
				else if (
					sscanf(
						line, "%f %f %f %f %f %f %f %f %f\n",
						&year, &mon, &day, &hour, &min, &sec, epc_lat, epc_lon, dep
					) == 9
				) {
				/* Keep the column of second to zero, then add it back after converting */
					_otime.tm_year  = (int)year - 1900;
					_otime.tm_mon   = (int)mon - 1;
					_otime.tm_mday  = (int)day;
					_otime.tm_hour  = (int)hour;
					_otime.tm_min   = (int)min;
					_otime.tm_sec   = 0;
					_otime.tm_isdst = 0;
				/* If we do so, we can have fraction of second */
					*dep    = -(*dep);
					*otime  = (double)timegm(&_otime);
					*otime += sec;

					fclose(fd);
					return 0;
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

/**
 * @brief
 *
 * @param snl_info
 * @param origin_time
 * @return int
 */
static int pick_pwave_arrival( SNL_INFO *snl_info, const double origin_time )
{
	int start_pick = (int)((origin_time - snl_info->starttime) / snl_info->delta);

/* */
	if ( start_pick < 0 )
		start_pick = 0;
/* */
	snl_info->snr = 0.0;
	if ( (snl_info->parrival_pos = pickwu_p_arrival_pick( snl_info->seis[0], snl_info->npts, snl_info->delta, 2, start_pick )) ) {
		if ( pickwu_p_trigger_check( snl_info->seis[0], snl_info->npts, snl_info->delta, snl_info->parrival_pos ) ) {
			if ( pickwu_p_arrival_quality_calc( snl_info->seis[0], snl_info->npts, snl_info->delta, snl_info->parrival_pos, &snl_info->snr ) < 4 ) {
				return 1;
			}
		}
	}
/* */
	snl_info->parrival_pos = start_pick;

	return 0;
}

/**
 * @brief
 *
 * @param snl_info
 * @param end_pos
 */
static void proc_acc( SNL_INFO *snl_info, const int end_pos )
{
	float vec_sum[snl_info->npts];

/* */
	snl_info->pga80_pos = snl_info->pga4_pos = end_pos;
/* */
	if ( VecSumSwitch ) {
		for ( int j = snl_info->parrival_pos; j < end_pos; j++ ) {
			vec_sum[j] =
				snl_info->seis[0][j] * snl_info->seis[0][j] +
				snl_info->seis[1][j] * snl_info->seis[1][j] +
				snl_info->seis[2][j] * snl_info->seis[2][j];
			vec_sum[j] = sqrtf(vec_sum[j]);
		/* */
			if ( vec_sum[j] > 4.0 ) {
			/* */
				if ( j < snl_info->pga4_pos )
					snl_info->pga4_pos = j;
			/* */
				if ( vec_sum[j] > 80.0 ) {
					if ( j < snl_info->pga80_pos )
						snl_info->pga80_pos = j;
				}
			}
		/* */
			if ( vec_sum[j] > snl_info->pga ) {
				snl_info->pga     = vec_sum[j];
				snl_info->pga_pos = j;
			}
		}
	/* */
		snl_info->pa3 = calc_peak_value( &vec_sum[snl_info->parrival_pos], end_pos, snl_info->delta, 3 );
	}
	else {
		for ( int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
			for ( int j = snl_info->parrival_pos; j < end_pos; j++ ) {
			/* */
				if ( fabs(snl_info->seis[i][j]) > 4.0 ) {
				/* */
					if ( j < snl_info->pga4_pos )
						snl_info->pga4_pos = j;
				/* */
					if ( fabs(snl_info->seis[i][j]) > 80.0 ) {
						if ( j < snl_info->pga80_pos )
							snl_info->pga80_pos = j;
					}
				}
			/* */
				if ( fabs(snl_info->seis[i][j]) > snl_info->pga ) {
					snl_info->pga     = fabs(snl_info->seis[i][j]);
					snl_info->pga_pos = j;
				}
			}
		}
	/* */
		snl_info->pa3 = calc_peak_value( snl_info->seis[0] + snl_info->parrival_pos, end_pos, snl_info->delta, 3 );
	}


	return;
}

/**
 * @brief
 *
 * @param snl_info
 * @param end_pos
 */
static void proc_vel( SNL_INFO *snl_info, const int end_pos )
{
	float vec_sum[snl_info->npts];

/* */
	if ( VecSumSwitch ) {
		for ( int j = snl_info->parrival_pos; j < end_pos; j++ ) {
			vec_sum[j] =
				snl_info->seis[0][j] * snl_info->seis[0][j] +
				snl_info->seis[1][j] * snl_info->seis[1][j] +
				snl_info->seis[2][j] * snl_info->seis[2][j];
		/* */
			snl_info->sum_vel[j] = vec_sum[j] + snl_info->sum_vel[j > 0 ? j - 1 : 0];
		/* */
			vec_sum[j] = sqrtf(vec_sum[j]);
		/* */
			if ( vec_sum[j] > snl_info->pgv ) {
				snl_info->pgv     = vec_sum[j];
				snl_info->pgv_pos = j;
			}
		}
	/* */
		snl_info->pv3 = calc_peak_value( &vec_sum[snl_info->parrival_pos], end_pos, snl_info->delta, 3 );
	}
	else {
		for ( int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
			for ( int j = snl_info->parrival_pos; j < end_pos; j++ ) {
			/* */
				if ( fabs(snl_info->seis[i][j]) > snl_info->pgv ) {
					snl_info->pgv     = fabs(snl_info->seis[i][j]);
					snl_info->pgv_pos = j;
				}
			}
		}
	/* */
		for ( int i = snl_info->parrival_pos; i < end_pos; i++ )
			snl_info->sum_vel[i] = snl_info->seis[0][i] * snl_info->seis[0][i] + snl_info->sum_vel[i > 0 ? i - 1 : 0];
	/* */
		snl_info->pv3 = calc_peak_value( snl_info->seis[0] + snl_info->parrival_pos, end_pos, snl_info->delta, 3 );
	}

	return;
}

/**
 * @brief
 *
 * @param snl_info
 * @param end_pos
 */
static void proc_disp( SNL_INFO *snl_info, const int end_pos )
{
	float vec_sum[snl_info->npts];

/* */
	snl_info->pd35_pos = end_pos;
/* */
	if ( VecSumSwitch ) {
		for ( int j = snl_info->parrival_pos; j < end_pos; j++ ) {
			vec_sum[j] =
				snl_info->seis[0][j] * snl_info->seis[0][j] +
				snl_info->seis[1][j] * snl_info->seis[1][j] +
				snl_info->seis[2][j] * snl_info->seis[2][j];
		/* */
			snl_info->sum_dis[j] = vec_sum[j] + snl_info->sum_dis[j > 0 ? j - 1 : 0];
		/* */
			vec_sum[j] = sqrtf(vec_sum[j]);
		/* */
			if ( vec_sum[j] > 0.35 ) {
				if ( j < snl_info->pd35_pos )
					snl_info->pd35_pos = j;
			}
		/* */
			if ( vec_sum[j] > snl_info->pgd ) {
				snl_info->pgd     = vec_sum[j];
				snl_info->pgd_pos = j;
			}
		}
	/* Computation of Pd at 3 seconds */
		snl_info->pd3 = calc_peak_value( &vec_sum[snl_info->parrival_pos], end_pos, snl_info->delta, 3 );
	}
	else {
	/* */
		for ( int i = 0; i < NUM_CHANNEL_SNL; i++ ) {
			for ( int j = snl_info->parrival_pos; j < end_pos; j++ ) {
			/* */
				if ( fabs(snl_info->seis[i][j]) > 0.35 ) {
					if ( j < snl_info->pd35_pos )
						snl_info->pd35_pos = j;
				}
			/* */
				if ( fabs(snl_info->seis[i][j]) > snl_info->pgd ) {
					snl_info->pgd     = fabs(snl_info->seis[i][j]);
					snl_info->pgd_pos = j;
				}
			}
		}
	/* */
		for ( int i = snl_info->parrival_pos; i < end_pos; i++ )
			snl_info->sum_dis[i] = snl_info->seis[0][i] * snl_info->seis[0][i] + snl_info->sum_dis[i > 0 ? i - 1 : 0];
	/* Computation of Pd at 3 seconds */
		snl_info->pd3 = calc_peak_value( snl_info->seis[0] + snl_info->parrival_pos, end_pos, snl_info->delta, 3 );
	}

	return;
}

/**
 * @brief
 *
 * @param snl_info
 */
static void proc_leadtime( SNL_INFO *snl_info )
{
/* */
	if ( !snl_info->pick_flag || (snl_info->pd35_pos <= 0 && snl_info->pga80_pos <= 0) ) {
		snl_info->pga_leadtime = snl_info->pgv_leadtime = -1.0;
	/* Reset the P-wave peak value 'cause there is not valid arrival time */
		if ( !snl_info->pick_flag )
			snl_info->pa3 = snl_info->pv3 = snl_info->pd3 = snl_info->tc = -1.0;
	}
	else {
	/* */
		if ( (snl_info->pga_pos - snl_info->pd35_pos) <= (snl_info->pga_pos - snl_info->pga80_pos) )
			snl_info->pga_leadtime = (snl_info->pga_pos - snl_info->pga80_pos) * snl_info->delta;
		else
			snl_info->pga_leadtime = (snl_info->pga_pos - snl_info->pd35_pos) * snl_info->delta;
	/* */
		if ( (snl_info->pgv_pos - snl_info->pd35_pos) <= (snl_info->pgv_pos - snl_info->pga80_pos) )
			snl_info->pgv_leadtime = (snl_info->pgv_pos - snl_info->pga80_pos) * snl_info->delta;
		else
			snl_info->pgv_leadtime = (snl_info->pgv_pos - snl_info->pd35_pos) * snl_info->delta;
	/* */
		if ( snl_info->pga_leadtime < 0.0 )
			snl_info->pga_leadtime = 0.0;
		if ( snl_info->pgv_leadtime < 0.0 )
			snl_info->pgv_leadtime = 0.0;
	}

	return;
}

/**
 * @brief
 *
 * @param snl_info
 * @param filter_sw
 */
static void integral_waveforms( SNL_INFO *snl_info, const _Bool filter_sw )
{
	for ( int i = 0; i < NUM_CHANNEL_SNL; i++ )
		_integral_waveform( snl_info->seis[i], snl_info->npts, snl_info->delta, filter_sw );

	return;
}

/**
 * @brief
 *
 * @param snl_info
 */
static void differential_waveform( SNL_INFO *snl_info )
{
	for ( int i = 0; i < NUM_CHANNEL_SNL; i++ )
		_differential_waveform( snl_info->seis[i], snl_info->npts, snl_info->delta );

	return;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param delta
 * @param filter_sw
 * @return float*
 */
static float *_integral_waveform( float *input, const int npts, const double delta, const _Bool filter_sw )
{
	const float half_delta = delta * 0.5;

	float last_seis  = 0.0;
	float this_pseis = 0.0;

/* */
	for ( int i = 0; i < npts; i++ ) {
		this_pseis = (input[i] + last_seis) * half_delta + this_pseis;
		last_seis  = input[i];
		input[i]   = this_pseis;
	}
/* */
	if ( filter_sw )
		highpass_filter( input, npts, delta, false );

	return input;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param delta
 * @return float*
 */
static float *_differential_waveform( float *input, const int npts, const double delta )
{
	float last_seis  = 0.0;
	float this_pseis = 0.0;

/* */
	for ( int i = 0; i < npts; i++ ) {
		this_pseis = (input[i] - last_seis) / delta;
		last_seis  = input[i];
		input[i]   = this_pseis;
	}

	return input;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param delta
 * @param zero_phase
 * @return float*
 */
static float *highpass_filter( float *input, const int npts, const double delta, const _Bool zero_phase )
{
	IIR_FILTER filter;
	IIR_STAGE *stage;

	filter = iirfilter_design( 2, IIR_HIGHPASS_FILTER, IIR_BUTTERWORTH, 0.075, 0.0, delta );
	stage  = (IIR_STAGE *)calloc(filter.nsects, sizeof(IIR_STAGE));

/* First time, forward filtering */
	memset(stage, 0, sizeof(IIR_STAGE) * filter.nsects);
	for ( int i = 0; i < npts; i++ )
		input[i] = iirfilter_apply( input[i], &filter, stage );
/* Second time, backward filtering */
	if ( zero_phase ) {
		memset(stage, 0, sizeof(IIR_STAGE) * filter.nsects);
		for ( int i = npts - 1; i >= 0; i-- )
			input[i] = iirfilter_apply( input[i], &filter, stage );
	}
/* */
	free(stage);

	return input;
}

/**
 * @brief
 *
 * @param sum_dis
 * @param sum_vel
 * @param npts
 * @param delta
 * @param sec
 * @return float
 */
static float calc_tau_c( const float *sum_dis, const float *sum_vel, const int npts, const float delta, const int sec )
{
	const int i_end = (int)(sec / delta);

	float result = 0.0;

/* */
	if ( i_end < npts ) {
		result = PI2 * sqrt(sum_dis[i_end] / sum_vel[i_end]);
	}

	return result > 10.0 ? 10.0 : result;
}

/**
 * @brief
 *
 * @param input
 * @param npts
 * @param delta
 * @param sec
 * @return float
 */
static float calc_peak_value( const float *input, const int npts, const float delta, const int sec )
{
	const int i_end = (int)(sec / delta);

	float result = 0.0;

/* */
	if ( i_end > npts )
		return 0.0;

	for ( int i = 1; i < i_end; i++ ) {
	/* */
		if ( fabs(input[i]) > result )
			result = fabs(input[i]);
	}

	return result;
}

/**
 * @brief  Transforms the coordinate(latitude & longitude) into distance(unit: km)
 *
 * @param elon
 * @param elat
 * @param slon
 * @param slat
 * @return double
 */
static double coor2distf( const double elon, const double elat, const double slon, const double slat )
{
	const double avlat = (elat + slat) * 0.5;

	double a = 1.840708 + avlat * (.0015269 + avlat * (-.00034 + avlat * (1.02337e-6)));
	double b = 1.843404 + avlat * (-6.93799e-5 + avlat * (8.79993e-6 + avlat * (-6.47527e-8)));

	a *= (slon - elon) * 60.0;
	b *= (slat - elat) * 60.0;

	return sqrt(a * a + b * b);
}
