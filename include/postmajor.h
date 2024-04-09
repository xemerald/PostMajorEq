/**
 * @file postmajor.h
 * @author Benjamin Yang @ National Taiwan University (b98204032@gmail.com)
 * @brief
 * @version 2.0.0
 * @date 2024-04-05
 *
 * @copyright Copyright (c) 2024
 *
 */

#pragma once

#include <sachead.h>
/* */
#define PROG_NAME       "postmajor"
#define VERSION         "2.0.0 - 2024-04-05"
#define AUTHOR          "Benjamin Ming Yang"
/* */
#define  OUTPUT_FILE_HEADER \
		"#SNL          PGA         PGV         PGD         PA3         PV3         PD3         TauC3       " \
		"PGA_LT      PGV_LT      E_Dist         SNR"
#define  OUTPUT_FILE_COOR_HEADER \
		"         LAT         LON          ELEV"
/* */
#define  OUTPUT_DATA_FORMAT \
		"%s.%s.%s %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %14.6lf"
#define  OUTPUT_DATA_COOR_FORMAT \
		" %11.6lf %11.6lf %8.2lf"
/* */
#define  NUM_CHANNEL_SNL  3
#define  EV_DURATION      180
#define  MAX_STR_SIZE     512
/* */
#define  PI  3.141592653589793238462643383279f
#define  PI2 6.283185307179586476925286766559f

/**
 * @brief
 *
 */
typedef struct {
/* */
	char sta[K_LEN + 1];
	char net[K_LEN + 1];
	char loc[K_LEN + 1];
	char chan[NUM_CHANNEL_SNL][K_LEN + 1];
/* */
	float gain[NUM_CHANNEL_SNL];
	float latitude;
	float longitude;
	float elevation;
/* */
	float *seis[3];       /* input trace buffer */
	int    npts;
	float  delta;
	double starttime;
/* Derived from waveform */
	int    parrival_pos;
	int    sarrival_pos;
	double snr;
/* */
	float pga;
	float pgv;
	float pgd;
	float pa3;
	float pv3;
	float pd3;
	float tc;
/* */
	int pga_pos;
	int pgv_pos;
	int pgd_pos;
	int pd35_pos;
	int pga80_pos;
	int pga4_pos;
/* */
	float pga_leadtime;
	float pgv_leadtime;
/* */
	float epic_dist;
} SNL_INFO;