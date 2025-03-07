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
#define VERSION         "2.1.0 - 2025-01-16"
#define AUTHOR          "Benjamin Ming Yang"
/* */
#define OUTPUT_FILE_HEADER \
		"#SNL          PGA         PGV         PGD         PA3         PV3         PD3         TauC3       " \
		"PGA_LT      PGV_LT      E_Dist         SNR"
#define OUTPUT_FILE_COOR_HEADER \
		"         LAT         LON          ELEV"
/* */
#define OUTPUT_DATA_FORMAT \
		"%s.%s.%s %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %14.6lf"
#define OUTPUT_DATA_COOR_FORMAT \
		" %11.6lf %11.6lf %8.2lf"
/* */
#define NUM_CHANNEL_SNL  3
#define EV_DURATION      180
#define MAX_STR_SIZE     512
/* */
#define PI  3.141592653589793238462643383279f
#define PI2 6.283185307179586476925286766559f
/* */
#define DEF_PD_WARN_THRESHOLD   0.35f
#define DEF_PGA_WARN_THRESHOLD  80.0f
#define DEF_PGA_WATCH_THRESHOLD 4.0f

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
	int    npts;
	int    gaps;
	float  delta;
	float *seis[3];       /* input trace buffer */
	double starttime;
/* Derived from waveform */
	int    pick_flag;
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
	float *sum_vel;
	float *sum_dis;
/* */
	int pga_pos;
	int pgv_pos;
	int pgd_pos;
	int pd_warn_pos;
	int pga_warn_pos;
	int pga_watch_pos;
/* */
	float pga_leadtime;
	float pgv_leadtime;
/* */
	float epic_dist;
} SNL_INFO;
