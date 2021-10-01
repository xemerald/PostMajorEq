/*
 *
 */

/* */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* */
#include <picker_wu.h>

/* Internal Function Prototypes */
static inline double characteristic_func_1( const double );
static inline double characteristic_func_2( const double, const double );

/*
 *
 */
int pickwu_p_arrival_pick(
	const float *input_z, const int np, const double delta, const int cf_flag, const int p_start
) {
	int   pos_trigger = 0;
	int   pos_tmp     = 0;
	int   result      = 0;
	int   i, ilta, ista;
	_Bool trigger = 0;

	double sum, ratio;
	double x_sta, x_lta;

/* Change P wave LTA & STA from seconds to sampling points */
	ilta = (int)(PWAVE_LTA / delta);
	ista = (int)(PWAVE_STA / delta);
/* Using ISTA+100 points to calculate initial STA */
	sum     = 0.0;
	pos_tmp = ista + 100;
	if ( pos_tmp > np )
		pos_tmp = np;
	for ( i = 0; i < pos_tmp; i++ ) {
	/* */
		switch ( cf_flag ) {
		case 2: default:
			sum += characteristic_func_2( input_z[i+1], input_z[i] );
			break;
		case 1:
			sum += characteristic_func_1( input_z[i] );
			break;
		}
	}
	x_sta = sum / pos_tmp;
	x_lta = x_sta * 1.25;

/* Start to detect P arrival, picking P arrival on V-component */
	for ( i = pos_tmp; i < np; i++ ) {
	/*  */
		switch ( cf_flag ) {
		case 2: default:
			sum = characteristic_func_2( input_z[i], input_z[i-1] );
			break;
		case 1:
			sum = characteristic_func_1( input_z[i] );
			break;
		}
	/* Update STA & LTA for each incoming data points */
		x_sta = (x_sta * (ista - 1) + sum) / (double)ista;
		/* x_sta += (sum - x_sta) / ista */
		x_lta = (x_lta * (ilta - 1) + sum) / (double)ilta;
	/* Set upper limit for LTA to avoid false trigger */
		if ( x_lta < 0.005 )
			x_lta = 0.005;
	/* Calculate STA/LTA ratio for check P arrival trigger */
		ratio = x_sta / x_lta;
	/*
	 * If STA/LTA ratio less than PWAVE_ARRIVE, keep this point,
	 * when trigger condition was met, this point define to P wave
	 * arrival.
	 */
		if ( ratio <= PWAVE_ARRIVE )
			result = i;
	/*
	 * If STA/LTA ratio bigger than PWAVE_TRIGGER, declare
	 * P wave trigger, and exit this loop, and to go to
	 * next step to calculate P arrival's quality.
	 */
		if ( i < p_start )
			continue;
		if ( ratio > PWAVE_TRIGGER ) {
			trigger     = 1;
			pos_trigger = i;
			break;
		}
	}

/*
 * If itrigger = 0, it mean P wave did not trigger
 * no phase was picked on this station, return
 */
	if ( !trigger )
		result = 0;

	return result;
}

/*
 *
 */
int pickwu_s_arrival_pick(
	const float *input_n, const float *input_e, const int np,
	const double delta, const int cf_flag, const int p_arrival
) {
	const int samprate = (int)(1.0 / delta);
/*
 * Setup picking range, picking S arrival from 2 to 42
 * seconds after P trigger
 */
	int   pos_2sec    = p_arrival + samprate * 2;
	int   pos_42sec   = pos_2sec + samprate * 40;
	int   pos_trigger = 0;
	int   result      = 0;
	int   i, ilta, ista;
	_Bool trigger = 0;

	double sum, ratio;
	double x_sta, x_lta;

/* */
	if ( input_n == NULL || input_e == NULL )
		return result;

/* Change S wave LTA & STA from seconds to sampling points */
	ilta = (int)(SWAVE_LTA / delta);
	ista = (int)(SWAVE_STA / delta);
/*
 * If Istart+2*ISTA bigger than total number, declare no
 * enough signal to pick S arrival,
 */
	if ( (pos_2sec + ista * 2) > np )
		return result;
/* Using ISTA points to calculate initial STA */
	sum = 0.0;
	for ( i = pos_2sec; i < pos_2sec + ista; i++ ) {
	/* */
		switch ( cf_flag ) {
		case 2: default:
			sum += characteristic_func_2( input_n[i], input_n[i-1] );
			sum += characteristic_func_2( input_e[i], input_e[i-1] );
			break;
		case 1:
			sum += characteristic_func_1( input_n[i] );
			sum += characteristic_func_1( input_e[i] );
			break;
		}
	}
/* Initialize the LTA, is setted to equal X_STA */
	x_sta = sum / (double)ista;
	x_lta = x_sta;

/* Start to Pick S wave arrival */
	for ( i = pos_2sec + ista; i < pos_42sec; i++ ) {
	/*  */
		switch ( cf_flag ) {
		case 2: default:
			sum += characteristic_func_2( input_n[i], input_n[i-1] );
			sum += characteristic_func_2( input_e[i], input_e[i-1] );
			break;
		case 1:
			sum += characteristic_func_1( input_n[i] );
			sum += characteristic_func_1( input_e[i] );
			break;
		}
	/* Update STA & LTA for each incoming data points */
		x_sta = (x_sta * (ista - 1) + sum) / (double)ista;
		/* x_sta += (sum - x_sta) / ista */
		x_lta = (x_lta * (ilta - 1) + sum) / (double)ilta;
	/* Set upper limit for LTA to avoid false trigger */
		if ( x_lta < 0.05 )
			x_lta = 0.05;
	/* Calculate STA/LTA ratio for check S arrival trigger */
		ratio = x_sta / x_lta;
	/*
	 * If STA/LTA ratio less than SWAVE_ARRIVE, keep this
	 * point, when trigger condition was met, this point
	 * define to S wave arrival.
	 */
		if ( ratio <= SWAVE_ARRIVE )
			result = i;
	/*
	 * If STA/LTA ratio bigger than PWAVE_TRIGGER, declare
	 * P wave trigger, and exit to this loop, and to go to
	 * next step to calculate P arrival's quality.
	 */
		if ( ratio > SWAVE_TRIGGER ) {
			trigger     = 1;
			pos_trigger = i;
			break;
		}
	}
/*
 * If itrigger=0, it mean S wave did not trigger
 * no S phase picks on this station, return
 */
	if ( !trigger )
		return result;
/*
 * if Isarr not settle when start picking
 * set Isarr at 0.1 seconds before S trigger
 */
	if ( !result )
		result = pos_trigger - (int)(0.1 * samprate);

	return result;
}

/*
 *
 */
int pickwu_p_arrival_quality_calc(
	const float *input_z, const int np, const double delta, const int p_arrival, double *snr
) {
	int   i, tmp;
	int   result;
	float sum0, sum1;

	const int samprate = (int)(1.0 / delta);

/* */
	*snr = 0.1;
/*
 * If (P arrival + 1 second) bigger than total data
 * length, it mean no enough data for calculating
 * picking quality, define P arrival's weighting = 3,
 * and return to main program
 */
	tmp = p_arrival + samprate;
	if ( tmp > np ) {
		result = 3;
		return result;
	}

/*
 * Calculate P arrival's quality, P arrival's quality
 * is defined by the ratio of sum of 1 sec amplitude
 * square after P arrival over sum of 1 sec amplitude
 * aquare before P arrival.
 */

/* calculating sum of 1 sec amplitude square after P arrival */
	sum0 = 0.0;
	for ( i = p_arrival; i < tmp; i++ )
		sum0 += characteristic_func_1( input_z[i] );
	sum0 /= (double)samprate;
/* calculating sum of 1 sec amplitude square before P arrival */
	sum1 = 0.0;
	tmp = p_arrival - samprate;
	if ( tmp < 0 )
		tmp = 0;
	for ( i = tmp; i < p_arrival; i++ )
		sum1 += characteristic_func_1( input_z[i] );
	sum1 /= (double)(p_arrival - tmp);

/*
	Calculating the ratio
	If Ratio  >  30 P arrival's weighting define 0
	30 > R >  15 P arrival's weighting define 1
	15 > R >   3 P arrival's weighting define 2
	3 > R > 1.5 P arrival's weighting define 3
	1.5 > R       P arrival's weighting define 4
*/
	if ( sum1 > 0.0001 )
		*snr = sum0 / sum1;

/* */
	if ( *snr > 30.0 )
		result = 0;
	else if ( *snr > 15.0 )
		result = 1;
	else if ( *snr > 3.0 )
		result = 2;
	else if ( *snr > 1.5 )
		result = 3;
	else
		result = 4;

	return result;
}

/*
 *
 */
int pickwu_s_arrival_quality_calc(
	const float *input_n, const float *input_e, const int np, const double delta, const int s_arrival, double *snr
) {
	int   i;
	int   pos_tmp;
	int   result;
	float sum0, sum1;

	const int samprate = (int)(1.0 / delta);

/* */
	*snr = 0.1;
/*
 * If (Istri+ 1 sec) bigger than total number, declare no
 * enough signal to calculate S quality, set S_wei=3
 * and return to main program
 */
	if ( (s_arrival + samprate) > np ) {
		result = 3;
		return result;
	}

/*
 * Calculate S arrival's quality, S arrival's quality
 * is defined by the ratio of sum of 1 sec amplitude
 * square after S arrival over sum of 1 sec amplitude
 * square before S arrival, using 2-horizontial component
 */
/* calculating sum of 1 sec amplitude square after S arrival, using 2 horizontial component */
	sum0 = 0.0;
	pos_tmp = s_arrival + samprate;
	for ( i = s_arrival; i < pos_tmp; i++ )
		sum0 += characteristic_func_1( input_n[i] ) + characteristic_func_1( input_e[i] );
	sum0 /= (double)samprate;
/* calculating sum of 1 s amplitude square before S arrival, using 2 horizontial component */
	sum1 = 0.0;
	for ( i = s_arrival - samprate; i < s_arrival; i++ )
		sum1 += characteristic_func_1( input_n[i] ) + characteristic_func_1( input_e[i] );
	sum1 /= (double)samprate;

/*
 * Calculate the ratio
 * If Ratio  >  30 S arrival's weighting define 0
 * 30 > R >  15 S arrival's weighting define 1
 * 15 > R >   5 S arrival's weighting define 2
 * 5 > R >   2 S arrival's weighting define 3
 * 2 > R       S arrival's weighting define 4
 */
	if ( sum1 > 0.001 )
		*snr = sum0 / sum1;

/* */
	if ( *snr > 30.0 )
		result = 0;
	else if ( *snr > 15.0 )
		result = 1;
	else if ( *snr > 5.0 )
		result = 2;
	else if ( *snr > 2.0 )
		result = 3;
	else
		result = 4;

	return result;
}

/*
 *
 */
int pickwu_p_trigger_check( const float *input_z, const int np, const double delta, const int p_arrival )
{
	int   i, tmp;
	float ratio;
	float sum0, sum1;

	const int samprate = (int)(1.0 / delta);
/*
 * Checking false trigger caused by transmission error!
 * Two type errors were occurred in our network,
 * spike    - Short time transmission error
 * DC drift - Long time transmission error
 *
 * In this program we uses the ratio of sum of 1 sec
 * amplitude square after P arrival pass 2 sec over the
 * sum of 1 sec amplitude aquare before P arrival, to
 * dectect SPIKE false trigger, if the ration less than 1.05
 * it may be SPIKE false trigger occurred, remove this P
 * arrival and return to main program.
 */


/* calculating sum of 1 sec amplitude square before P arrival */
	sum0 = 0.0;
	tmp = p_arrival - samprate;
	if ( tmp < 0 )
		tmp = 0;
	for ( i = tmp; i < p_arrival; i++ )
		sum0 += characteristic_func_1( input_z[i] );
	sum0 /= (double)(p_arrival - tmp);
/*
 * If (P arrival + 3 sec) bigger than total data length,
 * it mean no enough data for checking these two types false trigger.
 */
	tmp = p_arrival + samprate * 3;
	if ( tmp > np )
		return 0;
/* calculating sum of 1 sec amplitude square after P arrival pass 2 sec */
	sum1 = 0.0;
	for ( i = p_arrival + samprate * 2; i < tmp; i++ )
		sum1 += characteristic_func_1( input_z[i] );
	sum1 /= (double)samprate;
/*
 * Calculating the ratio of sum of 1 sec amplitude square after P arrival
 * pass 2 sec over the sum of 1 sec amplitude aquare before P arrival
 */
	if ( sum1 > 0.0001 )
		ratio = sum1 / sum0;
	else
		ratio = 0.1;
/*
 * If ratio less than 1.05, it may be spike type false trigger occurred,
 * remove this P arrival and return to main program.
 */
	if ( ratio < 1.05 )
		return 0;

/*
 * Detect DC drift type false trigger
 * in this program we check the 1 sec average after
 * P arrival minus 1 sec average before 1 sec of P
 * arrival If the difference bigger than 1.0 gal, it may be
 * DC drift false trigger, remove this pick and
 * return to main program
 */

/* Calculating mean after P arrival */
	sum0 = 0.0;
	for ( i = p_arrival; i < p_arrival + samprate; i++ )
		sum0 += input_z[i];
	sum0 /= (double)samprate;
/* Calculating mean 1 sec before P arrival */
	sum1 = 0.0;
	tmp  = p_arrival - samprate * 2;
	if ( tmp < 0 )
		tmp = 0;
/* Notice for ip_arr!!!! */
	for ( i = tmp; i < p_arrival; i++ )
		sum1 += input_z[i];
	sum1 /= (double)(p_arrival - tmp);
/*
 * If the difference bigger than 1.0, it may be the DC
 * drift false trigger occurred, remove this P arrival
 * and return to main program.
 */
	if ( fabs(sum0 - sum1) > 1.0 )
		return 0;
/* End of checking P arrival */

	return 1;
}

/*
 *
 */
static inline double characteristic_func_1( const double sample )
{
	return sample * sample;
}

/*
 *
 */
static inline double characteristic_func_2( const double sample, const double sample_prev )
{
	double result;

/* */
	result  = sample - sample_prev;
	result *= result;
	result += sample * sample;

	return result;
}
