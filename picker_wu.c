#include <math.h>

/* */
#define PWAVE_TRIGGER 2.85
#define PWAVE_ARRIVE  1.25
#define PWAVE_STA     0.4
#define PWAVE_LTA     40.0
/* */
#define SWAVE_TRIGGER 3.0
#define SWAVE_ARRIVE  1.25
#define SWAVE_STA     0.5
#define SWAVE_LTA     3.0

/* */
static double characteristic_func_1( const double );
static double characteristic_func_2( const double, const double );

/*
 *
 */
int picker_wu(
	const double *input_z, const double *input_n, const double *input_e,
	const int np, const double delta, const int cf_flag, const int is,
	double *p_arrive, double *s_arrive, double *p_weight, double *s_weight
) {
	_Bool trigger = 0;
	int   i, j, tmp;
	int   ip_arr, ip_tri;
	int   is_arr, is_tri;
	int   ilta, ista;
	const int samprate = (int)(1.0 / delta);

	double *_input[3] = {
		input_z, input_n, input_e
	};
	double sum0, sum1, sum2;
	double ratio;
	double x_sta = 0.0;
	double x_lta = 0.0;

/* */
	tmp = np * 0.1;
	if ( tmp <= 0 )
		tmp = np;
/* */
	for ( i = 0; i < 3; i++ ) {
	/* */
		sum0 = 0.0;
		for ( j = 0; j < tmp; j++ )
			sum0 += _input[i][j];
	/* */
		sum0 /= (double)tmp;
		for ( j = 0; j < np; j++ )
			_input[i][j] -= sum0;
	}

/* Initial the arrivals & weightings */
	*p_arrive = 0.0;
	*s_arrive = 0.0;
	*p_weight = 4.0;
	*s_weight = 4.0;

/* Change P wave LTA & STA from seconds to sampling points */
	ilta = (int)(PWAVE_LTA / delta);
	ista = (int)(PWAVE_STA / delta);
/* Using ISTA+100 points to calculate initial STA */
	sum0 = 0.0;
	tmp = ista + 100;
	if ( tmp > np )
		tmp = np;
	for ( i = 0; i < tmp; i++ ) {
	/* Might need a switch statement here!!! */
		sum0 += characteristic_func_2( _input[0][i+1], _input[0][i] );
	}
	x_sta = sum0 / tmp;
	x_lta = x_sta * 1.25;

/* Start to detect P arrival, picking P arrival on V-component */
	j = tmp;
	for ( i = tmp; i < np; i++ ) {
	/* Might need a switch statement here!!! */
		sum0 = characteristic_func_2( _input[0][i], _input[0][i-1] );
	/* Update STA & LTA for each incoming data points */
		x_sta = (x_sta * (ista - 1) + sum0) / (double)ista;
		/* x_sta += (sum - x_sta) / ista */
		x_lta = (x_lta * (ilta - 1) + sum0) / (double)ilta;
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
			ip_arr = i;
	/*
	 * If STA/LTA ratio bigger than PWAVE_TRIGGER, declare
	 * P wave trigger, and exit this loop, and to go to
	 * next step to calculate P arrival's quality.
	 */
		//if ( i < is )
			//continue;
		if ( ratio > PWAVE_TRIGGER ) {
			trigger = 1;
			ip_tri  = i;
			break;
		}
	}

/*
 * If itrigger=0, it mean P wave did not trigger
 * no phase was picked on this station, return
 */
	if ( !trigger )
		return 0;
/* Change P arrival from data points to seconds */
	*p_arrive = ip_arr * delta;

/*
 * If (P arrival + 1 second) bigger than total data
 * length, it mean no enough data for calculating
 * picking quality, define P arrival's weighting = 3,
 * and return to main program
 */
	tmp = ip_arr + samprate;
	if ( tmp > np ) {
		*p_weight = 3.0;
		return 0;
	}

/*
 * Calculate P arrival's quality, P arrival's quality
 * is defined by the ratio of sum of 1 sec amplitude
 * square after P arrival over sum of 1 sec amplitude
 * aquare before P arrival.
 */

/* calculating sum of 1 sec amplitude square after P arrival */
	sum0 = 0.0;
	for ( i = ip_arr; i < tmp; i++ )
		sum0 += characteristic_func_1( _input[0][i] );
	sum0 /= (double)samprate;
/* calculating sum of 1 sec amplitude square before P arrival */
	sum1 = 0.0;
	tmp = ip_arr - samprate;
	if ( tmp < 0 )
		tmp = 0;
	for ( i = tmp; i < ip_arr; i++ )
		sum1 += characteristic_func_1( _input[0][i] );
	sum1 /= (double)(ip_arr - tmp);

/*
	Calculating the ratio
	If Ratio  >  30 P arrival's weighting define 0
	30 > R >  15 P arrival's weighting define 1
	15 > R >   3 P arrival's weighting define 2
	3 > R > 1.5 P arrival's weighting define 3
	1.5 > R       P arrival's weighting define 4
*/
	if ( sum1 > 0.0001 )
		ratio = sum0 / sum1;
	else
		ratio = 0.1;
/* */
	if ( ratio > 30.0 )
		*p_weight = 0.0;
	else if ( ratio > 15.0 )
		*p_weight = 1.0;
	else if ( ratio > 3.0 )
		*p_weight = 2.0;
	else if ( ratio > 1.5 )
		*p_weight = 3.0;
	else
		*p_weight = 4.0;

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

/*
 * If (P arrival + 3 sec) bigger than total data length,
 * it mean no enough data for checking these two types false trigger.
 */
	tmp = ip_arr + samprate * 3;
	if ( tmp > np )
		return 0;
/* calculating sum of 1 sec amplitude square after P arrival pass 2 sec */
	sum2 = 0.0;
	for ( i = ip_arr + samprate * 2; i < tmp; i++ )
		sum2 += characteristic_func_1( _input[0][i] );
	sum2 /= (double)samprate;
/*
 * Calculating the ratio of sum of 1 sec amplitude square after P arrival
 * pass 2 sec over the sum of 1 sec amplitude aquare before P arrival
 */
	if ( sum1 > 0.0001 )
		ratio = sum2 / sum1;
	else
		ratio = 0.1;
/*
 * If ratio less than 1.05, it may be spike type false trigger occurred,
 * remove this P arrival and return to main program.
 */
	if ( ratio < 1.05 ) {
		*p_arrive = 0.0;
		*p_weight = 4.0;
		*s_weight = 4.0;
		return 0;
	}

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
	for ( i = ip_arr; i < ip_arr + samprate; i++ )
		sum0 += _input[0][i];
	sum0 /= (double)samprate;
/* Calculating mean 1 sec before P arrival */
	sum1 = 0.0;
	tmp  = ip_arr - samprate * 2;
	if ( tmp < 0 )
		tmp = 0;
/* Notice for ip_arr!!!! */
	for ( i = tmp; i < ip_arr; i++ )
		sum1 += _input[0][i];
	sum1 /= (double)(ip_arr - tmp);
/*
 * If the difference bigger than 1.0, it may be the DC
 * drift false trigger occurred, remove this P arrival
 * and return to main program.
 */
	if ( fabs(sum0 - sum1) > 1.0 ) {
		*p_arrive = 0.0;
		*p_weight = 4.0;
		*s_weight = 4.0;
		return 0;
	}
/* End of detecting P arrival */

/*
 * Setup picking range, picking S arrival from 2 to 42
 * seconds after P trigger
 */
	tmp = ip_tri + samprate * 2;
	j   = tmp + samprate * 40;
/* Change S wave LTA & STA from seconds to sampling points */
	ilta = (int)(SWAVE_LTA / delta);
	ista = (int)(SWAVE_STA / delta);
/*
 * If Istart+2*ISTA bigger than total number, declare no
 * enough signal to pick S arrival,
 */
	if ( (tmp + ista * 2) > np )
		return 1;
/* Using ISTA points to calculate initial STA */
	sum0 = 0.0;
	for ( i = tmp; i < tmp + ista; i++ ) {
	/* */
		sum0 += characteristic_func_2( _input[1][i], _input[1][i-1] );
		sum0 += characteristic_func_2( _input[2][i], _input[2][i-1] );
	}
/* Initialize the LTA, is setted to equal X_STA */
	x_sta = sum0 / (double)ista;
	x_lta = x_sta;

/* Start to Pick S wave arrival */
	trigger = 0;
	is_arr  = 0;
	for ( i = tmp + ista; i < j; i++ ) {
	/* Might need a switch statement here!!! */
		sum0  = characteristic_func_2( _input[1][i], _input[1][i-1] );
		sum0 += characteristic_func_2( _input[1][i], _input[1][i-1] );
	/* Update STA & LTA for each incoming data points */
		x_sta = (x_sta * (ista - 1) + sum0) / (double)ista;
		/* x_sta += (sum - x_sta) / ista */
		x_lta = (x_lta * (ilta - 1) + sum0) / (double)ilta;
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
			is_arr = i;
	/*
	 * If STA/LTA ratio bigger than PWAVE_TRIGGER, declare
	 * P wave trigger, and exit to this loop, and to go to
	 * next step to calculate P arrival's quality.
	 */
		if ( ratio > SWAVE_TRIGGER ) {
			trigger = 1;
			is_tri  = i;
			break;
		}
	}
/*
 * If itrigger=0, it mean S wave did not trigger
 * no S phase picks on this station, return
 */
	if ( !trigger )
		return 1;
/*
 * if Isarr not settle when start picking
 * set Isarr at 0.1 seconds before S trigger
 */
	if ( is_arr == 0 )
		is_arr = is_tri - (int)(0.1 * samprate);
/* Transfer picking result from points to seconds */
	*s_arrive = is_arr * delta;
/*
 * If (Istri+ 1 sec) bigger than total number, declare no
 * enough signal to calculate S quality, set S_wei=3
 * and return to main program
 */
	if ( (is_tri + samprate) > np ) {
		*s_weight = 3.0;
		return 1;
	}

/*
 * Calculate S arrival's quality, S arrival's quality
 * is defined by the ratio of sum of 1 sec amplitude
 * square after S arrival over sum of 1 sec amplitude
 * square before S arrival, using 2-horizontial component
 */
/* calculating sum of 1 sec amplitude square after S arrival, using 2 horizontial component */
	sum0 = 0.0;
	for ( i = is_arr; i < is_arr + samprate; i++ )
		sum0 += characteristic_func_1( _input[1][i] ) + characteristic_func_1( _input[2][i] );
	sum0 /= (double)samprate;
/* calculating sum of 1 s amplitude square before S arrival, using 2 horizontial component */
	sum1 = 0.0;
	for ( i = is_arr - samprate; i < is_arr; i++ )
		sum1 += characteristic_func_1( _input[1][i] ) + characteristic_func_1( _input[2][i] );
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
		ratio = sum0 / sum1;
	else
		ratio = 0.1;
/* */
	if ( ratio > 30.0 )
		*s_weight = 0.0;
	else if ( ratio > 15.0 )
		*s_weight = 1.0;
	else if ( ratio > 5.0 )
		*s_weight = 2.0;
	else if ( ratio > 2.0 )
		*s_weight = 3.0;
	else
		*s_weight = 4.0;

	return 2;
}


/*
 *
 */
static double characteristic_func_1( const double sample )
{
	return sample * sample;
}

/*
 *
 */
static double characteristic_func_2( const double sample, const double sample_prev )
{
	double result;

/* */
	result  = sample - sample_prev;
	result *= result;
	result += sample * sample;

	return result;
}
