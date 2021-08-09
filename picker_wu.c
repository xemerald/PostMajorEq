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


double picker_wu(
	const double *input_z, const double *input_n, const double *input_e,
	const int np, const double delta,
	double *p_arrive, double *s_arrive, double *p_weigh, double *s_weigh
) {
	int    i, j, tmp;
	int    ilta, ista;
	double *_input[3] = {
		input_z, input_n, input_e
	};
	double sum;

/* */
	tmp = np * 0.1;
	if ( tmp <= 0 )
		tmp = np;
/* */
	for ( i = 0; i < 3; i++ ) {
	/* */
		sum = 0.0;
		for ( j = 0; j < tmp; j++ )
			sum += _input[i][j];
	/* */
		sum /= (double)tmp;
		for ( j = 0; j < np; j++ )
			_input[i][j] -= sum;
	}

/* Initial the arrivals & weightings */
	*p_arrive = 0.0;
	*s_arrive = 0.0;
	*p_weigh  = 4.0;
	*s_weigh  = 4.0;

/* Change P wave LTA & STA from seconds to sampling points */
	ilta = (int)(PWAVE_LTA / delta);
	ista = (int)(PWAVE_STA / delta);
/* Using ISTA+100 points to calculate initial STA */
	sum = 0.0;
	tmp = ista + 100;
	if ( tmp > np )
		tmp = np;
	for ( i = 0; i < tmp; i++ ) {
		sum += _input[0][i+1] * _input[0][i+1];
		sum += (_input[0][i+1] - _input[0][i]) * (_input[0][i+1] - _input[0][i]);
	}
	
}
