/*
 *
 */
#pragma once

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
int pickwu_p_arrival_pick( const float *, const int, const double, const int, const int );
int pickwu_s_arrival_pick( const float *, const float *, const int, const double, const int, const int );
int pickwu_p_arrival_quality_calc( const float *, const int, const double, const int, double * );
int pickwu_s_arrival_quality_calc( const float *, const float *, const int, const double, const int, double * );
int pickwu_p_trigger_check( const float *, const int, const double, const int );
