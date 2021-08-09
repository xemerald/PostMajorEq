/*
 * Standalone program to read SAC data files and write
 * earthworm TRACE_BUF2 messages.
 * That file can then be made into a tankplayer file using remux_tbuf.
 *
 * Pete Lombard; May 2001
 */

#define VERSION_NUM  "0.0.2 2013-06-14"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <errno.h>
#include "trace_buf.h"
#include "sachead.h"
#include "swap.h"
#include "time_ew.h"

#define DEF_MAX_STA 700

/* Internal Structure Prototypes */
double delaz( double, double, double, double );
int locpt( float, float, float *, float *, int, int, int * );

typedef struct {
	/* Station profile */
	char station[8];
	float longitude;
	float latitude;
	/*double altitude;*/ /* Unused for now */
	float pga;
	float pgv;
} STAINFO;

int main( int argc, char **argv )
{
	int i, j, count;
	int zonecount;
	float x[1500], y[1500];
	FILE *fp;

	STAINFO stainfo[DEF_MAX_STA];
	STAINFO tmp;

	/* Check command line arguments
	******************************/
	if ( argc != 2 )
	{
		fprintf( stderr, "Usage: <Inputfile>\n" );
		exit( 0 );
	}

	/* Initialize the station info array
	************************************/
	memset( stainfo, 0, sizeof(stainfo) );
	memset( &tmp, 0, sizeof(tmp) );
	count = 0;

	if ( (fp = fopen("taiwan.txt","r")) == NULL )
	{
		printf("Error opening zone file!\n");
		exit(-1);
	}
	else
	{
		zonecount = 0;
		while( fscanf(fp,"%f %f", &x[zonecount], &y[zonecount]) == 2 ) zonecount++;
		fclose(fp);
		printf("Reading map zone file finish!\n");
	}

	/* Open the result file
	***********************/
	if ( (fp = fopen( argv[1], "rb" )) == (FILE *)NULL)
	{
		fprintf(stderr, "Error opening %s\n", argv[1]);
		exit( 1 );
	}

	while ( fread( &tmp, sizeof(STAINFO), 1, fp ) == 1 )
	{
		for ( i = 0; i < DEF_MAX_STA; i++ )
		{
			if ( strcmp( tmp.station, stainfo[i].station ) == 0 )
			{
				if ( tmp.pga > stainfo[i].pga ) stainfo[i].pga = tmp.pga;
				if ( tmp.pgv > stainfo[i].pgv ) stainfo[i].pgv = tmp.pgv;
				break;
			}
			else if ( stainfo[i].station[0] == 0 && stainfo[i].pga == 0.0 )
			{
				memcpy( &stainfo[i], &tmp, sizeof(STAINFO) );
				count++;
				break;
			}
			else continue;
		}
	}

	fclose(fp);

	FILE *Event_File;
	char outputname[128];
	//char *dotpos = strchr(argv[1], '.');

	int m, k, l, mm;
	int mcount, kcount;
	int count_25 = 0, count_80 = 0;
	int count_250 = 0, count_400 = 0;

	double lonmax = 0.0, lonmin = 180.0;
	double latmax = 0.0, latmin = 90.0;
	double tmplat = 0.0, tmplon = 0.0;

	double velsum = 0.0, voverd = 0.0;
	double galsum = 0.0, goverd = 0.0;
	double dists = 0.0, dissum = 0.0;
	double c_lat = 0.0, c_lon = 0.0, max_gal = 0.0;

	double m25 = 0.0, m80 = 0.0, m250 = 0.0, m400 = 0.0;

	for ( i = 0; i < count; i++ )
	{
		if ( stainfo[i].latitude <= latmin ) latmin = stainfo[i].latitude;
		if ( stainfo[i].latitude >= latmax ) latmax = stainfo[i].latitude;

		if ( stainfo[i].longitude <= lonmin ) lonmin = stainfo[i].longitude;
		if ( stainfo[i].longitude >= lonmax ) lonmax = stainfo[i].longitude;
	}

	latmin -= 0.05;
	lonmin -= 0.05;
	latmax += 0.05;
	lonmax += 0.05;

	kcount = (int)((lonmax - lonmin) / 0.01);
	mcount = (int)((latmax - latmin) / 0.01);

	//*dotpos = '\0';
	sprintf(outputname, "%s_Grid.rep", argv[1]);

	Event_File = fopen( outputname, "w" );
	fprintf( Event_File, "# lon lat acc vel\n" );

	for ( k = 0; k < kcount; k++ )   /*longitude*/
	{
		for ( m = 0; m < mcount; m++ )    /*latitude*/
		{
			tmplat = latmin + 0.01*m;
			tmplon = lonmin + 0.01*k;
			if ( locpt((float)tmplon, (float)tmplat, x, y, zonecount, l, &mm) > -1 )
			{
				dissum = 0.0;
				galsum = 0.0;
				goverd = 0.0;
				velsum = 0.0;
				voverd = 0.0;

				for ( i = 0; i < count; i++ )
				{
					dists = pow(delaz(tmplat, tmplon, stainfo[i].latitude, stainfo[i].longitude), 2);

					if ( dists <= 3600.0 )
					{
						if ( dists <= 1.0 )
						{
							galsum = stainfo[i].pga;
							velsum = stainfo[i].pgv;

							dissum = 1.0;
							break;
						}
						else
						{
							galsum += stainfo[i].pga / dists;
							velsum += stainfo[i].pgv / dists;

							dissum += 1.0 / dists;
						}
					}
				}

				goverd = galsum/dissum;
				voverd = velsum/dissum;


				if( goverd > max_gal )
				{
					max_gal = goverd;
					c_lat = tmplat;
					c_lon = tmplon;
				}

				if ( goverd >= 80.0 )
				{
					count_25++;
					if ( goverd >= 110.0 )
					{
						count_80++;
						if ( goverd >= 370.0 )
						{
							count_250++;
							if ( goverd >= 400.0 ) count_400++;
						}
					}
				}

				fprintf( Event_File, "%6.2f %6.2f %.4f %.4f\n", tmplon, tmplat, goverd, voverd );
			}
		}	/*for (m=0;m<mcount;m++)*/
	}	/*for (k=0;k<kcount;k++)*/

	fclose( Event_File );
/*
	sprintf(outputname, "%s_Mag.rep", argv[1]);
	Event_File = fopen( outputname, "w" );

	if ( count_25 > 0 )
	{
		m25 = (0.002248 * 80 + 0.279229) * log10(count_25*25) + 4.236343;

		if ( count_80 > 0 )
		{
			m80 = (0.002248 * 110 + 0.279229) * log10(count_80*25) + 4.236343;

			if ( count_250 > 0 )
			{
				m250 = (0.002248 * 370 + 0.279229) * log10(count_250*25) + 4.236343;

				if ( count_400 > 0 )
				{
					m400 = (0.002248 * 400 + 0.279229) * log10(count_400*25) + 4.236343;
				}
			}
		}
	}

	fprintf( Event_File, "%4.1f %4.1f %4.1f %4.1f %4d %4d %4d %4d %3d %6.2f %5.2f %6.2f\n",
			 m25, m80, m250, m400, count_25*25, count_80*25, count_250*25, count_400*25,
			 count, c_lon, c_lat, max_gal );

	fclose( Event_File );
*/
	sprintf(outputname, "%s_Sta.rep", argv[1]);
	Event_File = fopen( outputname, "w" );

	for ( i = 0; i < count; i++ )
	{
		fprintf( Event_File, "%8s %6.2f %6.2f %6.2f %6.2f\n",
				 stainfo[i].station, stainfo[i].latitude, stainfo[i].longitude,
				 stainfo[i].pga, stainfo[i].pgv );
	}

	fclose( Event_File );

	return 0;
}

/***********************************************************************************
 * delaz() Transforms the coordinate(latitude & lontitude) into distance(unit: km) *
 ***********************************************************************************/
double delaz(double elat,double elon,double slat, double slon)
{
	double delta;
	double avlat,a,b,dlat,dlon,dx,dy;

	avlat=0.5*(elat+slat);
	a=1.840708+avlat*(.0015269+avlat*(-.00034+avlat*(1.02337e-6)));
	b=1.843404+avlat*(-6.93799e-5+avlat*(8.79993e-6+avlat*(-6.47527e-8)));
	dlat=slat-elat;
	dlon=slon-elon;
	dx=a*dlon*60.0;
	dy=b*dlat*60.0;
	delta=sqrt(dx*dx+dy*dy);
	return delta;
}

/*****************************************************************************
 * locpt() Checks the position is whether inside the polygonal curve or not  *
 *****************************************************************************/
int locpt (float ex, float ey, float *x, float *y, int n, int l, int *m)
{
	int i;
	int n0 = n;
	double eps = 0.000001;
	double angle, pi, pi2, sum, theta[1500], tol, u, v;


	if ( *x == *(x+n-1) && *y == *(y+n-1) ) n0 = n - 1;

	pi = atan2(0.0, -1.0);
	pi2 = 2.0*pi;
	tol = 4.0*eps*pi;
	l = -1;
	*m = 0;

	u = *x - ex;
	v = *y - ey;

	if ( u == 0.0 && v == 0.0 ) 
	{
		l = 0;
		return l;
	} 
	
	if ( n0 < 2 ) return l;
	
	theta[1] = atan2(v, u);
	sum = 0.0;
	theta[0] = theta[1];

	for( i=1; i<n0; i++ )
	{
		u = *(x+i) - ex;
		v = *(y+i) - ey;

		if (u == 0.0 && v == 0.0) 
		{
			l = 0;
			return l;
		}

		theta[i] = atan2(v, u);
  
		angle = fabs(theta[i] - theta[0]);
		
		if (fabs(angle - pi) < tol)
		{
			l = 0;
			return l;
		}

		if (angle > pi) angle = angle - pi2;
		if (theta[0] > theta[i]) angle = -angle;
		
		sum = sum + angle;
		theta[0] = theta[i];
	}

	angle = fabs(theta[1] - theta[0]);
	
	if (fabs(angle - pi) < tol)
	{
		l = 0;
		return l;
	}
	
	if (angle > pi) angle = angle - pi2;
	if (theta[0] > theta[1]) angle = -angle;

	sum = sum + angle;
	*m = (int)(fabs(sum)/pi2 + 0.2);
	
	if (*m == 0) return l;
	l = 1;
	if (sum < 0.0) *m = -(*m);
	return l;
}
