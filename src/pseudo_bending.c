/*
	3D velocity earthquake location program
	Created 2009/07/26 by Yih-Min Wu
	Modifed 2017/08/11 by Benjamin Ming Yang
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

/* Used constant */
#define MSG_NUMBER   16384
#define PI           3.14159265358979323846
#define RAD_PER_DEG  0.01745329251994329577
#define DEG_PER_RAD  57.2957795130823208768
#define B2A_SQ       0.993305521
#define EPS          1.0e-8  /* Maybe 1.0e-10 */
#define RNULL        0.0e+10

#define HALFCYCLE    PI
#define FULLCYCLE    PI*2.0

#define GRID_NUMBER  6880000

/*
	parameters for calculation
		XFAC    = enhancement factor (see um and thurber, 1987)
		NLOOP   = number of bending iterations
		N1, N2  = min & max of ray segments
		MINS    = min. length of segment (km)
*/
#define XFAC         1.9
#define N1           2
#define N2           MSG_NUMBER
#define NLOOP        12800
#define FLIMIT       1.e-6
#define MINS         2.0

typedef struct {
	double rsq, sina, cosa;
	double x, y, z;
} NODE_INFO;

typedef struct {
	double r, a, b, v;
} RAY_INFO;

typedef struct {
	double vel[8];
} VEL_GRID;

/* For 3D-velocity model maping */
/*
static double Lat_c[100], Lon_c[100], Dep_c[100];
static double Vel_p[100][100][100];
static double Vel_s[100][100][100];
*/
static double *Lat_c, *Lon_c, *Dep_c;
static double *Vel_p, *Vel_s;
static VEL_GRID *Vel_grid;

static int *Ilon_c, *Ilat_c, *Idep_c;

/* Constant */
static double Ro, Rs;
static int Ilat1_c, Ilon1_c, Idep1_c;
static int Ilonmax, Ilatmax, Idepmax;
static int Ibld3, Ibld4;
static int Nlat_c, Nlon_c, Ndep_c;
static int Nxyz_c, Nxy_c, Nx_c;

/*
 *
 */
void pbr(double evla, double evlo, double evdp, double stla, double stlo, double stel, RAY_INFO *record, int *np, double *tk) {
	int ni, i, j, k, l;

	RAY_INFO *ray = NULL;

	RAY_INFO *ray_now = NULL;
	RAY_INFO *ray_prev = NULL;
	RAY_INFO *ray_next = NULL;
	RAY_INFO *ray_end = NULL;
	RAY_INFO ray_mid, calc_tmp;

	NODE_INFO node[2];

	double xfac;
	double shiftlo;
	double to, tn;

	double acosa, sina, rsina;
	//double cosa;

	double x1, y1, z1;
	double x2, y2, z2;
	double x3, y3, z3;

	double dn, dr, da, db;
	double dseg, ddseg;
	double upz, dwz;
	double vr, vb, va;
	double rvs, cc, rcur;

/*
	parameters for calculation
	initialization
*/
	xfac = XFAC;
/*
	ni : number of ray segments
*/
	ni = N1;

	ray = (RAY_INFO *)calloc(MSG_NUMBER + 1, sizeof(RAY_INFO)); /* Need to be optimize */

/* Random algorism used */
	//srand(time(NULL));

/* Check coordinates */
	if ( evla < -90.0 || evla > 90.0 ) {
		printf("Latitude of source is out of range!\n");
		exit(-1);
	}
	if ( stla < -90.0 || stla > 90.0 ) {
		printf("Latitude of station is out of range!\n");
		exit(-1);
	}

	if ( evlo < -180.0 || evlo > 180.0 ) {
		printf("Longitude of source is out of range!\n");
		exit(-1);
	}
	if ( stlo < -180.0 || stlo > 180.0 ) {
		printf("Longitude of station is out of range!\n");
		exit(-1);
	}

/*
	Longitude and latitude range from 0 to 180.
	This program does not work with angles
	greater than 180.

	Pass from latitude to colatitude
*/
	evla = (90.0 - evla) * RAD_PER_DEG;
	stla = (90.0 - stla) * RAD_PER_DEG;

/* da = bre, db = bso, dr = dlo */
	if ( stlo < 0.0 ) da = 360.0 + stlo;
	else da = stlo;

	if ( evlo < 0.0 ) db = 360.0 + evlo;
	else db = evlo;

	dr = fabs(db - da);
	shiftlo = 0.0;

	if ( dr < 180.0 ) {
		if ( db < da ) {
			evlo = (180.0 - dr) * 0.5;
			stlo = evlo + dr;
			shiftlo = db - evlo;
		}
		else {
			stlo = (180.0 - dr) * 0.5;
			evlo = stlo + dr;
			shiftlo = da - stlo;
		}
	}
	else {
		dr = 360.0 - dr;
		if ( db < da ) {
			evlo = (180.0 - dr) * 0.5 + dr;
			stlo = evlo - dr;
			shiftlo = db - evlo;
		}
		else {
			stlo = (180.0 - dr) * 0.5 + dr;
			evlo = stlo - dr;
			shiftlo = da - stlo;
		}
	}

	evlo *= RAD_PER_DEG;
	stlo *= RAD_PER_DEG;

	evdp = Ro - evdp;
	stel = Ro - stel;

/*
	Initial straight ray
*/
/* Epc. coordinates */
	node[0].x = evdp*sin(evla);
	node[0].y = node[0].x*sin(evlo);
	node[0].x *= cos(evlo);
	node[0].z = evdp*cos(evla);

/* Rec. coordinates */
	node[1].x = stel*sin(stla);
	node[1].y = node[1].x*sin(stlo);
	node[1].x *= cos(stlo);
	node[1].z = stel*cos(stla);

	dr = (node[1].x - node[0].x)/ni;
	da = (node[1].y - node[0].y)/ni;
	db = (node[1].z - node[0].z)/ni;

	ray_now = ray;
	ray_end = ray + ni;

	do {
		i = ray_now - ray;

		x1 = node[0].x + dr*(double)i;
		y1 = node[0].y + da*(double)i;
		z1 = node[0].z + db*(double)i;

		dn = x1*x1 + y1*y1 + z1*z1;
		ray_now->r = sqrt(dn + EPS);

		acosa = z1/ray_now->r;
		if ( acosa <= -1.0 )
			ray_now->a = PI;
		else if ( acosa >= 1.0 )
			ray_now->a = 0.0;
		else
			ray_now->a = acos(acosa);

		sina = sin(ray_now->a);

		acosa = x1/(sina*ray_now->r);
		if ( acosa <= -1.0 )
			ray_now->b = PI;
		else if ( acosa >= 1.0 )
			ray_now->b = 0.0;
		else
			ray_now->b = acos(acosa);

		if ( y1 < 0.0 ) ray_now->b = FULLCYCLE - ray_now->b;

		ray_now->v = velocity(*ray_now, shiftlo);

		if ( ray_now == ray || ray_now == ray_end ) {
			l = ray_now > ray ? 1 : 0;

			node[l].rsq = dn;
			node[l].sina = sina;
			node[l].cosa = cos(ray_now->a);
		}
	} while ( ++ray_now <= ray_end );

	tn = rtim(ni, ray, node);

/* iteration loop */
	while ( ni <= N2 ) {
		xfac = XFAC;
		l = ni - 1;
		for ( k = 0; k < NLOOP; k++ ) {
			if ( ni > 2 || k == 0 ) {
				for( j = 0; j < l; j++ ) {
				/* See Um & Thurber (1987) p.974. */
					if ( !(j & 0x01) )
						ray_now = ray + (j >> 1) + 1;
					else
						ray_now = ray_end - ((j + 1) >> 1);

					ray_prev = ray_now - 1;
					ray_next = ray_now + 1;

				/* Check if the first point of ray */
					if ( ray_prev > ray ) {
						x1 = ray_prev->r*sin(ray_prev->a);
						y1 = x1*sin(ray_prev->b);
						x1 *= cos(ray_prev->b);
						z1 = ray_prev->r*cos(ray_prev->a);
					}
					else {
						x1 = node[0].x;
						y1 = node[0].y;
						z1 = node[0].z;
					}

				/* Check if the last point of ray */
					if ( ray_next < ray_end ) {
						x3 = ray_next->r*sin(ray_next->a);
						y3 = x3*sin(ray_next->b);
						x3 *= cos(ray_next->b);
						z3 = ray_next->r*cos(ray_next->a);
					}
					else {
						x3 = node[1].x;
						y3 = node[1].y;
						z3 = node[1].z;
					}

				/* Check if the first time of loop & the first two perturbation,
				   if ture, we can just skip the calculation of cartesian coordinates. */
					if ( k > 0 || j >= 2 ) {
						x2 = x1 + x3;
						y2 = y1 + y3;
						z2 = z1 + z3;

						ray_mid.r = sqrt(x2*x2 + y2*y2 + z2*z2 + EPS);

						acosa = z2/ray_mid.r;
						if ( acosa <= -1.0 )
							ray_mid.a = PI;
						else if ( acosa >= 1.0 )
							ray_mid.a = 0.0;
						else
							ray_mid.a = acos(acosa);

						sina = sin(ray_mid.a);
						rsina = sina * ray_mid.r;
				        /* cosa = cos(ray_mid.a); */

						acosa = x2/rsina;
						if ( acosa <= -1.0 )
							ray_mid.b = PI;
						else if ( acosa >= 1.0 )
							ray_mid.b = 0.0;
						else
							ray_mid.b = acos(acosa);

						if ( y2 < 0.0 ) ray_mid.b = FULLCYCLE - ray_mid.b;
						ray_mid.r *= 0.5;
						rsina *= 0.5;

					/*
						Determine velocity at 3 points
					*/
						/* v1 = ray_prev->v; */
						ray_mid.v = velocity(ray_mid, shiftlo);
						/* v3 = ray_next->v; */
					}
					else {
						memcpy(&ray_mid, ray_now, sizeof(RAY_INFO));
						sina = sin(ray_mid.a);
						rsina = sina * ray_mid.r;
					}

					dr = x3 - x1;
					da = y3 - y1;
					db = z3 - z1;

					dn = dr * dr + da * da + db * db;
					dseg = sqrt(dn + EPS);

					dr = (ray_next->r - ray_prev->r); /* Orig. dr = (ray_next->r - ray_prev->r)/dseg */
					da = (ray_next->a - ray_prev->a); /* Orig. da = (ray_next->a - ray_prev->a)/dseg */
					db = (ray_next->b - ray_prev->b); /* Orig. db = (ray_next->b - ray_prev->b)/dseg */

				/*
					Now ddseg will be a distance to find dV
					along the coordinates
					Begin find the gradients and velocities
					first find the length of segment
				*/
					dseg *= 0.5;
					ddseg = dseg * 0.5;

				/*
					Begin to determine coordinates
					of pints surroundibg point a2, b2, r2
					at the distance ddseg
				*/
					upz = ray_mid.r + ddseg;
					dwz = ray_mid.r - ddseg;

					if ( upz > Rs ) {
					/* I guess it should be Ro + 10.0(Rs) */
						upz = Rs;
						dwz = upz - dseg;
					}

					if ( dwz <= EPS ) {
						dwz = EPS;
						upz = dwz + dseg;
					/*
						Set to Ro, the old mistake?
						upz = Ro;
					*/
					}

				/*
					The following if-endif is just for P & S, thus comment out for SKS & PKP!
					This gives the lowermost mantle Vp in the outer core
				*/
					memset(&calc_tmp, 0, sizeof(RAY_INFO));

					calc_tmp.a = ray_mid.a;
					calc_tmp.b = ray_mid.b;

					calc_tmp.r = upz;
					vr = velocity(calc_tmp, shiftlo);
					calc_tmp.r = dwz;
					vr -= velocity(calc_tmp, shiftlo);
					//vr = calc_tmp.v; /* Orig. vr = calc_tmp.v/dseg */

					km2deg(ray_mid, ddseg, RNULL, &calc_tmp);
					vb = velocity(calc_tmp, shiftlo);
					km2deg(ray_mid, -ddseg, RNULL, &calc_tmp);
					vb -= velocity(calc_tmp, shiftlo);
					//vb = calc_tmp.v; /* Orig. vb = calc_tmp.v/dseg */

					km2deg(ray_mid, RNULL, ddseg, &calc_tmp);
					va = velocity(calc_tmp, shiftlo);
					km2deg(ray_mid, RNULL, -ddseg, &calc_tmp);
					va -= velocity(calc_tmp, shiftlo);
					//va = calc_tmp.v; /* Orig. va = calc_tmp.v/dseg */

				/*
					spherical
					velocity gradient
					va = va / r2
					vb = vb / r2 / sina
					(tangential vector) = (slowness vector) / s
				*/
					/* dr *= 1 */
					da *= ray_mid.r;
					db *= rsina;
					rcur = dr*vr + da*va + db*vb;
					rcur /= dn;

					dr = vr - rcur*dr;
					da = va - rcur*da;
					db = vb - rcur*db;

					rvs = dr*dr + da*da + db*db;

					if ( rvs <= EPS ) {
					/* Special condition: velocity gradient equal to zero */
						memcpy(ray_now, &ray_mid, sizeof(RAY_INFO));
					} else {
						rvs = 1.0 / sqrt(rvs); /* sqrt(rvs + EPS) */

						cc = 0.5 / ray_prev->v + 0.5 / ray_next->v; /* cc = (1.0/v1 + 1.0/v3)/2.0 */
						rcur = vr * dr + va * da + vb * db;
						rcur /= dseg;
						rcur *= rvs;

					/*
						Tut esli rcur < 0.0 proishodit hernia
						poetomu postavlen abs. Ne yasno mozhno li eto delat
						ili net no rabotaet. Obichno oshibka poyavliaetsia
						ochen redko v nekotorih tochkah
						v etom sluchae abs prosto ne daet oshibki y posledniaya iteraciya
						uzhe ne imeet rcur negativnim y podgoniaet normalno reshenie
						(mozhet bit)
					*/
						if ( rcur < 0.0 ) {
							printf("Got negative Rc!\n");
							rcur = fabs(rcur);
						}

						ray_mid.v *= cc;
						rcur = -(ray_mid.v + 1.0)/(4.0*cc*rcur);
						rcur += sqrt(rcur*rcur + dn/(8.0*ray_mid.v) + EPS);

						rcur *= rvs;
						dr *= rcur;
						rcur /= rsina;
						da *= rcur * sina;
						db *= rcur;

						dr += ray_mid.r;
						da += ray_mid.a;
						db += ray_mid.b;

						ray_now->r = (dr - ray_now->r)*xfac + ray_now->r;
					/* if ray_now->r > 6371 then force it to the surface. */
						if ( ray_now->r > Rs )
							ray_now->r = Rs;
						ray_now->a = (da - ray_now->a)*xfac + ray_now->a;
						ray_now->b = (db - ray_now->b)*xfac + ray_now->b;
						ray_now->v = velocity(*ray_now, shiftlo);
					}
				}
			} else {
			/* ray_now = ray + (ni >> 1); */
				ray_now->r = (dr - ray_now->r)*xfac + ray_now->r;
			/* if ray_now->r > 6371 then force it to the surface. */
				if ( ray_now->r > Rs )
					ray_now->r = Rs;
				ray_now->a = (da - ray_now->a)*xfac + ray_now->a;
				ray_now->b = (db - ray_now->b)*xfac + ray_now->b;
				ray_now->v = velocity(*ray_now, shiftlo);
			}

			to = tn;
			tn = rtim(ni, ray, node);

			if ( fabs(to - tn) <= to*FLIMIT ) break;

		/* Random algorism, it will be faster but the result is not stable! */
			//xfac = (double)(rand()%10 + 10)/10.0;

			xfac *= 0.98;
			if ( xfac < 1.0 )
				xfac = 1.0;
			/* printf("%d %lf\n", k, xfac); */
		}

	/*
		Skip increasing of segment number if minimum length
		of segment is exceed or maximum number of segments
		was reached
	*/
		if ( dseg < MINS || ni >= N2 ) {
			/* igood = 1; */
			break;
		}

	/* Double the number of points. */
		ray_next = ray_end + 1;
		memcpy(ray_next, ray + 1, sizeof(RAY_INFO) * ni);
		for ( i=0; i<ni; i++ ) ray[(i+1)<<1] = ray_next[i];

		ni <<= 1;

		ray_end = ray + ni;

		x1 = node[0].x;
		y1 = node[0].y;
		z1 = node[0].z;

		ray_now = ray + 1;
		ray_next = ray_now + 1;

		do {
			if ( ray_next < ray_end ) {
				x3 = ray_next->r*sin(ray_next->a);
				y3 = x3*sin(ray_next->b);
				x3 *= cos(ray_next->b);
				z3 = ray_next->r*cos(ray_next->a);
			}
			else {
				x3 = node[1].x;
				y3 = node[1].y;
				z3 = node[1].z;
			}

			x2 = x3 + x1;
			y2 = y3 + y1;
			z2 = z3 + z1;

			ray_now->r = sqrt(x2*x2 + y2*y2 + z2*z2 + EPS);

			acosa = z2/ray_now->r;
			if ( acosa <= -1.0 ) ray_now->a = PI;
			else if ( acosa >= 1.0 ) ray_now->a = 0.0;
			else ray_now->a = acos(acosa);

			sina = sin(ray_now->a);

			acosa = x2/(sina*ray_now->r);
			if ( acosa <= -1.0 ) ray_now->b = PI;
			else if ( acosa >= 1.0 ) ray_now->b = 0.0;
			else ray_now->b = acos(acosa);

			if ( y2 < 0.0 ) ray_now->b = FULLCYCLE - ray_now->b;

			ray_now->r *= 0.5;

			ray_now->v = velocity(*ray_now, shiftlo);

			x1 = x3;
			y1 = y3;
			z1 = z3;

			ray_now += 2;
			ray_next += 2;
		} while ( ray_next <= ray_end );

		to = tn;
		tn = rtim(ni, ray, node);

		if ( fabs(to - tn) <= to*FLIMIT ) {
			/* igood = 1; */
			break;
		}
	}

/* Return coordinates to the origin */
/*
	ray_now = ray;
	for(i=0; i<=ni; i++) {
		ray_now.r = Ro - ray_now.r;
		ray_now.a = ray_now.a * DEG_PER_RAD;
		ray_now.a = geoc_to_geog(90.0 - ray_now.a);
		ray_now.b = ray_now.b * DEG_PER_RAD + shiftlo;
		if ( ray_now.b < 0.0 ) ray_now.b = 360.0 + ray_now.b;
		ray_now++;
	}
*/

/* Output the result */
	*np = ni + 1;
	*tk = tn;

/* Free all the used memory */
	ray_now = NULL;
	ray_prev = NULL;
	ray_next = NULL;
	ray_end = NULL;

	free(ray);

	return;
}

/*
 *
 */
static int input_vel( const char *modelfile )
{
	int     i, j, k;
	double  bld3, bld4;
	double  average;
	double *ptrtmp = NULL;
	char    fracline[1024] = { 0 };

	FILE   *fp = NULL;

/* */
	if ( (fp = fopen(modelfile, "r")) == NULL ) {
		printf("Reading VPVSMOD error; exiting!\n" );
		exit(-1);
	}

	if ( fscanf(fp, "%lf %lf %d %d %d\n", &bld3, &bld4, &Nlon_c, &Nlat_c, &Ndep_c) != 5 ) {
		printf("Reading VpVs Model header error; exiting!\n");
		exit(-1);
	}

	Lon_c = calloc(Nlon_c, sizeof(double));
	Lat_c = calloc(Nlat_c, sizeof(double));
	Dep_c = calloc(Ndep_c, sizeof(double));

	Nx_c = Nlon_c;
	Nxy_c = Nx_c * Nlat_c;
	Nxyz_c = Nxy_c * Ndep_c;

	Vel_p = calloc(Nxyz_c, sizeof(double));
	Vel_s = calloc(Nxyz_c, sizeof(double));
	Vel_grid = calloc(Nxyz_c, sizeof(VEL_GRID));

	if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
		ptrtmp = Lon_c;
		for ( i = 0; i < Nlon_c; i++ ) {
			if ( i < Nlon_c - 1 ) {
				if ( sscanf( fracline, " %lf %[^\n]", ptrtmp, fracline) != 2 ) {
					printf("Reading VpVs Model lon_c error; exiting!\n");
					exit(-1);
				}
			}
			else {
				if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
					printf("Reading VpVs Model lon_c error; exiting!\n");
					exit(-1);
				}
			}
			ptrtmp++;
		}
	}

	if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
		ptrtmp = Lat_c;
		for ( i = 0; i < Nlat_c; i++ ) {
			if ( i < Nlat_c - 1 ) {
				if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
					printf("Reading VpVs Model lat_c error; exiting!\n");
					exit(-1);
				}
			}
			else {
				if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
					printf("Reading VpVs Model lat_c error; exiting!\n");
					exit(-1);
				}
			}
			ptrtmp++;
		}
	}

	if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
		ptrtmp = Dep_c;
		for ( i = 0; i < Ndep_c; i++ ) {
			if ( i < Ndep_c - 1 ) {
				if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
					printf("Reading VpVs Model dep_c error; exiting!\n");
					exit(-1);
				}
			}
			else {
				if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
					printf("Reading VpVs Model dep_c error; exiting!\n");
					exit(-1);
				}
			}
			ptrtmp++;
		}
	}

/* Read P velocity model */
	ptrtmp = Vel_p;
	for ( k = 0; k < Ndep_c; k++ ) {
		for ( j = 0; j < Nlat_c; j++ ) {
			if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
				for ( i = 0; i < Nlon_c; i++ ) {
					if ( i < Nlon_c - 1 ) {
						if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
							printf("Reading VpVs Model Vel_p error; exiting!\n");
							exit(-1);
						}
					}
					else {
						if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
							printf("Reading VpVs Model Vel_p error; exiting!\n");
							exit(-1);
						}
					}
					ptrtmp++;
				}
			}
		}
	}

/* Read S velocity model */
	ptrtmp = Vel_s;
	for ( k = 0; k < Ndep_c; k++ ) {
		for ( j = 0; j < Nlat_c; j++ ) {
			if ( fgets( fracline, sizeof(fracline) - 1, fp ) != NULL ) {
				for ( i = 0; i < Nlon_c; i++ ) {
					if (i < Nlon_c - 1) {
						if ( sscanf( fracline, "%lf %[^\n]", ptrtmp, fracline ) != 2 ) {
							printf("Reading VpVs Model Vel_s error; exiting!\n");
							exit(-1);
						}
					}
					else {
						if ( sscanf( fracline, "%lf", ptrtmp ) != 1 ) {
							printf("Reading VpVs Model Vel_s error; exiting!\n");
							exit(-1);
						}
					}
					ptrtmp++;
				}
			}
		}
	}

	bldmap(bld3, bld4);

/*
	Nx2_c = Nlon_c - 2;
	Nxy2_c = Nx2_c * (Nlat_c-2);
	Nxyz2_c = Nxy2_c * (Ndep_c-2);
*/
	average = 0.0;
	for ( i = 0; i < Nlat_c; i++ ) {
		average += Lat_c[i];
	}
	average /= (double)Nlat_c;
	Ro = earthr(average);
	Rs = Ro + 10.0;

	fclose(fp);
	ptrtmp = NULL;

	return 0;
}

/*
 *
 */
int vel_point_to_grid( const double *vel_p, VEL_GRID *vel_g, const int nxyz, const int nxy, const int nx ) {
	int i;
	const double *velptr;
	const double *velprt_end = vel_p + nxyz;
	VEL_GRID *velgptr;

	for ( i=0; i<nxyz; i++ ) {
		velptr = vel_p + i;
		velgptr = vel_g + i;

		velgptr->vel[0] = *(velptr);
		if ( velptr + nx > velprt_end ) velgptr->vel[2] = *(velptr);
		else velgptr->vel[2] = *(velptr + nx);
		if ( velptr + 1 > velprt_end ) velgptr->vel[4] = *(velptr);
		else velgptr->vel[4] = *(velptr + 1);
		if ( velptr + 1 + nx > velprt_end ) velgptr->vel[6] = *(velptr);
		else velgptr->vel[6] = *(velptr + 1 + nx);

		if ( velptr + nxy > velprt_end ) {
			velgptr->vel[1] = velgptr->vel[0];
			velgptr->vel[3] = velgptr->vel[2];
			velgptr->vel[5] = velgptr->vel[4];
			velgptr->vel[7] = velgptr->vel[6];
		} else {
			velptr += nxy;
			velgptr->vel[1] = *(velptr);
			if ( velptr + nx > velprt_end ) velgptr->vel[3] = *(velptr);
			else velgptr->vel[3] = *(velptr + nx);
			if ( velptr + 1 > velprt_end ) velgptr->vel[5] = *(velptr);
			else velgptr->vel[5] = *(velptr + 1);
			if ( velptr + 1 + nx > velprt_end ) velgptr->vel[7] = *(velptr);
			else velgptr->vel[7] = *(velptr + 1 + nx);
		}
	}

	return 0;
}

/*
 *
 */
static void bldmap( const double bld3, const double bld4 )
{
	int i;
	double lat1, lon1, dep1;
	double lon_now, lat_now, dep_now;

    /* For crustal velocity */
	lon1 = -Lon_c[0];
	lat1 = -Lat_c[0];
	dep1 = -Dep_c[0];

	Ilonmax = (int)(EPS + (Lon_c[Nlon_c-1] + lon1)/bld3);
	Ilatmax = (int)(EPS + (Lat_c[Nlat_c-1] + lat1)/bld3);
	Idepmax = (int)(EPS + (Dep_c[Ndep_c-1] + dep1)/bld4);

	if ( Ilonmax > MSG_NUMBER || Ilatmax > MSG_NUMBER || Idepmax > MSG_NUMBER ) {
		printf("Error; model dimension too big!\n");
		exit(-1);
	}

	Ilon_c = calloc(Ilonmax, sizeof(int));
	Ilat_c = calloc(Ilatmax, sizeof(int));
	Idep_c = calloc(Idepmax, sizeof(int));

	int matrix_index = 0;
	for ( i=0; i<Ilonmax; i++ ) {
		lon_now = (double)i * bld3 - lon1;
		if(lon_now >= Lon_c[matrix_index + 1]) matrix_index++;
		Ilon_c[i] = matrix_index;
	}

	matrix_index = 0;
	for ( i=0; i<Ilatmax; i++ ) {
		lat_now = (double)i * bld3 - lat1;
		if(lat_now >= Lat_c[matrix_index + 1]) matrix_index++;
		Ilat_c[i] = matrix_index;
	}

	matrix_index = 0;
	for ( i=0; i<Idepmax; i++ ) {
		dep_now = (double)i * bld4 - dep1;
		if(dep_now >= Dep_c[matrix_index + 1]) matrix_index++;
		Idep_c[i] = matrix_index;
	}

	/* Change coordinates to integer */
	Ilon1_c = (int)(lon1/bld3);
	Ilat1_c = (int)(lat1/bld3);
	Idep1_c = (int)(dep1/bld4 + EPS);
	Ibld3 = (int)(1/bld3 + EPS);
	Ibld4 = (int)(1/bld4 + EPS);

	//printf("%d %d %d %lf\n", Ilon1_c, Ilat1_c, Idep1_c);

	return;
}

/*
 *
 */
static double velocity( const RAY_INFO ray, const double shiftlon )
{
    double lat, lon, dep;

    lat = geoc_to_geog(90.0 - ray.a * DEG_PER_RAD);
    lon = ray.b * DEG_PER_RAD + shiftlon;
    dep = Ro - ray.r;

    return vel3(lon, lat, dep);
}

/*
 *
 */
static double vel3( const double lon, const double lat, const double dep )
{
	double lonf, latf, depf;
	double wv[8];
	int ip, jp, kp;

	VEL_GRID velg;

	ip = (int)(lon * Ibld3);
	jp = (int)(lat * Ibld3);
	kp = (int)(dep * Ibld4);

	intmap_3d(&ip, &jp, &kp);

	lonf = Lon_c[ip];
	lonf = (lon - lonf) / (Lon_c[ip + 1] - lonf);
	latf = Lat_c[jp];
	latf = (lat - latf) / (Lat_c[jp + 1] - latf);
	depf = Dep_c[kp];
	depf = (dep - depf) / (Dep_c[kp + 1] - depf);

	//memcpy(&velg, Vel_grid + ip + jp*Nx_c + kp*Nxy_c, sizeof(VEL_GRID));
	velg = *(Vel_grid + ip + jp * Nx_c + kp * Nxy_c);

	/*
		lonf1 = 1.0 - lonf
		latf1 = 1.0 - latf
		depf1 = 1.0 - depf
	*/

	wv[6] = wv[7] = lonf * latf;         /* lonf * latf; 1,1,0 & 1,1,1   */
	wv[4] = wv[5] = lonf - wv[7];        /* lonf * latf1; 1,0,0 & 1,0,1  */
	wv[2] = wv[3] = latf - wv[7];        /* lonf1 * latf; 0,1,0 & 0,1,1  */
	wv[0] = wv[1] = 1.0 - lonf - wv[3];  /* lonf1 * latf1; 0,0,0 & 0,0,1 */

	wv[0] *= velg.vel[0];
	wv[2] *= velg.vel[2];
	wv[4] *= velg.vel[4];
	wv[6] *= velg.vel[6];
	wv[0] += wv[2] + wv[4] + wv[6];

	if ( depf <= EPS ) {
		wv[1] = 0.0;
	}
	else {
		wv[1] *= velg.vel[1];
		wv[3] *= velg.vel[3];
		wv[5] *= velg.vel[5];
		wv[7] *= velg.vel[7];

		wv[1] += wv[3] + wv[5] + wv[7];
		wv[1]  = wv[0] - wv[1];
		wv[1] *= depf;

		wv[0] -= wv[1];
	}

	return wv[0];
}

/*
 *
 */
static void intmap_3d( int *ip, int *jp, int *kp )
{
	const int lon = *ip;
	const int lat = *jp;
	const int dep = *kp;

/*
	*ip = (int)((lon + Lon1_c)/Bld3-1.0);
	*jp = (int)((lat + Lat1_c)/Bld3-1.0);
	*kp = (int)((dep + Dep1_c)/Bld4-1.0);
*/

	*ip += Ilon1_c;
	*jp += Ilat1_c;
	*kp += Idep1_c;

	if (
		*ip < 0 || *jp < 0 || *kp < 0 ||
		*ip >= (Ilonmax - 2) || *jp >= (Ilatmax - 2) || *kp >= (Idepmax - 2)
	) {
		printf("Error, lon, lat and dep out of range!\n");
		printf("lon=%lf, lat=%lf, dep=%lf\n", (double)lon/Ibld3, (double)lat/Ibld3, (double)dep/Ibld4);
		printf("ip=%d, jp=%d, kp=%d\n", *ip, *jp, *kp);
		exit(-1);
	}

	*ip = Ilon_c[*ip];
	*jp = Ilat_c[*jp];
	*kp = Idep_c[*kp];

	return;
}

/*
 * This subroutine calculate position of new point
 * in polar coordinates basing on the coordinates
 * of main point in radians (la is colatitude) and dx and dy in kilometers
 */
static void km2deg( const RAY_INFO origin, const double dx, const double dy, RAY_INFO *next )
{
	double dps;

	dps     = origin.r * sin(origin.a);
	next->b = origin.b + atan2(dx, dps);
	next->a = origin.a + atan2(dy, origin.r);

	if ( next->a > HALFCYCLE ) {
		next->a = FULLCYCLE - next->a;
		next->b = HALFCYCLE + next->b;
	}
	if ( next->a < 0.0 ) {
		next->a = fabs(next->a);
		next->b = HALFCYCLE + next->b;
	}

	if ( next->b < 0.0 )
		next->b = FULLCYCLE + next->b;
	if ( next->b > FULLCYCLE )
		next->b = next->b - FULLCYCLE;

	next->r = sqrt(origin.r*origin.r + dx * dx + dy * dy + EPS);

	return;
}

/*
 *
 */
static double rtim( const int m, const RAY_INFO *ray, const NODE_INFO *node )
{
	const RAY_INFO *ray_end = ray + m;
	RAY_INFO ray_now = *ray;

	double dl;
	double r1, r2, sin1, sin2, cos1, cos2, cosa;
	double r1sq, r2sq;
	double rv1, rv2, sm;
	double result = 0.0;

/* */
	if ( m > MSG_NUMBER ) {
		printf("rtim: ERROR! Node number(ni) is larger than %d!\n", MSG_NUMBER);
		return result;
	}

/* */
	rv1  = 0.5/ray_now.v; /* 1.0/v[0]/2.0 */
	r1   = ray_now.r;
	r1sq = node->rsq;
	sin1 = node->sina;
	cos1 = node->cosa;
	cosa = ray_now.b;
/* */
	node++;
	ray++;

/* */
	do {
		ray_now = *ray; /* Cache */
		r2 = ray_now.r;
		if ( ray < ray_end ) {
			r2sq = r2 * r2;
			sin2 = sin(ray_now.a);
			cos2 = cos(ray_now.a);
		}
		else {
			r2sq = node->rsq;
			sin2 = node->sina;
			cos2 = node->cosa;
		}
		cosa = cos(cosa - ray_now.b);
		cosa = sin1 * sin2 * cosa + cos1 * cos2;

		dl = r1sq + r2sq - 2.0 * r1 * r2 * cosa;

		rv2 = 0.5 / ray_now.v; /* 1.0/v[i]/2.0 */
		sm = rv1 + rv2;
		result += sqrt(dl + EPS) * sm;

		rv1  = rv2;
		r1   = r2;
		r1sq = r2sq;
		sin1 = sin2;
		cos1 = cos2;
		cosa = ray_now.b;
	} while ( ++ray <= ray_end );

	return result;
}

/*
 * this routine establishes the short distance conversion factors
 * given the origin of coordinates
 * the rotation angle is converted to radians also
 */
static double earthr( const double xlat )
{
	double dlt1, dxlt;
	const double re = 6378.163;
	const double ell = 298.26;

	dxlt = xlat * RAD_PER_DEG;
/* conversion factor for latitude, SDC. conversion factor is 0.99330647 */
	dlt1 = atan(0.99330647 * tan(dxlt));
/* dlt1 = atan(B2A_SQ * tan(dxlt * RAD_PER_DEG / 60.0)); */
	dlt1 = sin(dlt1);
	dlt1 *= dlt1;

	return re * (1.0 - dlt1 / ell);
}

/*
 *
 */
static inline double geog_to_geoc( const double xla )
{
	return atan(B2A_SQ * tan(RAD_PER_DEG * xla)) * DEG_PER_RAD;
}

/*
 *
 */
static inline double geoc_to_geog( const double xla )
{
	return atan(tan(RAD_PER_DEG * xla) / B2A_SQ) * DEG_PER_RAD;
}
