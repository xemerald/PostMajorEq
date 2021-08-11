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
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "sachead.h"

#define DEF_MAX_SAMPS 100


/* Internal Function Prototypes */
static int read_sac_header( FILE *, struct SAChead * );
static void swap_order_4byte( void * );

int picker_wu(
	const float *, const float *, const float *, const int, const double, const int, const int,
	double *, double *, double *, double *
);

/*
 *
 */
int main( int argc, char **argv )
{
	struct SAChead sh;

	FILE  *fp;
	float *seis[3];       /* input trace buffer */

	int   i, j;
	int   npts, datalen;
	int   _flag = 0;
	float baseline;
	float gain_factor;

/* Final output data */
	float pga, pgv, pgd;
	float pga_time, pgv_time, pgd_time;
	float pd35_time, pga80_time, pga4_time, pga2_time;
	float pga_leadtime, pgv_leadtime;
	float pd3, tc;
	double p_arrival, s_arrival, p_weight, s_weight;

/* Check command line arguments */
	if ( argc != 4 ) {
		fprintf(stderr, "Usage: %s <Z Component File> <N Component File> <E Component File>\n", argv[0]);
		exit(0);
	}

/* Open all the three SAC files */
	for ( i = 0; i < 3; i++ ) {
	/* Opening the SAC files */
		if ( (fp = fopen(argv[i+1], "rb")) == (FILE *)NULL ) {
			fprintf(stderr, "Error opening %s\n", argv[i+1]);
			exit(1);
		}
		if ( (_flag = read_sac_header(fp, &sh)) < 0 ) {
			fclose(fp);
			exit(1);
		}
	/* */
		npts    = (int)sh.npts;
		datalen = npts * sizeof(float);
		if ( (seis[i] = (float *)malloc((size_t)datalen)) == (float *)NULL ) {
			fprintf(stderr, "Out of memory for %d float samples\n", npts);
			fclose(fp);
			exit(1);
		}
	/* Read the sac data into a buffer */
		if ( (datalen = (int)fread(seis[i], sizeof(float), sh.npts, fp)) != npts ) {
			fprintf(stderr, "Error reading SAC data: %s\n", strerror(errno));
			fclose(fp);
			exit(1);
		}
		fclose(fp);
	/* */
		if ( _flag == 1 )
			for ( j = 0; j < npts; j++ )
				swap_order_4byte( &(seis[i][j]) );

	/* */
		switch ( sh.kcmpnm[2] ) {
		case 'Z':
			gain_factor = 0.059814;
			break;
		case 'N':
			gain_factor = 0.059814;
			break;
		case 'E':
			gain_factor = 0.059814;
			break;
		default:
			//fprintf(stderr, "The SAC header do not define the component or not clear at all, we can't define the direction!\n");
			gain_factor = 0.059814;
			break;
		}
	/* */
		_flag    = 0;
		baseline = 0.0;
		for ( j = 0; j < npts; j++ ) {
			if ( seis[i][j] == SACUNDEF ) {
				seis[i][j] = 0.0;
				_flag++;
				continue;
			}
			seis[i][j] *= gain_factor;
			baseline   += seis[i][j];
		}
		baseline /= (float)(npts - _flag);
	/* */
		for ( j = 0; j < npts; j++ )
			if ( seis[i][j] != 0.0 )
				seis[i][j] -= baseline;
	}

/*
 *
 */
	_flag = picker_wu(
		seis[0], seis[1], seis[2], npts, sh.delta, 2, 0,
		&p_arrival, &s_arrival, &p_weight, &s_weight
	);

	if ( _flag == 0 ) {
		fprintf(stderr, "Can't find the P arrival time, just skip this station.\n");
		exit(0);
	}

/*
 *
 */
	datalen    = (int)(p_arrival / sh.delta);
	pga        = 0.0;
	pga4_time  = npts * sh.delta;
	pga80_time = npts * sh.delta;
	for ( i = 0; i < 3; i++ ) {
		for ( j = datalen; j < npts; j++ ) {
		/* */
			if ( fabs(seis[i][j]) > 4.0 ) {
			/* */
				if ( (j * sh.delta) < pga4_time )
					pga4_time = j * sh.delta;
			/* */
				if ( fabs(seis[i][j]) > 80.0 ) {
					if ( (j * sh.delta) < pga80_time )
						pga80_time = j * sh.delta;
				}
			}
		/* */
			if ( fabs(seis[i][j]) > pga ) {
				pga      = fabs(seis[i][j]);
				pga_time = j * sh.delta;
			}
		}
	}

/*
 *
 */
	for ( i = 0; i < 3; i++ ) {
		float pseis      = 0.0;
		float seis_prev  = 0.0;
		float pseis_prev = 0.0;
		for ( j = 0; j < npts; j++ ) {
			pseis      = (seis[i][j] + seis_prev) * sh.delta * 0.5 + pseis_prev;
			seis_prev  = seis[i][j];
			pseis_prev = pseis;
		/* */
			seis[i][j] = 0.0;
		}
	}

/*
 *
 */
	pgv = 0.0;
	for ( i = 0; i < 3; i++ ) {
		for ( j = datalen; j < npts; j++ ) {
		/* */
			if ( fabs(seis[i][j]) > pgv ) {
				pgv      = fabs(seis[i][j]);
				pgv_time = j * sh.delta;
			}
		}
	}

/*
 *
 */
	for ( i = 0; i < 3; i++ ) {
		float pseis      = 0.0;
		float seis_prev  = 0.0;
		float pseis_prev = 0.0;
		for ( j = 0; j < npts; j++ ) {
			pseis      = (seis[i][j] + seis_prev) * sh.delta * 0.5 + pseis_prev;
			seis_prev  = seis[i][j];
			pseis_prev = pseis;
		/* */
			seis[i][j] = 0.0;
		}
	}

/*
 *
 */
	pgd       = 0.0;
	pd35_time = npts * sh.delta;
	for ( i = 0; i < 3; i++ ) {
		for ( j = datalen; j < npts; j++ ) {
		/* */
			if ( fabs(seis[i][j]) > 0.35 ) {
				if ( (j * sh.delta) < pd35_time )
					pd35_time = j * sh.delta;
			}
		/* */
			if ( fabs(seis[i][j]) > pgd ) {
				pgd      = fabs(seis[i][j]);
				pgd_time = j * sh.delta;
			}
		}
	}

	pd3 = 0.0;
	for ( j = datalen; j < datalen + 300; j++ ) {
	/* */
		if ( fabs(seis[i][j]) > pd3 )
			pd3 = fabs(seis[i][j]);
	}

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

	return 0;
}


/*
 * read_sac_header: read the header portion of a SAC file into memory.
 *  arguments: file pointer: pointer to an open file from which to read
 *             filename: pathname of open file, for logging.
 * returns: 0 on success
 *          1 on success and if byte swapping is needed
 *         -1 on error reading file
 *     The file is left open in all cases.
 */
static int read_sac_header(FILE *fp, struct SAChead *psh)
{
	int i;
	struct SAChead2 *psh2;
	int fileSize;
	int ret = 0;

/* obtain file size */
	fseek(fp, 0, SEEK_END);
	fileSize = ftell(fp);
	rewind(fp);

	psh2 = (struct SAChead2 *)psh;

	if ( fread(psh, sizeof(struct SAChead2), 1, fp) != 1 ) {
		fprintf(stderr, "read_sac_header: error reading SAC file: %s\n", strerror(errno));
		return -1;
	}

/* mtheo 2007/10/19
* Guessing if byte swapping is needed
* fileSize should be equal to sizeof(struct SAChead) + (psh->npts * sizeof(float)) */
	if( fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
		ret = 1;
		fprintf(stderr, "WARNING: Swapping is needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
		for ( i=0; i<NUM_FLOAT; i++ )
			swap_order_4byte( &(psh2->SACfloat[i]) );
		for ( i=0; i<MAXINT; i++ )
			swap_order_4byte( &(psh2->SACint[i]) );
		if( fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
			fprintf(stderr, "ERROR: Swapping is needed again! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
			ret = -1;
		}
	} else {
		fprintf(stderr, "Swapping is not needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
	}

	return ret;
}

/* Do byte swapping on the given 4-byte integer or float. */
static void swap_order_4byte( void *data )
{
   uint8_t temp;

   union {
      uint8_t c[4];
   } dat;

   memcpy( &dat, data, sizeof(uint32_t) );
   temp     = dat.c[0];
   dat.c[0] = dat.c[3];
   dat.c[3] = temp;
   temp     = dat.c[1];
   dat.c[1] = dat.c[2];
   dat.c[2] = temp;
   memcpy( data, &dat, sizeof(uint32_t) );
   return;
}
