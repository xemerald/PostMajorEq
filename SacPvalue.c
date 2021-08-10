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
static int readSACHdr( FILE *, struct SAChead * );
static void Swap4ByteOrder( void * );
int picker_wu(
	const float *, const float *, const float *, const int, const double, const int, const int,
	double *, double *, double *, double *
);
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
	struct SAChead sh;

	FILE *fp, *fp_station;

	int i;
	int npts, datalen;
	float *seis, *pseis;       /* input trace buffer */
	float max_seis, max_pseis;
	float baseline;
	float acc0, vel0;
	double sTime;

	STAINFO stainfo;

	float alt;
	float gain[3], gain_factor;
	short type;
	short flag = 0;
	short swapIsNeeded = 0;

	double	pa, sa, pw, sw;



/* Check command line arguments */
	if ( argc != 2 ) {
		fprintf(stderr, "Usage: <Inputfile>\n");
		exit(0);
	}

	memset(&stainfo, 0, sizeof(stainfo));

/* Open the SAC file */
	if ( (fp = fopen(argv[1], "rb")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening %s\n", argv[1]);
		exit(1);
	}

	if ( (swapIsNeeded = readSACHdr(fp, &sh)) < 0 ) {
		fclose(fp);
		exit(1);
	}

	npts = (int)sh.npts;
	datalen = npts * sizeof(float);

	if ( (seis = (float *) malloc((size_t)datalen)) == (float *)NULL ) {
		fprintf(stderr, "Out of memory for %d SAC samples\n", npts);
		fclose(fp);
		exit(1);
	}

	if ( (pseis = (float *) malloc((size_t)datalen)) == (float *)NULL ) {
		fprintf(stderr, "Out of memory for %d SAC samples\n", npts);
		fclose(fp);
		exit(1);
	}

/* Read the sac data into a buffer */
	if ( (i = (int)fread(seis, sizeof(float), sh.npts, fp)) != npts ) {
		fprintf(stderr, "Error reading SAC data: %s\n", strerror(errno));
		fclose(fp);
		exit(1);
	}
	fclose(fp);

	if ( swapIsNeeded == 1 ) {
		for ( i = 0; i < npts; i++ )
			Swap4ByteOrder( &(seis[i]) );
	}

	do {
		if ( sh.kcmpnm[2] == 'Z' ) {
			gain_factor = 0.059814;
			flag = 0;
		}
		else if ( sh.kcmpnm[2] == 'N' ) {
			gain_factor = 0.059814;
			flag = 0;
		}
		else if ( sh.kcmpnm[2] == 'E' ) {
			gain_factor = 0.059814;
			flag = 0;
		}
		else {
			if ( flag > 1 ) {
				printf("Incorrect! Please enter the direction again (Z, N or E):\n");
			}
			else {
				printf("The SAC header do not define the component or not clear at all, we can't define the direction!\n");
				printf("Enter the direction (Z, N or E) for %s.%s.%s:\n", sh.kstnm, sh.kcmpnm, sh.knetwk);
			}
			scanf("%c", &sh.kcmpnm[2]);
			flag = 2;
		}
	} while ( flag );

	flag = 0;
	acc0 = 0.0;
	vel0 = 0.0;
	baseline = 0.0;
/* Find the baseline */
	for ( i = 0; i < npts; i++ ) {
		if ( seis[i] == SACUNDEF ) {
			seis[i] = 0.0;
			flag++;
			continue;
		}
		seis[i]  *= gain_factor;
		baseline += seis[i];
	}
	baseline /= (float)(npts - flag);
	max_seis  = 0.0;
	max_pseis = 0.0;

	for ( i = 0; i < npts; i++ ) {
		if ( seis[i] != 0.0 )
			seis[i] -= baseline;

		if ( fabs(seis[i]) > max_seis )
			max_seis = fabs(seis[i]);
		if ( fabs(pseis[i]) > max_pseis )
			max_pseis = fabs(pseis[i]);
	}

	picker_wu(
		seis, NULL, NULL, npts, sh.delta, 2, 0,
		&pa, &sa, &pw, &sw
	);

	printf("%lf %lf %lf %lf\n", pa, sa, pw, sw);

	free(seis);
	free(pseis);

	return 0;
}


/*
 * readSACHdr: read the header portion of a SAC file into memory.
 *  arguments: file pointer: pointer to an open file from which to read
 *             filename: pathname of open file, for logging.
 * returns: 0 on success
 *          1 on success and if byte swapping is needed
 *         -1 on error reading file
 *     The file is left open in all cases.
 */
static int readSACHdr(FILE *fp, struct SAChead *psh)
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
		fprintf(stderr, "readSACHdr: error reading SAC file: %s\n", strerror(errno));
		return -1;
	}

/* mtheo 2007/10/19
* Guessing if byte swapping is needed
* fileSize should be equal to sizeof(struct SAChead) + (psh->npts * sizeof(float)) */
	if( fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
		ret = 1;
		fprintf(stderr, "WARNING: Swapping is needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
		for ( i=0; i<NUM_FLOAT; i++ )
			Swap4ByteOrder( &(psh2->SACfloat[i]) );
		for ( i=0; i<MAXINT; i++ )
			Swap4ByteOrder( &(psh2->SACint[i]) );
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
static void Swap4ByteOrder( void *data )
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
