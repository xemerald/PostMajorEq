/*
 * Standalone program to read SAC data files and compute
 * the peak values include acceleration, velocity & displacement.
 *
 * It can also derive the warning time for each station.
 *
 * Benjamin Yang; Feb 2018
 */

#define VERSION_NUM  "1.0.0 2018-02-14"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "sachead.h"
#include "swap.h"
#include "time_ew.h"

#define DEF_MAX_SAMPS 100

/* Internal Function Prototypes */
static int readSACHdr(FILE *, struct SAChead *);
static double sacRefTime( struct SAChead * );
static int strib( char *string );

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

	double	B0, B1, B2;
	double	A1, A2;
	double	OUTPUT;
	double	X1, X2, Y1, Y2;


	/* Check command line arguments
	******************************/
	if ( argc != 2 )
	{
		fprintf( stderr, "Usage: <Inputfile>\n" );
		exit( 0 );
	}

	memset( &stainfo, 0, sizeof(stainfo) );

	/* Open the SAC file
	********************/
	if ( (fp = fopen( argv[1], "rb" )) == (FILE *)NULL)
	{
		fprintf(stderr, "Error opening %s\n", argv[1]);
		exit( 1 );
	}

	if ( (swapIsNeeded = readSACHdr(fp, &sh)) < 0)
	{
		fclose(fp);
		exit( 1 );
	}

	if ( (fp_station = fopen( "StationList","r" )) == (FILE *)NULL )
	{
		fprintf(stderr, "Error opening Station list\n");
		fclose(fp);
		exit( 1 );
	}
	else
	{
		flag = 0;

		while ( fscanf( fp_station, "%s %f %f %f %f %f %f %hd",
				stainfo.station, &stainfo.latitude, &stainfo.longitude, &alt, &gain[0], &gain[1], &gain[2], &type ) == 8 )
		{
			if ( strncmp(sh.kstnm, stainfo.station, strlen(stainfo.station)) == 0 )
			{
				flag = 1;
				break;
			}
		}

	}

	fclose(fp_station);

	if ( !flag )
	{
		fprintf(stderr, "The station %s is not in the list!\n", sh.kstnm );
		fclose(fp);
		exit( 1 );
	}

	npts = (int)sh.npts;
	datalen = npts * sizeof(float);

	if ( (seis = (float *) malloc((size_t)datalen)) == (float *)NULL )
	{
		fprintf(stderr, "Out of memory for %d SAC samples\n", npts);
		fclose(fp);
		exit( 1 );
	}

	if ( (pseis = (float *) malloc((size_t)datalen)) == (float *)NULL )
	{
		fprintf(stderr, "Out of memory for %d SAC samples\n", npts);
		fclose(fp);
		exit( 1 );
	}

	sTime = (double)sh.b + sacRefTime( &sh );
	fprintf(stdout, "Start %4.4d,%3.3d,%2.2d:%2.2d:%2.2d.%4.4d %f\n", sh.nzyear,
			sh.nzjday, sh.nzhour, sh.nzmin, sh.nzsec, sh.nzmsec, sTime);

	/* Read the sac data into a buffer */
	if ( (i = (int)fread(seis, sizeof(float), sh.npts, fp)) != npts)
	{
		fprintf(stderr, "Error reading SAC data: %s\n", strerror(errno));
		fclose(fp);
		exit( 1 );
	}
	fclose(fp);

	if ( swapIsNeeded == 1 )
	{
		for (i = 0; i < npts; i++) SwapFloat( &(seis[i]) );
	}

	do
	{
		if ( sh.kcmpnm[2] == 'Z' )
		{
			gain_factor = gain[0];
			flag = 0;
		}
		else if ( sh.kcmpnm[2] == 'N' )
		{
			gain_factor = gain[1];
			flag = 0;
		}
		else if ( sh.kcmpnm[2] == 'E' )
		{
			gain_factor = gain[2];
			flag = 0;
		}
		else
		{
			if ( flag > 1 )
			{
				printf("Incorrect! Please enter the direction again (Z, N or E):\n");
			}
			else
			{
				printf("The SAC header do not define the component or not clear at all, we can't define the direction!\n");
				printf("Enter the direction (Z, N or E) for %s.%s.%s:\n", sh.kstnm, sh.kcmpnm, sh.knetwk);
			}

			scanf("%c", &sh.kcmpnm[2]);
			flag = 2;
		}
	} while ( flag );

	switch ( (int)(1.0/sh.delta) ) {
		case 100:
			B0 = 0.9966734;
			B1 = -1.993347;
			B2 = 0.9966734;
			A1 = -1.993336;
			A2 = 0.9933579;
			break;
		case 50:
			B0 = 0.993357897;
			B1 = -1.98671579;
			B2 = 0.993357897;
			A1 = -1.98667157;
			A2 = 0.986759841;
			break;
		case 40:
			B0 = 0.9917042;
			B1 = -1.983408;
			B2 = 0.9917042;
			A1 = -1.983340;
			A2 = 0.9834772;
		case 20:
			B0 = 0.983477175;
			B1 = -1.96695435;
			B2 = 0.983477175;
			A1 = -1.96668136;
			A2 = 0.967227399;
			break;
		default:
			fprintf(stderr, "SAC sample rate not in the list: %f\n", 1./sh.delta);
			exit( 1 );
	}

	flag = 0;
	acc0 = 0.0;
	vel0 = 0.0;
	baseline = 0.0;
	X1 = 0.0;
	X2 = 0.0;
	Y1 = 0.0;
	Y2 = 0.0;

	for ( i = 0; i < npts; i++ ) /* Find the baseline */
	{
		if ( seis[i] == SACUNDEF ) /* Turn the "-12345" to zero */
		{
			seis[i] = 0.0;
			flag++;

			continue;
		}

		seis[i] *= gain_factor;
		baseline += seis[i];
	}

	baseline /= (float)(npts - flag);
	max_seis = 0.0;
	max_pseis = 0.0;

	for ( i = 0; i < npts; i++ )
	{
		if ( seis[i] != 0.0 ) seis[i] -= baseline;

		if ( type == 1 ) /* Acceleration to Velocity */
		{
			pseis[i] = ( seis[i] + acc0 )*sh.delta/2.0 + vel0;
			acc0 = seis[i];
			vel0 = pseis[i];

		/* Recursive Filter, high pass 2 poles at 0.075 Hz */
			OUTPUT = B0*pseis[i] + B1*X1 + B2*X2;
			OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 );
			Y2 = Y1;
			Y1 = OUTPUT;
			X2 = X1;
			X1 = pseis[i];
			pseis[i] = OUTPUT;
		}
		else if ( type == 2 ) /* Velocity to Acceleration */
		{
			pseis[i] = ( seis[i] - vel0 )/sh.delta;
			vel0 = seis[i];
			acc0 = pseis[i];
		}

		if ( fabs(seis[i]) > max_seis ) max_seis = fabs(seis[i]);
		if ( fabs(pseis[i]) > max_pseis ) max_pseis = fabs(pseis[i]);
	}


	fp = fopen( "result_tmp", "ab" );
	if ( fp == (FILE *)NULL )
	{
		fprintf(stderr, "Error opening result file!\n");
		exit( 1 );
	}

	if ( type == 1 )
	{
		stainfo.pga = max_seis;
		stainfo.pgv = max_pseis;
	}
	else if ( type == 2 )
	{
		stainfo.pga = max_pseis;
		stainfo.pgv = max_seis;
	}

	if ( fwrite( &stainfo, 1, sizeof(STAINFO), fp) != sizeof(STAINFO) )
	{
		fprintf(stderr, "Error writing result_tmp file: %s\n", strerror(errno));
		exit( 1 );
	}

	fclose(fp);
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

  // obtain file size:
  fseek (fp , 0 , SEEK_END);
  fileSize = ftell (fp);
  rewind (fp);

  psh2 = (struct SAChead2 *)psh;

  if (fread( psh, sizeof(struct SAChead2), 1, fp) != 1)
  {
    fprintf(stderr, "readSACHdr: error reading SAC file: %s\n",
            strerror(errno));
    return -1;
  }

  /* mtheo 2007/10/19
   * Guessing if byte swapping is needed
   * fileSize should be equal to sizeof(struct SAChead) + (psh->npts * sizeof(float)) */
  if(fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
      ret = 1;
      fprintf(stderr, "WARNING: Swapping is needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
      for (i = 0; i < NUM_FLOAT; i++)
	  SwapFloat( &(psh2->SACfloat[i]) );
      for (i = 0; i < MAXINT; i++)
	  SwapInt( &(psh2->SACint[i]) );
      if(fileSize != (sizeof(struct SAChead) + (psh->npts * sizeof(float))) ) {
	  fprintf(stderr, "ERROR: Swapping is needed again! (fileSize %d, psh.npts %d)\n",
		  fileSize, psh->npts);
	  ret = -1;
      }
  } else {
      fprintf(stderr, "Swapping is not needed! (fileSize %d, psh.npts %d)\n", fileSize, psh->npts);
  }

// /* SAC files are always in "_SPARC" byte order; swap if necessary */
// #ifdef _INTEL
//   for (i = 0; i < NUM_FLOAT; i++)
//     SwapLong( (long *) &(psh2->SACfloat[i]));
//   for (i = 0; i < MAXINT; i++)
//     SwapLong( (long *) &(psh2->SACint[i]));
// #endif

  return ret;
}

/*
 * sacRefTime: return SAC reference time as a double.
 *             Uses a trick of mktime() (called by timegm_ew): if tm_mday
 *             exceeds the normal range for the month, tm_mday and tm_mon
 *             get adjusted to the correct values. So while mktime() ignores
 *             tm_yday, we can still handle the julian day of the SAC header.
 *             This routine does NOT check for undefined values in the
 *             SAC header.
 *  Returns: SAC reference time as a double.
 */
static double sacRefTime( struct SAChead *pSH )
{
  struct tm tms;
  double sec;

  tms.tm_year = pSH->nzyear - 1900;
  tms.tm_mon = 0;    /* Force the month to January */
  tms.tm_mday = pSH->nzjday;  /* tm_mday is 1 - 31; nzjday is 1 - 366 */
  tms.tm_hour = pSH->nzhour;
  tms.tm_min = pSH->nzmin;
  tms.tm_sec = pSH->nzsec;
  tms.tm_isdst = 0;
  sec = (double)timegm_ew(&tms);

  return (sec + (pSH->nzmsec / 1000.0));
}

/*
 * strib: strips trailing blanks (space, tab, newline)
 *    Returns: resulting string length.
 */
static int strib( char *string )
{
  int i, length = 0;

  if ( string && (length = strlen( string )) > 0)
  {
    for ( i = length-1; i >= 0; i-- )
    {
      switch ( string[i] )
      {
		  case ' ':
		  case '\n':
		  case '\t':
		    string[i] = '\0';
		    break;
		  default:
		    return ( i+1 );
      }
    }
  }

  return length;
}
