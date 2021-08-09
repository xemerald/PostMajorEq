#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include <sachead.h>
#include <trace_buf.h>
#include <swap.h>
#include <tb2sac.h>

#define VERSION         "0.0.1 - 2021-08-06"
#define AUTHOR          "Benjamin Ming Yang"

static int  extract2sac( void const *, TRACE_NODE * );

extern volatile unsigned int TotalDuty;
extern volatile unsigned int CompDuty;

static char         OutDir[128]   = { 0 };
static char         OutFormat[16] = { 0 };

static void        *Root      = NULL;         /* Root of binary tree */
static TRACE_NODE **TraceList = NULL;         /* List for remaping */

int main( int argc, char *argv[] )
{
	int    i;
	int    ifd;         /* file of waveform data to read from   */
	struct stat fs;

	void *tankstart;
	void *tankend;
	int   totaltrace = 0;

	struct timespec tt1, tt2;  /* Nanosecond Timer */

	if ( argc != 2 ) {
		fprintf(stderr, "Usage: postmajor <tankfile>\n" );
		fprintf(stderr, " version: %s\n", VERSION );
		fprintf(stderr, " author:  %s\n", AUTHOR );
		return 0;
	}

/* Nanosecond Timer */
	clock_gettime(CLOCK_MONOTONIC, &tt1);

/* Open a waveform files
***********************/
	ifd = open(argv[1], O_RDONLY, 0);
	if ( ifd < 0 ) {
		fprintf(stderr, "%s Cannot open tankfile <%s>\n", nowprog(), argv[1] );
		return -1;
	}

/* For future function... */
	sprintf(OutFormat, "SAC");

/* Create the output directory
******************************/
	sprintf(OutDir, "%s_%s", argv[1], OutFormat);
/* Move the index to the head of real filename, just skip the path */
	if ( strrchr(OutDir, '/') != NULL )
		i = strrchr(OutDir, '/') - OutDir;
	else i = 0;
/* Replace the '.' in the file name to '_' */
	for ( ; i < (int)strlen(OutDir); i++ )
		if ( OutDir[i] == '.' ) OutDir[i] = '_';
/* Check if the directory is existing or not */
	if ( stat(OutDir, &fs) == -1 ) {
		mkdir(OutDir, 0775);
	}

	fstat(ifd, &fs);
	fprintf(stdout, "%s Open the tankfile %s, size is %ld bytes.\n", nowprog(), argv[1], (size_t)fs.st_size);
	fprintf(stdout, "%s Mapping the tankfile %s into memory...\n", nowprog(), argv[1]);
	tankstart = mmap(NULL, (size_t)fs.st_size, PROT_READ | PROT_WRITE, MAP_PRIVATE, ifd, 0);
	tankend   = (uint8_t *)tankstart + (size_t)fs.st_size;

	if ( (totaltrace = scantracebuf( tankstart, tankend )) <= 0 ) {
		fprintf(stderr, "%s Cannot mark all the tracebuf from <%s>\n", nowprog(), argv[1] );
		return -1;
	}

	TotalDuty = totaltrace;
	fprintf(stdout, "%s Estimation complete, total %d traces.\n", nowprog(), totaltrace );

	for ( i = 0; i < totaltrace; i++ ) {
		extract2sac( tankstart, TraceList[i] );
	}

	munmap(tankstart, (size_t)fs.st_size);
	free(TraceList);
	close(ifd);

/* Nanosecond Timer */
	clock_gettime(CLOCK_MONOTONIC, &tt2);
	fprintf(stdout, "%s Convertion complete! Total processing time: %.3f sec.\n", nowprog(),
		(float)(tt2.tv_sec - tt1.tv_sec) + (float)(tt2.tv_nsec - tt1.tv_nsec)* 1e-9);

	return 0;
}

/*
*/
static int extract2sac( void const *tankstart, TRACE_NODE *tnode )
{
	int    i, j;
	int    gapcount  = 0;
	int    nfill_max = 0;
	int    nsamp_trace;
	double starttime, endtime; /* times for current scnl         */
	double samprate = 1.0;

	TRACE2_HEADER  *trh2     = NULL;  /* Tracebuf message read from file      */
	uint8_t        *tankbyte = NULL;

	uint8_t        *outbuf       = NULL;
	uint8_t        *outbufend    = NULL;
	size_t          buffersiz    = sizeof(struct SAChead) + MAX_NUM_TBUF * 100 * sizeof(SACWORD);
	char            sacfile[256] = { 0 };
	struct SAChead *sachead      = NULL;
	SACWORD        *seis         = NULL;
	SACWORD         multiplier   = 1.0;
	SACWORD         fill         = (SACWORD)SACUNDEF;

	fprintf(stdout, "%s Extracting %s.%s.%s.%s to SAC file...\n",
		nowprog(), tnode->sta, tnode->chan, tnode->net, tnode->loc );

	outbuf    = (uint8_t *)malloc(buffersiz);
	outbufend = outbuf + buffersiz;
	sachead   = (struct SAChead *)outbuf;
	seis      = (SACWORD *)(sachead + 1);

	memset(outbuf, 0, buffersiz);

/* Initialize all the columns in the SAC header */
	sacinit( sachead );
/* Set some columns to the default values that every one should be the same */
	sacdefault( sachead );

/* Copy the SCNL into the header and blank the trailing chars */
/* Station name */
	strcpy(sachead->kstnm, tnode->sta);
	for ( i = (int)strlen(tnode->sta); i < K_LEN; i++ ) sachead->kstnm[i] = ' ';
/* Channel code */
	strcpy(sachead->kcmpnm, tnode->chan);
	for ( i = (int)strlen(tnode->chan); i < K_LEN; i++ ) sachead->kcmpnm[i] = ' ';
/* Network code */
	strcpy(sachead->knetwk, tnode->net);
	for ( i = (int)strlen(tnode->net); i < K_LEN; i++ ) sachead->knetwk[i] = ' ';
/* Location code */
	strcpy(sachead->khole, tnode->loc);
	for ( i = (int)strlen(tnode->loc); i < K_LEN; i++ ) sachead->khole[i] = ' ';

/* Orientation of seismometer -	determine the orientation based on the third character
of the component name */
	switch ( tnode->chan[2] ) {
	/* Vertical component */
		case 'Z' :
		case 'z' :
			sachead->cmpaz  = 0.0;
			sachead->cmpinc = 0.0;
			break;
	/* North-south component */
		case 'N' :
		case 'n' :
			sachead->cmpaz  = 0.0;
			sachead->cmpinc = 90.0;
			break;
	/* East-west component */
		case 'E' :
		case 'e' :
			sachead->cmpaz  = 90.0;
			sachead->cmpinc = 90.0;
			break;
	/* Anything else */
		default :
			sachead->cmpaz  = (float)SACUNDEF;
			sachead->cmpinc = (float)SACUNDEF;
			break;
	} /* switch */

	fprintf(stdout, "%s SAC header of %s.%s.%s.%s preparation complete!\n",
		nowprog(), tnode->sta, tnode->chan, tnode->net, tnode->loc );

	nsamp_trace = 0;
	endtime     = time(NULL);

/* Go through the tracebuf list of this trace */
	for ( i = 0; i < tnode->ntbuf; i++ ) {
		tankbyte  = (uint8_t *)tankstart + tnode->tlist[i].offset;
		trh2      = (TRACE2_HEADER *)tankbyte;
		starttime = trh2->starttime;
		samprate  = trh2->samprate;

		int32_t   *idata           = (int32_t *)(trh2 + 1);
		int16_t   *sdata           = (int16_t *)idata;
		float     *fdata           = (float *)idata;
		double    *ddata           = (double *)idata;
		const char byte_order      = trh2->datatype[0];         /* Byte order of this TYPE_TRACEBUF2 msg */
		const int  byte_per_sample = atoi(&trh2->datatype[1]);  /* For TYPE_TRACEBUF2 msg                */

		if ( i > 0 ) {
		/* Starttime is set for new packet; endtime is still set for old packet */
			if ( endtime + (2.0 / samprate) < starttime ) {
			/* There's a gap, so fill it */
				int    nfill   = (int)((float)samprate * (float)(starttime - endtime));
				size_t newsize = (nsamp_trace + nfill)*sizeof(SACWORD) + sizeof(struct SAChead);

				if ( newsize > buffersiz ) {
					if ( newsize > buffersiz*2 ) {
						fprintf(stderr, "%s *** Bogus gap (%d); skipping! ***\n", nowprog(), nfill);
						continue;
					}
					else {
					/* Still need to be revised, 'cause this block is repeated */
						buffersiz += MAX_NUM_TBUF * 100 * sizeof(SACWORD);
						outbuf = realloc(outbuf, buffersiz);
						if ( outbuf == NULL ) {
							fprintf(stderr, "%s *** Could not realloc output buffer to %ld bytes ***\n",
								nowprog(), buffersiz);
							return -1;
						}
					/* */
						seis      = (SACWORD *)((struct SAChead *)outbuf + 1) + (seis - (SACWORD *)(sachead + 1));
						sachead   = (struct SAChead *)outbuf;
						outbufend = outbuf + buffersiz;
					}
				}
			/* Do the filling */
				for ( j = 0; j < nfill && seis < (SACWORD *)outbufend; j++, seis++ ) *seis = fill;
				nsamp_trace += nfill;
			/* Keep track of how many gaps and the largest one */
				gapcount++;
				if ( nfill_max < nfill ) nfill_max = nfill;
			}
			else if ( endtime + (1.0 / samprate) > trh2->endtime ) {
			/* This is a duplicated tracebuf, we'll just skip it */
				continue;
			}
		}
		else {
		/* gmttime makes juldays starting with 0 */
			time_t ltime    = (time_t)starttime;
			struct tm *ltm  = gmtime(&ltime);

			sachead->nzyear = (int32_t)ltm->tm_year + (int32_t)1900; /* Calendar year of reference time */
			sachead->nzjday = (int32_t)ltm->tm_yday + (int32_t)1;    /* Julian day, 0 - 365 */
			sachead->nzhour = (int32_t)ltm->tm_hour;
			sachead->nzmin  = (int32_t)ltm->tm_min;
			sachead->nzsec  = (int32_t)ltm->tm_sec;
			sachead->nzmsec = (int32_t)((starttime - (int32_t)starttime) * 1000.0);
		}

	/* Depends on the data type chose the pointer */
		if ( byte_order == 'i' || byte_order == 's' ) {
			if ( byte_per_sample == 4 ) {
			/* Integer 4 bytes */
				for ( j = 0; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, idata++ )
					*seis = *idata * multiplier;
			}
			else if ( byte_per_sample == 2 ) {
			/* Short integer 2 bytes */
				for ( j = 0; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, sdata++ )
					*seis = *sdata * multiplier;
			}
			else continue;
		}
		else if ( byte_order == 'f' || byte_order == 't' ) {
		/* Following would not work now... */
			if ( byte_per_sample == 4 ) {
			/* Float 4 bytes */
				for ( j = 0; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, fdata++ )
					*seis = *fdata * multiplier;
			}
			else if ( byte_per_sample == 8 ) {
			/* Double float 8 bytes */
				for ( j = 0; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, ddata++ )
					*seis = *ddata * multiplier;
			}
			else continue;
		}
		else continue;

	/* Allocate more space if necessary
	**********************************/
		if ( seis >= (SACWORD *)outbufend ) {
			buffersiz += MAX_NUM_TBUF * 100 * sizeof(SACWORD);
			outbuf = realloc(outbuf, buffersiz);
			if ( outbuf == NULL ) {
				fprintf(stderr, "%s *** Could not realloc output buffer to %ld bytes ***\n",
					nowprog(), buffersiz);
				return -1;
			}
		/* */
			seis      = (SACWORD *)((struct SAChead *)outbuf + 1) + (seis - (SACWORD *)(sachead + 1));
			sachead   = (struct SAChead *)outbuf;
			outbufend = outbuf + buffersiz;

			if ( byte_order == 'i' || byte_order == 's' ) {
				if ( byte_per_sample == 4 ) {
				/* Integer 4 bytes */
					for ( ; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, idata++ )
						*seis = *idata * multiplier;
				}
				else if ( byte_per_sample == 2 ) {
				/* Short integer 2 bytes */
					for ( ; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, sdata++ )
						*seis = *sdata * multiplier;
				}
			}
			else if ( byte_order == 'f' || byte_order == 't' ) {
			/* Following would not work now... */
				if ( byte_per_sample == 4 ) {
				/* Float 4 bytes */
					for ( ; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, fdata++ )
						*seis = *fdata * multiplier;
				}
				else if ( byte_per_sample == 8 ) {
				/* Double float 8 bytes */
					for ( ; j < trh2->nsamp && seis < (SACWORD *)outbufend; j++, seis++, ddata++ )
						*seis = *ddata * multiplier;
				}
			}
		}
	/* Advance endtime to the new packet; process this packet in the next iteration */
		endtime = trh2->endtime;
		nsamp_trace += trh2->nsamp;
	}

/*  All trace data fed into SAC data section.  Now fill in the rest of header */
	sachead->npts  = (int32_t)nsamp_trace;                 /* Samples in trace */
	sachead->delta = (float)(1.0/samprate);                /* Sample period */
	sachead->e     = (float)nsamp_trace * sachead->delta;  /* End time */

/* Output to the SAC file...
****************************/
	sprintf(sacfile, "%s/%s_%s_%s_%s.sac", OutDir, tnode->sta, tnode->chan, tnode->net, tnode->loc);
/* Compute the total size of the SAC file */
	size_t totalbyte = sizeof(struct SAChead) + sachead->npts * sizeof(SACWORD);
/* Open the file and write all buffer data into the file */
	FILE *ofp = fopen(sacfile, "wb");
	if ( fwrite(outbuf, 1, totalbyte, ofp) != totalbyte ) {
		fprintf(stderr, "%s *** Error writing sacfile: %s ***\n", nowprog(), strerror(errno));
		return -1;;
	}
	fclose(ofp);
	free(outbuf);

/* If there is any gap in this trace try to notice the user */
	if ( gapcount ) {
		fprintf(stderr, "%s *** %d gaps; largest %d for <%s.%s.%s.%s> ***\n",
			nowprog(), gapcount, nfill_max, tnode->sta, tnode->chan, tnode->net, tnode->loc);
	}

/* */
	CompDuty++;
	fprintf(stdout, "%s Finish the extraction of %s.%s.%s.%s to SAC file!\n",
		nowprog(), tnode->sta, tnode->chan, tnode->net, tnode->loc );

	return 0;
}
