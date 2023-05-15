
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
/* */
#include <sachead.h>
#include <sac_proc.h>

/* Internal Function Prototypes */
void usage( char * );

/*
 *
 */
int main(int argc, char **argv)
{
	FILE  *fp_out = NULL;;
	struct SAChead sh;
	char  *SACfile, *outbuf;
	float *seis, *sp;       /* input trace buffer */
	float gain_factor = 1.0;
	float acc0, vel0, ave0;
	double sTime;
	int arg, max_samps;
	int i, npts, datalen;
	char net[3];
	int appendOutput = 0;

	double B0, B1, B2;
	double A1, A2;
	double OUTPUT;
	double X1, X2, Y1, Y2;

/* */
	net[0] = 0;
/* */
	if ( argc < 2 ) {
		usage( argv[0] );
		exit(0);
	}

	arg = 1;
	while ( arg < argc && argv[arg][0] == '-' ) {
		switch ( argv[arg][1] ) {
		case 'N':
			arg++;
			strcpy(net, argv[arg]);
			break;
		case 'n':
			arg++;
			max_samps = atoi(argv[arg]);
			break;
		case 'g':
			arg++;
			gain_factor = atof(argv[arg]);
			break;
		case 'a':
			appendOutput = 1;
			break;
		default:
			usage( argv[0] );
			exit(0);
		}
		arg++;
	}
/* */
	if ( argc - arg == 1 ) {
#ifdef _WINNT
		usage( argv[0] );
		exit( 1 );
#else
		SACfile = argv[arg];
		arg++;
		fp_out = stdout;
#endif
	}
	else if ( argc - arg == 2 ) {
		SACfile = argv[arg];
		arg++;
	}
	else {
		usage( argv[0] );
		exit(0);
	}

	if ( sac_proc_sac_load( SACfile, &sh, &seis ) < 0)
		return -1;

	if ( fp_out == NULL ) {
		fp_out = fopen( argv[arg], appendOutput ? "ab" : "wb" );
		if ( fp_out == (FILE *)NULL ) {
			fprintf(stderr, "%s: error opening %s\n", argv[0], argv[arg]);
			exit(-1);
		}
	}

/* */
	sTime = (double)sh.b + sac_proc_reftime_fetch( &sh );
	fprintf(
		stderr, "%s start at %4.4d,%3.3d,%2.2d:%2.2d:%2.2d.%4.4d %f\n",
		sac_proc_scnl_print( &sh ),
		sh.nzyear, sh.nzjday, sh.nzhour, sh.nzmin, sh.nzsec, sh.nzmsec, sTime
	);
/* */
	npts    = (int)sh.npts;
	datalen = npts * sizeof(float);
	outbuf  = malloc(sizeof(struct SAChead) + datalen);
	memcpy(outbuf, (char *)&sh, sizeof(struct SAChead));
	sp = (float *)(outbuf + sizeof(struct SAChead));
/* */
	sac_proc_data_preprocess( seis, npts, gain_factor );
	memcpy(sp, seis, datalen);
/* */
	if ( fwrite(outbuf, 1, sizeof(struct SAChead) + datalen, fp_out) != sizeof(struct SAChead) + datalen ) {
		fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
		exit(-1);
	}
	fprintf(stderr, "%s SAC data convert finished!\n", sac_proc_scnl_print( &sh ) );
	fclose(fp_out);

	return 0;
}

/*
 *
 */
void usage( char *argv0 )
{
	fprintf(stderr, "Usage: %s [-N NN] [-n max-samples] infile >> outfile\n", argv0);
	fprintf(stderr, "    or %s [-N NN] [-n max-samples] [-a] infile outfile\n", argv0);
	fprintf(stderr, "    -N NN is the network code to use from the cmdline instead of SAC file\n");
	//fprintf(stderr, "    Version: %s\n", VERSION_NUM);

	return;
}
