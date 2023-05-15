#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <errno.h>
/* */
#include <sachead.h>
#include <sac_proc.h>

/*
 *
 */
int main( int argc, char **argv )
{
	struct SAChead sh0, sh1;
	char          *outfilename = NULL;
	float         *seis0, *seis1, *gapseis;           /* input trace buffer */
	double         starttime0, starttime1;  /* input trace buffer */
	double         gaptime;  /* input trace buffer */
	int            size0, size1;
	int            gapsamp = 0;
	int            _sh0samp = 0;
	FILE          *fd;
	char           output_path[512] = { 0 };

/* Check command line arguments */
	if ( argc != 4 ) {
		fprintf(stderr, "Usage: %s <Input SAC File 0> <Input SAC File 1> <Output Path>\n", argv[0]);
		exit(0);
	}

	if ( (size0 = sac_proc_sac_load( argv[1], &sh0, &seis0 )) < 0)
		return -1;
	starttime0 = sac_proc_reftime_fetch( &sh0 );
	fprintf(stderr, "First input SAC file ref. time is %.3f, end at %.3f. Total %d samples with %.3f delta.\n", starttime0, starttime0 + sh0.e, sh0.npts, sh0.delta);


	if ( (size1 = sac_proc_sac_load( argv[2], &sh1, &seis1 )) < 0)
		return -1;
	starttime1 = sac_proc_reftime_fetch( &sh1 );
	fprintf(stderr, "Second input SAC file ref. time is %.3f, end at %.3f. Total %d samples with %.3f delta.\n", starttime1, starttime1 + sh1.e, sh1.npts, sh1.delta);

	if ( fabs(sh0.delta - sh1.delta) > 0.000001 ) {
		fprintf(stderr, "The delta between these two SAC files are different(%f & %f). Just exit!\n", sh0.delta, sh1.delta);
		return -1;
	}

/* */
	gaptime = starttime1 - (starttime0 + sh0.e) + sh0.delta * 0.1;
	gapsamp = gaptime / sh1.delta - 1;
	if ( gapsamp >= 0 ) {
		fprintf(stderr, "Gap is %.3f seconds(total %d samples).\n", gaptime, gapsamp);
		fprintf(stderr, "Filling the gap with %.6f...\n", (double)SACUNDEF);
		gapseis = (float *)calloc(gapsamp, sizeof(float));
		for ( int i = 0; i < gapsamp; i++ )
			gapseis[i] = SACUNDEF;
	/* */
		_sh0samp = sh0.npts;
		sh0.npts += gapsamp + sh1.npts;
		sh0.e += gaptime + sh1.e;

	/* */
		if ( (outfilename = strrchr(argv[1], '/')) )
			outfilename++;
		else
			outfilename = argv[1];

		sprintf(output_path, "%s/%s", argv[3], outfilename);
		if ( (fd = fopen(output_path, "wb")) == (FILE *)NULL ) {
			fprintf(stderr, "Error opening %s\n", output_path);
			return -1;
		}
	/* */
		if ( fwrite(&sh0, 1, sizeof(struct SAChead), fd) != sizeof(struct SAChead) ) {
			fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
			return -1;
		}
	/* */
		if ( fwrite(seis0, sizeof(float), _sh0samp, fd) != _sh0samp ) {
			fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
			return -1;
		}
		if ( fwrite(gapseis, sizeof(float), gapsamp, fd) != gapsamp ) {
			fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
			return -1;
		}
		if ( fwrite(seis1, sizeof(float), sh1.npts, fd) != sh1.npts ) {
			fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
			return -1;
		}
	/* */
		free(gapseis);
		fclose(fd);
		fprintf(stderr, "%s processing finished!\n", output_path);
	}
	else {
		fprintf(stderr, "There is not any gap between these two SAC files. Just exit!\n");
	}
	free(seis0);
	free(seis1);

	return 0;
}
