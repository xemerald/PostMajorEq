/*
 *
 */

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
	struct SAChead sh;
	float         *seis;       /* input trace buffer */
	int            size = 0;
	FILE          *fd;
	char           output_path[512] = { 0 };

/* Check command line arguments */
	if ( argc != 4 ) {
		fprintf(stderr, "Usage: %s <Input SAC File> <New Channel Code> <Output Path>\n", argv[0]);
		exit(0);
	}

	if ( (size = sac_proc_sac_load( argv[1], &sh, &seis )) < 0)
		return -1;

	sac_proc_scnl_modify( &sh, NULL, argv[2], NULL, NULL );

/* */
	sprintf(output_path, "%s/%s", argv[3], sac_proc_scnl_print( &sh ));
	if ( (fd = fopen(output_path, "wb")) == (FILE *)NULL ) {
		fprintf(stderr, "Error opening %s\n", output_path);
		return -1;
	}
/* */
	if ( fwrite(&sh, 1, sizeof(struct SAChead), fd) != sizeof(struct SAChead) ) {
		fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
		return -1;
	}
/* */
	size -= sizeof(struct SAChead);
	if ( fwrite(seis, 1, size, fd) != size ) {
		fprintf(stderr, "Error writing sacfile: %s\n", strerror(errno));
		return -1;
	}
	fclose(fd);
	fprintf(stderr, "%s processing finished!\n", output_path);

	return 0;
}
