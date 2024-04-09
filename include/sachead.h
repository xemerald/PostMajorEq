
/*
 *   THIS FILE IS UNDER RCS - DO NOT MODIFY UNLESS YOU HAVE
 *   CHECKED IT OUT USING THE COMMAND CHECKOUT.
 *
 *    $Id: sachead.h 1625 2004-07-16 20:30:33Z lombard $
 *
 *    Revision history:
 *     $Log$
 *     Revision 1.11  2004/07/16 20:30:33  lombard
 *     Modified to support SCNL (includes SEED location code)
 *
 *     Revision 1.10  2002/02/19 18:23:40  lucky
 *     Added blessed valus for IDEP variable
 *
 *     Revision 1.9  2001/11/26 22:27:50  lucky
 *     Added SAC_IO #define to match what localmag/gmew are using
 *
 *     Revision 1.8  2001/04/12 03:03:08  lombard
 *     reorgainized to put the format-specific header items in *head.h
 *     and the putaway-specific items in *putaway.h, of which there
 *     currently is only sacputaway.h
 *     Function protoypes used by putaway.c moved to pa_subs.h
 *
 *     Revision 1.7  2000/12/04 20:11:23  lucky
 *     Added cOnset to SAC_ArrivalStruct
 *
 *     Revision 1.6  2000/11/15 17:03:13  lucky
 *     Added SACPABase_next_ev_review
 *
 *     Revision 1.5  2000/08/28 15:39:19  lucky
 *     *** empty log message ***
 *
 *     Revision 1.4  2000/07/24 18:24:22  lucky
 *     Removed MAX_PHS_PER_EQ -- it is set globally in earthworm.h
 *
 *     Revision 1.3  2000/04/26 20:34:09  lucky
 *     Added comments to known fields in sac header
 *
 *     Revision 1.2  2000/03/09 21:58:09  davidk
 *     removed some unneccessary #defines.
 *      Added a lot of new #define constants, typedefs, and functions prototypes
 *     for the SACPA_XXX routines and the SACPABase_XXX routines.
 *
 *     Revision 1.1  2000/02/14 20:05:54  lucky
 *     Initial revision
 *
 *
 */

#pragma once

/*
 *      SAC Header Include File
 *      This file contains structures and definitions for SAC file headers.
 *      If you are looking for SAC putaway routines, look in sacputaway.h
 *
 *      This include file should be used when ever the SAC header
 *      is accessed.
 */

/************************************/
/** #define CONSTANTS              **/
/************************************/

#define SACHEADERSIZE   632     /* # of bytes in SAC header */

#define SACVERSION      6       /* version number of SAC */
#define SAC_I_D         2.0     /* SAC I.D. number */

#define SACWORD         float   /* SAC data is an array of floats */
#define SACUNDEF        -12345  /* undefined value */
#define SACSTRUNDEF     "-12345  "      /* undefined value for strings */

#define SACTEMP4        11

#define SACTEMP6        3
#define SACTEMP7        10

#define KEVNMLEN        16      /* length of event name */
#define K_LEN           8       /* length of all other string fields */

/* header value  constants */
#define SAC_IUNKN       5l
#define SAC_IDISP       6l
#define SAC_IVEL        7l
#define SAC_IACC        8l
#define SAC_IVOLTS      50l
#define SAC_IBEGINTIME  9l
#define SAC_ITIME       1l
#define SAC_IO         11l


#define NUM_FLOAT       70      /* number of floats in header */
#define MAXINT          40      /* number of ints in header */
#define MAXSTRING       24      /* number of strings in header */

/************************************/
/** TYPE DECLARATIONS              **/
/************************************/

/*      This is the actual data structure which is used to enter the header
        variables. */

struct SAChead {
	/* floating point fields */
	float delta;    /* nominal increment between evenly spaced data (sec) */
	float depmin;   /* min value of trace data (dep var) */
	float depmax;   /* max value of trace data (dep var) */
	float scale;    /* scale factor - NOT USED */
	float odelta;   /* Observed delta if different from delta */
	float b;        /* beginning value of time (indep var)  */
	float e;        /* ending value of time (indep var) */
	float o;        /* event origin time, relative to reference time */
	float a;        /* first arrival time, relative to reference time */
	float internal1;
	float t0;       /* user-defined times */
	float t1;
	float t2;
	float t3;
	float t4;
	float t5;
	float t6;
	float t7;
	float t8;
	float t9;
	float f;        /* End of event time, relative to reference time */
	float resp0;
	float resp1;
	float resp2;
	float resp3;
	float resp4;
	float resp5;
	float resp6;
	float resp7;
	float resp8;
	float resp9;
	float stla;     /* station latitude (deg north) */
	float stlo;     /* station intitude (deg east) */
	float stel;     /* station elevation (meters) */
	float stdp;
	float evla;     /* event location, latitude (deg north) */
	float evlo;     /* event location, intitude (deg east) */
	float evel;
	float evdp;     /* event depth (km) */
	float blank1;
	float user0;    /* user defined values */
	float user1;
	float user2;
	float user3;
	float user4;
	float user5;
	float user6;
	float user7;
	float user8;
	float user9;
	float dist;     /* event-station epicentral distance (km) */
	float az;       /* event to station azimuth (deg) */
	float baz;      /* event to station back-azimuth (deg) */
	float gcarc;    /* event to station arc distance (deg) */
	float internal2;
	float internal3;
	float depmen;
	float cmpaz;    /* component azimuth (deg) */
	float cmpinc;   /* component inclination (deg) */
	float blank4[SACTEMP4];

	/* integer fields */
	int nzyear;    /* Reference time = Trace beginning time (year) */
	int nzjday;    /*  (julian day) */
	int nzhour;    /*  (hour) */
	int nzmin;     /*  (minute) */
	int nzsec;     /*  (second) */
	int nzmsec;    /*  (millisecond) */
	int internal4;
	int internal5;
	int internal6;
	int npts;       /* number of points in trace */
	int internal7;
	int internal8;
	int blank6[SACTEMP6];
	int iftype;      /* Type of data:  1 for time series */
	int idep;        /* Type of dependent data =UNKNOWN */
	int iztype;      /* zero time equivalence  =1 for beginning */
	int iblank6a;
	int iinst;
	int istreg;
	int ievreg;
	int ievtyp;      /* event type; IUNKN */
	int iqual;
	int isynth;
	int blank7[SACTEMP7];
	unsigned int leven;    /* =1 for evenly spaced data */
	unsigned int lpspol;   /* =1 for correct polarity, 0 for reversed */
	unsigned int lovrok;
	unsigned int lcalda;
	unsigned int lblank1;

	/* character string fields */
	char kstnm[K_LEN];      /* station name (blank padded) */
	char kevnm[KEVNMLEN];   /* event name */
	char khole[K_LEN];      /* for SEED location code, per Peter Goldstein */
	char ko[K_LEN];         /* Origin time label */
	char ka[K_LEN];         /* First arrival time label */
	char kt0[K_LEN];        /* labels for user-defined times */
	char kt1[K_LEN];
	char kt2[K_LEN];
	char kt3[K_LEN];
	char kt4[K_LEN];
	char kt5[K_LEN];
	char kt6[K_LEN];
	char kt7[K_LEN];
	char kt8[K_LEN];
	char kt9[K_LEN];
	char kf[K_LEN];
	char kuser0[K_LEN];     /* labels for user-defined quantities */
	char kuser1[K_LEN];
	char kuser2[K_LEN];
	char kcmpnm[K_LEN];
	char knetwk[K_LEN];
	char kdatrd[K_LEN];
	char kinst[K_LEN];
};

/* This is a structure of the same size as 'SAChead', but in a form
   that is easier to initialize. */


struct SAChead2 {
	float SACfloat[NUM_FLOAT];
	int SACint[MAXINT];
	char SACstring[MAXSTRING][K_LEN];
};

/* eof sachead.h */
