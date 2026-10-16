#ifndef GPSSIM_H
#define GPSSIM_H

#include <stdio.h>

//#define FLOAT_CARR_PHASE // For RKT simulation. Higher computational load, but smoother carrier phase.

#define TRUE	(1)
#define FALSE	(0)

/*! \brief Maximum length of a line in a text file (RINEX, motion) */
#define MAX_CHAR (100)

#define GPS_MAX_PRN (32)
#define QZSS_PRN_START (193)
#define QZSS_PRN_END (202)
#define QZSS_MAX_PRN (10)

/*! \brief Maximum number of satellites in RINEX file (GPS + QZSS) */
#define MAX_SAT (GPS_MAX_PRN + QZSS_MAX_PRN)

static inline int prn2sat(int prn)
{
	if (prn >= 1 && prn <= GPS_MAX_PRN)
		return prn - 1;
	if (prn >= QZSS_PRN_START && prn <= QZSS_PRN_END)
		return GPS_MAX_PRN + (prn - QZSS_PRN_START);
	return -1;
}

static inline int sat2prn(int sv)
{
	if (sv >= 0 && sv < GPS_MAX_PRN)
		return sv + 1;
	if (sv >= GPS_MAX_PRN && sv < MAX_SAT)
		return QZSS_PRN_START + (sv - GPS_MAX_PRN);
	return 0;
}

/*! \brief Maximum number of channels we simulate */
#define MAX_CHAN (16)

/*! \brief Maximum number of user motion points */
#ifndef USER_MOTION_SIZE
#define USER_MOTION_SIZE (3000) // max duration at 10Hz
#endif

/*! \brief Maximum duration for static mode*/
#define STATIC_MAX_DURATION (86400) // second

/*! \brief Number of subframes */
#define N_SBF (5) // 5 subframes per frame

/*! \brief Number of words per subframe */
#define N_DWRD_SBF (10) // 10 word per subframe

/*! \brief Number of words */
#define N_DWRD ((N_SBF+1)*N_DWRD_SBF) // Subframe word buffer size

/*! \brief C/A code sequence length */
#define CA_SEQ_LEN (1023)

#define SECONDS_IN_WEEK 604800.0
#define SECONDS_IN_HALF_WEEK 302400.0
#define SECONDS_IN_DAY 86400.0
#define SECONDS_IN_HOUR 3600.0
#define SECONDS_IN_MINUTE 60.0

#define POW2_M5  0.03125
#define POW2_M19 1.907348632812500e-6
#define POW2_M29 1.862645149230957e-9
#define POW2_M31 4.656612873077393e-10
#define POW2_M33 1.164153218269348e-10
#define POW2_M43 1.136868377216160e-13
#define POW2_M55 2.775557561562891e-17

#define POW2_M50 8.881784197001252e-016
#define POW2_M30 9.313225746154785e-010
#define POW2_M27 7.450580596923828e-009
#define POW2_M24 5.960464477539063e-008

// Conventional values employed in GPS ephemeris model (ICD-GPS-200)
#define GM_EARTH 3.986005e14
#define OMEGA_EARTH 7.2921151467e-5
#define PI 3.1415926535898

#define WGS84_RADIUS	6378137.0
#define WGS84_ECCENTRICITY 0.0818191908426

#define R2D 57.2957795131

#define SPEED_OF_LIGHT 2.99792458e8

/* ========================================================================= */
/* GNSS Band Definitions & Parameters                                        */
/* ========================================================================= */

typedef enum {
	BAND_L1 = 0,
	BAND_L5 = 1,
	BAND_DUAL = 2
} gnss_band_t;

/* Carrier Frequencies */
#define CARR_FREQ_L1 (1575.42e6)
#define CARR_FREQ_L5 (1176.45e6)
#define CARR_FREQ    CARR_FREQ_L1 /* Backward compatibility */

/* Chipping Rates */
#define CODE_FREQ_L1 (1.023e6)
#define CODE_FREQ_L5 (10.23e6)
#define CODE_FREQ    CODE_FREQ_L1 /* Backward compatibility */

/* Carrier to Code Ratios */
#define CARR_TO_CODE_L1 (1.0 / 1540.0)
#define CARR_TO_CODE_L5 (1.0 / 115.0)
#define CARR_TO_CODE    CARR_TO_CODE_L1 /* Backward compatibility */

/* Code Sequence Lengths (per 1 ms) */
#define CA_SEQ_LEN_L1 (1023)
#define CA_SEQ_LEN_L5 (10230)
#define MAX_CODE_SEQ_LEN (10230)

/* Carrier Wavelengths */
#define LAMBDA_L1 0.190293672798365
#define LAMBDA_L5 (SPEED_OF_LIGHT / CARR_FREQ_L5) /* ~0.254828048790854 */

/* Relative Ionospheric Delay Scale: (f_L1 / f)^2 */
#define IONO_SCALE_L1 (1.0)
#define IONO_SCALE_L5 ((CARR_FREQ_L1 * CARR_FREQ_L1) / (CARR_FREQ_L5 * CARR_FREQ_L5)) /* ~1.79327 */

typedef struct {
	gnss_band_t band;
	const char *name;           /* "L1" or "L5" */
	double carr_freq;           /* Carrier frequency [Hz] */
	double lambda;              /* Carrier wavelength [m] */
	double code_freq;           /* Code chipping rate [Hz] */
	double carr_to_code;        /* Ratio: code_freq / carr_freq */
	int code_seq_len;           /* Code chips per 1 ms (1023 or 10230) */
	double iono_scale;          /* Relative ionospheric scale factor */
	double default_samp_freq;   /* Default sampling frequency [Hz] */
} band_param_t;

static inline band_param_t get_band_param(gnss_band_t band)
{
	band_param_t p;
	if (band == BAND_L5) {
		p.band = BAND_L5;
		p.name = "L5";
		p.carr_freq = CARR_FREQ_L5;
		p.lambda = LAMBDA_L5;
		p.code_freq = CODE_FREQ_L5;
		p.carr_to_code = CARR_TO_CODE_L5;
		p.code_seq_len = CA_SEQ_LEN_L5;
		p.iono_scale = IONO_SCALE_L5;
		p.default_samp_freq = 25.0e6;
	} else {
		p.band = BAND_L1;
		p.name = "L1";
		p.carr_freq = CARR_FREQ_L1;
		p.lambda = LAMBDA_L1;
		p.code_freq = CODE_FREQ_L1;
		p.carr_to_code = CARR_TO_CODE_L1;
		p.code_seq_len = CA_SEQ_LEN_L1;
		p.iono_scale = IONO_SCALE_L1;
		p.default_samp_freq = 2.6e6;
	}
	return p;
}

// Sampling data format
#define SC01 (1)
#define SC08 (8)
#define SC16 (16)

#define EPHEM_ARRAY_SIZE (15) // for daily GPS broadcast ephemers file (brdc)

/*! \brief Structure representing GPS time */
typedef struct
{
	int week;	/*!< GPS week number (since January 1980) */
	double sec; 	/*!< second inside the GPS \a week */
} gpstime_t;

/*! \brief Structure repreenting UTC time */
typedef struct
{
	int y; 		/*!< Calendar year */
	int m;		/*!< Calendar month */
	int d;		/*!< Calendar day */
	int hh;		/*!< Calendar hour */
	int mm;		/*!< Calendar minutes */
	double sec;	/*!< Calendar seconds */
} datetime_t;

/*! \brief Structure representing ephemeris of a single satellite */
typedef struct
{
	int vflg;	/*!< Valid Flag */
	datetime_t t;
	gpstime_t toc;	/*!< Time of Clock */
	gpstime_t toe;	/*!< Time of Ephemeris */
	int iodc;	/*!< Issue of Data, Clock */
	int iode;	/*!< Isuse of Data, Ephemeris */
	double deltan;	/*!< Delta-N (radians/sec) */
	double cuc;	/*!< Cuc (radians) */
	double cus;	/*!< Cus (radians) */
	double cic;	/*!< Correction to inclination cos (radians) */
	double cis;	/*!< Correction to inclination sin (radians) */
	double crc;	/*!< Correction to radius cos (meters) */
	double crs;	/*!< Correction to radius sin (meters) */
	double ecc;	/*!< e Eccentricity */
	double sqrta;	/*!< sqrt(A) (sqrt(m)) */
	double m0;	/*!< Mean anamoly (radians) */
	double omg0;	/*!< Longitude of the ascending node (radians) */
	double inc0;	/*!< Inclination (radians) */
	double aop;
	double omgdot;	/*!< Omega dot (radians/s) */
	double idot;	/*!< IDOT (radians/s) */
	double af0;	/*!< Clock offset (seconds) */
	double af1;	/*!< rate (sec/sec) */
	double af2;	/*!< acceleration (sec/sec^2) */
	double tgd;	/*!< Group delay L2 bias */
	int svhlth;
	int codeL2;
	// Working variables follow
	double n; 	/*!< Mean motion (Average angular velocity) */
	double sq1e2;	/*!< sqrt(1-e^2) */
	double A;	/*!< Semi-major axis */
	double omgkdot; /*!< OmegaDot-OmegaEdot */
} ephem_t;

typedef struct
{
	int enable;
	int vflg;
	double alpha0,alpha1,alpha2,alpha3;
	double beta0,beta1,beta2,beta3;
	double A0,A1;
	int dtls,tot,wnt;
	int dtlsf,dn,wnlsf;
	// enable custom leap event
	int leapen;
} ionoutc_t;

typedef struct
{
	gpstime_t g;
	double range; /*!< Pseudorange for target band [m] */
	double rate;  /*!< Pseudorange rate [m/s] */
	double d;     /*!< Geometric distance [m] */
	double azel[2]; /*!< Azimuth and elevation [rad] */
	double iono_delay; /*!< Ionospheric delay for target band [m] */
	double iono_delay_l1; /*!< Base L1 ionospheric delay [m] */
	double clk_bias; /*!< SV clock bias c * clk [m] */
} range_t;

#include "l5_code.h"
#include "cnav.h"

/*! \brief Structure representing a Channel */
typedef struct
{
	int prn;	/*< PRN Number */
	int ca[CA_SEQ_LEN]; /*< C/A Sequence */
	double f_carr;	/*< Carrier frequency */
	double f_code;	/*< Code frequency */
#ifdef FLOAT_CARR_PHASE
	double carr_phase;
#else
	unsigned int carr_phase; /*< Carrier phase */
	int carr_phasestep;	/*< Carrier phasestep */
#endif
	double code_phase; /*< Code phase */
	gpstime_t g0;	/*!< GPS time at start */
	unsigned long sbf[5][N_DWRD_SBF]; /*!< current subframe */
	unsigned long dwrd[N_DWRD]; /*!< Data words of sub-frame */
	int iword;	/*!< initial word */
	int ibit;	/*!< initial bit */
	int icode;	/*!< initial code */
	int dataBit;	/*!< current data bit */
	int codeCA;	/*!< current C/A code */
	double azel[2];
	range_t rho0;

	/* L5 signal parameters */
	int l5i[CA_SEQ_LEN_L5];      /*!< L5 I5 primary PRN code (+1 or -1) */
	int l5q[CA_SEQ_LEN_L5];      /*!< L5 Q5 primary PRN code (+1 or -1) */
	cnav_channel_t cnav;         /*!< CNAV channel state */
	int cnav_symbol;             /*!< Current CNAV symbol (+1 or -1) */
	int nh10_chip;               /*!< Current NH10 chip (+1 or -1) */
	int nh20_chip;               /*!< Current NH20 chip (+1 or -1) */
	int codeL5I;                 /*!< Current L5 I5 code chip (+1 or -1) */
	int codeL5Q;                 /*!< Current L5 Q5 code chip (+1 or -1) */
	int ms_epoch;                /*!< 1 ms epoch index (0..19) */
	int ms_in_msg;               /*!< Millisecond within 6s CNAV message (0..5999) */
} channel_t;

/*! \brief Structure representing a synchronized simulation stream (e.g. L1 or L5) */
typedef struct
{
	gnss_band_t band;
	band_param_t param;
	char outfile[MAX_CHAR];
	FILE *fp;
	double samp_freq;
	int iq_buff_size;       /*!< Complex samples per 0.1sec block */
	double delt;            /*!< 1.0 / samp_freq */
	short *iq_buff;         /*!< 16-bit complex I/Q buffer */
	signed char *iq8_buff;  /*!< 8-bit or 1-bit packed buffer */
	channel_t chan[MAX_CHAN];
	int gain[MAX_CHAN];
	int use_cuda;
} sim_stream_t;

#endif
