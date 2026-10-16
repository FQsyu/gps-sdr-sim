#ifndef CNAV_H
#define CNAV_H

/* Note: Requires types ephem_t and ionoutc_t from gpssim.h */
#ifndef GPSSIM_H
#include "gpssim.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Total number of bits in one CNAV message (6 seconds at 50 bps) */
#define CNAV_MSG_BITS 300

/*! \brief Number of payload data bits per CNAV message */
#define CNAV_PAYLOAD_BITS 238

/*! \brief Number of CRC-24Q parity check bits */
#define CNAV_CRC_BITS 24

/*! \brief Header bits: Preamble(8) + PRN(6) + Type(6) + TOW(17) + Alert(1) */
#define CNAV_HEADER_BITS 38

/*! \brief CNAV message preamble: 10001011 binary (0x8B) */
#define CNAV_PREAMBLE 0x8B

/*! \brief Total number of FEC symbols per 6-second message (100 sps at Rate 1/2) */
#define CNAV_MSG_SYMBOLS (CNAV_MSG_BITS * 2)

/* Supported CNAV Message Types per IS-GPS-705 and IS-QZSS-PNT */
#define CNAV_TYPE_10 10 /* Ephemeris 1 */
#define CNAV_TYPE_11 11 /* Ephemeris 2 */
#define CNAV_TYPE_30 30 /* Clock, IONO, Group Delay */

/*! \brief Standard broadcast cycle: Type 10 -> Type 11 -> Type 30 (18 seconds total) */
#define CNAV_CYCLE_MSGS 3
#define CNAV_CYCLE_SEC (CNAV_CYCLE_MSGS * 6)

/*! \brief Calculate CRC-24Q parity over an array of bits (each 0 or 1).
 *  \param[in] bits Array of bits (each 0 or 1)
 *  \param[in] nbits Number of bits (typically 276)
 *  \return 24-bit CRC value
 */
unsigned int cnav_crc24q(const int *bits, int nbits);

/*! \brief Build CNAV Message Type 10 (Ephemeris 1).
 *  \param[in] eph Satellite ephemeris from RINEX
 *  \param[in] prn Satellite PRN number
 *  \param[in] tow_count TOW count (in 6-second units)
 *  \param[out] msg_bits Output array of 300 bits (including preamble and CRC-24Q)
 */
void cnav_build_type10(const ephem_t *eph, int prn, unsigned int tow_count, int *msg_bits);

/*! \brief Build CNAV Message Type 11 (Ephemeris 2).
 *  \param[in] eph Satellite ephemeris from RINEX
 *  \param[in] prn Satellite PRN number
 *  \param[in] tow_count TOW count (in 6-second units)
 *  \param[out] msg_bits Output array of 300 bits (including preamble and CRC-24Q)
 */
void cnav_build_type11(const ephem_t *eph, int prn, unsigned int tow_count, int *msg_bits);

/*! \brief Build CNAV Message Type 30 (Clock, Ionosphere & Group Delay).
 *  \param[in] eph Satellite ephemeris from RINEX
 *  \param[in] ionoutc Ionospheric and UTC parameters from RINEX
 *  \param[in] prn Satellite PRN number
 *  \param[in] tow_count TOW count (in 6-second units)
 *  \param[out] msg_bits Output array of 300 bits (including preamble and CRC-24Q)
 */
void cnav_build_type30(const ephem_t *eph, const ionoutc_t *ionoutc, int prn, unsigned int tow_count, int *msg_bits);

/*! \brief Rate 1/2 Convolutional FEC Encoder (K=7, G1=171 oct, G2=133 oct) */
typedef struct {
	unsigned int state; /* 6 memory shift register stages */
} cnav_fec_t;

/*! \brief Initialize or reset the convolutional FEC encoder to zero state */
void cnav_fec_init(cnav_fec_t *fec);

/*! \brief Encode a single data bit (50 bps) into 2 FEC symbols (100 sps).
 *  G1 is output first (0..10ms), followed by G2 (10..20ms).
 *  \param[in,out] fec Encoder state
 *  \param[in] bit Input data bit (0 or 1)
 *  \param[out] sym1 First symbol (G1, 0 or 1)
 *  \param[out] sym2 Second symbol (G2, 0 or 1)
 */
void cnav_fec_encode_bit(cnav_fec_t *fec, int bit, int *sym1, int *sym2);

/*! \brief Encode an array of bits into symbols using rate 1/2 convolutional code.
 *  \param[in,out] fec Encoder state (can be NULL for clean zero-initialized state)
 *  \param[in] bits Array of nbits (0 or 1)
 *  \param[in] nbits Number of input bits
 *  \param[out] symbols Array of 2 * nbits symbols (0 or 1)
 */
void cnav_fec_encode_msg(cnav_fec_t *fec, const int *bits, int nbits, int *symbols);

/*! \brief Structure representing CNAV generator state for a satellite channel */
typedef struct {
	int prn;
	ephem_t eph;
	ionoutc_t ionoutc;
	cnav_fec_t fec;

	double start_tow;          /* Simulation start TOW in seconds */
	unsigned int tow_count;    /* Current message TOW count (6-second epoch) */
	int cycle_idx;             /* 0: MT10, 1: MT11, 2: MT30 */

	int msg_bits[CNAV_MSG_BITS];       /* Current 300-bit message */
	int msg_symbols[CNAV_MSG_SYMBOLS]; /* Current 600 symbols (100 sps, 10 ms each) */

	int symbol_idx;            /* 0..599 (current 10 ms symbol index in 6s message) */
	int ms_in_symbol;          /* 0..9 (current 1 ms chip within 10 ms symbol) */
} cnav_channel_t;

/*! \brief Initialize a CNAV channel stream for a satellite.
 *  \param[out] chan CNAV channel
 *  \param[in] eph Satellite ephemeris
 *  \param[in] ionoutc Iono/UTC parameters
 *  \param[in] prn PRN number
 *  \param[in] start_tow GPS Time of Week in seconds
 */
void cnav_channel_init(cnav_channel_t *chan, const ephem_t *eph, const ionoutc_t *ionoutc, int prn, double start_tow);

/*! \brief Step the CNAV channel forward by 1 ms.
 *  Handles 10 ms symbol boundaries and 6-second message packet transitions.
 *  \param[in,out] chan CNAV channel
 */
void cnav_channel_step_ms(cnav_channel_t *chan);

/*! \brief Get current CNAV symbol (+1 or -1) for I5 channel.
 *  \param[in] chan CNAV channel
 *  \return Modulated sign (+1 or -1)
 */
int cnav_channel_get_symbol(const cnav_channel_t *chan);

/*! \brief Self-test and verification of CNAV message generation, CRC-24Q, and FEC encoding.
 *  \param[in] verbose If non-zero, print detailed verification logs to stdout.
 *  \return 0 on success, non-zero on failure.
 */
int test_cnav_verification(int verbose);

#ifdef __cplusplus
}
#endif

#endif /* CNAV_H */
