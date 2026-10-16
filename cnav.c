#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gpssim.h"
#include "cnav.h"

#define CRC24Q_POLY 0x864CFB

/* Bit-packing utilities */
static void pack_uint(int *bits, int offset, int nbits, unsigned long long val)
{
	for (int i = 0; i < nbits; i++) {
		bits[offset + i] = (int)((val >> (nbits - 1 - i)) & 1ULL);
	}
}

static void pack_sint(int *bits, int offset, int nbits, long long val)
{
	unsigned long long uval = (unsigned long long)val;
	unsigned long long mask = (nbits == 64) ? ~0ULL : ((1ULL << nbits) - 1ULL);
	uval &= mask;
	pack_uint(bits, offset, nbits, uval);
}

/* Calculate CRC-24Q parity */
unsigned int cnav_crc24q(const int *bits, int nbits)
{
	unsigned int crc = 0;
	for (int i = 0; i < nbits; i++) {
		int bit = bits[i] & 1;
		int msb = ((crc >> 23) & 1) ^ bit;
		crc = ((crc << 1) & 0xFFFFFF);
		if (msb) {
			crc ^= CRC24Q_POLY;
		}
	}
	return crc & 0xFFFFFF;
}

/* Common CNAV Header builder (Bits 0..37) */
static void build_cnav_header(int *bits, int prn, int type_id, unsigned int tow_count)
{
	/* Preamble (8 bits, 0x8B = 10001011) */
	pack_uint(bits, 0, 8, CNAV_PREAMBLE);

	/* PRN ID (6 bits): 1..32 for GPS, 1..10 for QZSS (PRN 193..202 -> 1..10) */
	int prn_id = (prn <= 32) ? prn : (prn - 192);
	pack_uint(bits, 8, 6, (unsigned int)(prn_id & 0x3F));

	/* Message Type ID (6 bits) */
	pack_uint(bits, 14, 6, (unsigned int)(type_id & 0x3F));

	/* Message TOW count (17 bits) */
	pack_uint(bits, 20, 17, tow_count & 0x1FFFF);

	/* Alert Flag (1 bit, 0 = normal) */
	bits[37] = 0;
}

/* Append 24-bit CRC-24Q to bits 276..299 */
static void append_cnav_crc(int *bits)
{
	unsigned int crc = cnav_crc24q(bits, 276);
	pack_uint(bits, 276, 24, crc);
}

/* Build CNAV Message Type 10 (Ephemeris 1) */
void cnav_build_type10(const ephem_t *eph, int prn, unsigned int tow_count, int *msg_bits)
{
	memset(msg_bits, 0, sizeof(int) * CNAV_MSG_BITS);

	/* 1. Header (Bits 0..37) */
	build_cnav_header(msg_bits, prn, CNAV_TYPE_10, tow_count);

	/* 2. Payload (Bits 38..275, 238 bits) */
	/* Week Number (13 bits, bits 38..50) */
	pack_uint(msg_bits, 38, 13, (unsigned long long)(eph->toe.week & 0x1FFF));

	/* Signal Health L1, L2, L5 (3 bits, bits 51..53) */
	msg_bits[51] = (eph->svhlth != 0) ? 1 : 0;
	msg_bits[52] = (eph->svhlth != 0) ? 1 : 0;
	msg_bits[53] = (eph->svhlth != 0) ? 1 : 0;

	/* Top: Data predict time of week (11 bits, scale factor 300s, bits 54..64) */
	unsigned int top = (unsigned int)floor(eph->toe.sec / 300.0) & 0x7FF;
	pack_uint(msg_bits, 54, 11, top);

	/* URAoe Index (5 bits signed, bits 65..69): nominal 2 */
	pack_sint(msg_bits, 65, 5, 2);

	/* Delta A: Semi-major axis difference relative to 26,559,710m (26 bits signed, scale 2^-9, bits 70..95) */
	double delta_A_m = eph->A - 26559710.0;
	long long delta_A = (long long)round(delta_A_m / 0.001953125); /* / 2^-9 */
	pack_sint(msg_bits, 70, 26, delta_A);

	/* A_DOT: Change rate in semi-major axis (25 bits signed, scale 2^-21, bits 96..120): 0 */
	pack_sint(msg_bits, 96, 25, 0);

	/* Delta n0: Mean motion difference (17 bits signed, scale 2^-44 semi-circles/s, bits 121..137) */
	long long delta_n0 = (long long)round((eph->deltan / PI) / pow(2.0, -44));
	pack_sint(msg_bits, 121, 17, delta_n0);

	/* Delta n0_DOT: Rate of mean motion difference (23 bits signed, scale 2^-57, bits 138..160): 0 */
	pack_sint(msg_bits, 138, 23, 0);

	/* M0: Mean anomaly at reference time (33 bits signed, scale 2^-32 semi-circles, bits 161..193) */
	long long m0 = (long long)round((eph->m0 / PI) / pow(2.0, -32));
	pack_sint(msg_bits, 161, 33, m0);

	/* e: Eccentricity (33 bits unsigned, scale 2^-34, bits 194..226) */
	unsigned long long ecc = (unsigned long long)round(eph->ecc / pow(2.0, -34));
	pack_uint(msg_bits, 194, 33, ecc);

	/* omega: Argument of perigee (33 bits signed, scale 2^-32 semi-circles, bits 227..259) */
	long long aop = (long long)round((eph->aop / PI) / pow(2.0, -32));
	pack_sint(msg_bits, 227, 33, aop);

	/* Toe: Ephemeris data reference time (11 bits unsigned, scale 300s, bits 260..270) */
	unsigned int toe = (unsigned int)round(eph->toe.sec / 300.0) & 0x7FF;
	pack_uint(msg_bits, 260, 11, toe);

	/* Reserved: 5 bits (bits 271..275) */
	pack_uint(msg_bits, 271, 5, 0);

	/* 3. Parity: CRC-24Q (bits 276..299) */
	append_cnav_crc(msg_bits);
}

/* Build CNAV Message Type 11 (Ephemeris 2) */
void cnav_build_type11(const ephem_t *eph, int prn, unsigned int tow_count, int *msg_bits)
{
	memset(msg_bits, 0, sizeof(int) * CNAV_MSG_BITS);

	/* 1. Header (Bits 0..37) */
	build_cnav_header(msg_bits, prn, CNAV_TYPE_11, tow_count);

	/* 2. Payload (Bits 38..275, 238 bits) */
	/* Toe: Ephemeris reference time (11 bits unsigned, scale 300s, bits 38..48) */
	unsigned int toe = (unsigned int)round(eph->toe.sec / 300.0) & 0x7FF;
	pack_uint(msg_bits, 38, 11, toe);

	/* Omega0: Reference right ascension angle (33 bits signed, scale 2^-32 semi-circles, bits 49..81) */
	long long omg0 = (long long)round((eph->omg0 / PI) / pow(2.0, -32));
	pack_sint(msg_bits, 49, 33, omg0);

	/* Delta Omega_DOT: Difference from reference -2.6e-9 (17 bits signed, scale 2^-44, bits 82..98) */
	long long delta_omgdot = (long long)round(((eph->omgdot / PI) - (-2.6e-9)) / pow(2.0, -44));
	pack_sint(msg_bits, 82, 17, delta_omgdot);

	/* i0: Inclination angle at reference time (33 bits signed, scale 2^-32 semi-circles, bits 99..131) */
	long long inc0 = (long long)round((eph->inc0 / PI) / pow(2.0, -32));
	pack_sint(msg_bits, 99, 33, inc0);

	/* i_DOT: Rate of inclination angle (15 bits signed, scale 2^-44 semi-circles/s, bits 132..146) */
	long long idot = (long long)round((eph->idot / PI) / pow(2.0, -44));
	pack_sint(msg_bits, 132, 15, idot);

	/* Cis: Amplitude of sine harmonic correction to inclination (16 bits signed, scale 2^-30 rad, bits 147..162) */
	long long cis = (long long)round(eph->cis / pow(2.0, -30));
	pack_sint(msg_bits, 147, 16, cis);

	/* Cic: Amplitude of cosine harmonic correction to inclination (16 bits signed, scale 2^-30 rad, bits 163..178) */
	long long cic = (long long)round(eph->cic / pow(2.0, -30));
	pack_sint(msg_bits, 163, 16, cic);

	/* Crs: Amplitude of sine correction to orbit radius (24 bits signed, scale 2^-8 m, bits 179..202) */
	long long crs = (long long)round(eph->crs / pow(2.0, -8));
	pack_sint(msg_bits, 179, 24, crs);

	/* Crc: Amplitude of cosine correction to orbit radius (24 bits signed, scale 2^-8 m, bits 203..226) */
	long long crc = (long long)round(eph->crc / pow(2.0, -8));
	pack_sint(msg_bits, 203, 24, crc);

	/* Cus: Amplitude of sine correction to argument of latitude (21 bits signed, scale 2^-30 rad, bits 227..247) */
	long long cus = (long long)round(eph->cus / pow(2.0, -30));
	pack_sint(msg_bits, 227, 24 - 3, cus);

	/* Cuc: Amplitude of cosine correction to argument of latitude (21 bits signed, scale 2^-30 rad, bits 248..268) */
	long long cuc = (long long)round(eph->cuc / pow(2.0, -30));
	pack_sint(msg_bits, 248, 21, cuc);

	/* Reserved: 7 bits (bits 269..275) */
	pack_uint(msg_bits, 269, 7, 0);

	/* 3. Parity: CRC-24Q (bits 276..299) */
	append_cnav_crc(msg_bits);
}

/* Build CNAV Message Type 30 (Clock, Ionosphere & Group Delay) */
void cnav_build_type30(const ephem_t *eph, const ionoutc_t *ionoutc, int prn, unsigned int tow_count, int *msg_bits)
{
	memset(msg_bits, 0, sizeof(int) * CNAV_MSG_BITS);

	/* 1. Header (Bits 0..37) */
	build_cnav_header(msg_bits, prn, CNAV_TYPE_30, tow_count);

	/* 2. Payload (Bits 38..275, 238 bits) */
	/* Top: Data predict time of week (11 bits unsigned, scale 300s, bits 38..48) */
	unsigned int top = (unsigned int)floor(eph->toc.sec / 300.0) & 0x7FF;
	pack_uint(msg_bits, 38, 11, top);

	/* Toc: Clock reference time of week (11 bits unsigned, scale 300s, bits 49..59) */
	unsigned int toc = (unsigned int)floor(eph->toc.sec / 300.0) & 0x7FF;
	pack_uint(msg_bits, 49, 11, toc);

	/* URAoc Index (5 bits signed, bits 60..64): 2 */
	pack_sint(msg_bits, 60, 5, 2);

	/* URAoc1 Index (3 bits unsigned, bits 65..67): 0 */
	pack_uint(msg_bits, 65, 3, 0);

	/* URAoc2 Index (3 bits unsigned, bits 68..70): 0 */
	pack_uint(msg_bits, 68, 3, 0);

	/* af0: Clock bias (26 bits signed, scale 2^-35 s, bits 71..96) */
	long long af0 = (long long)round(eph->af0 / pow(2.0, -35));
	pack_sint(msg_bits, 71, 26, af0);

	/* af1: Clock drift (20 bits signed, scale 2^-48 sec/sec, bits 97..116) */
	long long af1 = (long long)round(eph->af1 / pow(2.0, -48));
	pack_sint(msg_bits, 97, 20, af1);

	/* af2: Clock drift rate (10 bits signed, scale 2^-60 sec/sec^2, bits 117..126) */
	long long af2 = (long long)round(eph->af2 / pow(2.0, -60));
	pack_sint(msg_bits, 117, 10, af2);

	/* TGD: Group delay differential (13 bits signed, scale 2^-35 s, bits 127..139) */
	long long tgd = (long long)round(eph->tgd / pow(2.0, -35));
	pack_sint(msg_bits, 127, 13, tgd);

	/* ISC_L1CA, ISC_L2C, ISC_L5I5, ISC_L5Q5 (4 x 13 bits signed, bits 140..191): nominal 0 */
	pack_sint(msg_bits, 140, 13, 0);
	pack_sint(msg_bits, 153, 13, 0);
	pack_sint(msg_bits, 166, 13, 0);
	pack_sint(msg_bits, 179, 13, 0);

	/* Klobuchar Ionosphere Parameters (alpha0..3, beta0..3: 8 x 8 bits signed, bits 192..255) */
	long a0 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->alpha0 / pow(2.0, -30)) : 0);
	long a1 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->alpha1 / pow(2.0, -27)) : 0);
	long a2 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->alpha2 / pow(2.0, -24)) : 0);
	long a3 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->alpha3 / pow(2.0, -24)) : 0);
	long b0 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->beta0 / 2048.0) : 0);
	long b1 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->beta1 / 16384.0) : 0);
	long b2 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->beta2 / 65536.0) : 0);
	long b3 = (long)round((ionoutc && ionoutc->enable) ? (ionoutc->beta3 / 65536.0) : 0);

	pack_sint(msg_bits, 192, 8, a0);
	pack_sint(msg_bits, 200, 8, a1);
	pack_sint(msg_bits, 208, 8, a2);
	pack_sint(msg_bits, 216, 8, a3);
	pack_sint(msg_bits, 224, 8, b0);
	pack_sint(msg_bits, 232, 8, b1);
	pack_sint(msg_bits, 240, 8, b2);
	pack_sint(msg_bits, 248, 8, b3);

	/* Reserved: 20 bits (bits 256..275) */
	pack_uint(msg_bits, 256, 20, 0);

	/* 3. Parity: CRC-24Q (bits 276..299) */
	append_cnav_crc(msg_bits);
}

/* Convolutional FEC Encoder */
void cnav_fec_init(cnav_fec_t *fec)
{
	if (fec)
		fec->state = 0;
}

void cnav_fec_encode_bit(cnav_fec_t *fec, int bit, int *sym1, int *sym2)
{
	/* Register has 7 bits: [bit (MSB, bit 6), state (6 bits, bits 5..0)] */
	unsigned int reg = (((unsigned int)(bit & 1)) << 6) | (fec->state & 0x3F);

	/* G1 = 171 octal = 0b1111001 (taps: 6, 5, 4, 3, 0) */
	/* G2 = 133 octal = 0b1011011 (taps: 6, 4, 3, 1, 0) */
	unsigned int g1_val = reg & 0x79;
	unsigned int g2_val = reg & 0x5B;

	g1_val ^= g1_val >> 1;
	g1_val ^= g1_val >> 2;
	g1_val ^= g1_val >> 4;
	*sym1 = (int)(g1_val & 1);

	g2_val ^= g2_val >> 1;
	g2_val ^= g2_val >> 2;
	g2_val ^= g2_val >> 4;
	*sym2 = (int)(g2_val & 1);

	/* Advance shift register: shift right */
	fec->state = reg >> 1;
}

void cnav_fec_encode_msg(cnav_fec_t *fec, const int *bits, int nbits, int *symbols)
{
	cnav_fec_t local_fec;
	if (!fec) {
		cnav_fec_init(&local_fec);
		fec = &local_fec;
	}

	for (int i = 0; i < nbits; i++) {
		int s1, s2;
		cnav_fec_encode_bit(fec, bits[i], &s1, &s2);
		symbols[i * 2]     = s1;
		symbols[i * 2 + 1] = s2;
	}
}

/* Channel streaming */
void cnav_channel_build_msg(cnav_channel_t *chan)
{
	if (chan->cycle_idx == 0) {
		cnav_build_type10(&chan->eph, chan->prn, chan->tow_count, chan->msg_bits);
	} else if (chan->cycle_idx == 1) {
		cnav_build_type11(&chan->eph, chan->prn, chan->tow_count, chan->msg_bits);
	} else {
		cnav_build_type30(&chan->eph, &chan->ionoutc, chan->prn, chan->tow_count, chan->msg_bits);
	}

	/* Encode 300 data bits into 600 symbols */
	cnav_fec_encode_msg(&chan->fec, chan->msg_bits, CNAV_MSG_BITS, chan->msg_symbols);
}

void cnav_channel_init(cnav_channel_t *chan, const ephem_t *eph, const ionoutc_t *ionoutc, int prn, double start_tow)
{
	memset(chan, 0, sizeof(cnav_channel_t));
	chan->prn = prn;
	if (eph) chan->eph = *eph;
	if (ionoutc) chan->ionoutc = *ionoutc;
	chan->start_tow = start_tow;

	cnav_fec_init(&chan->fec);

	/* Calculate initial 6-second message TOW count */
	chan->tow_count = (unsigned int)floor(start_tow / 6.0);
	chan->cycle_idx = chan->tow_count % CNAV_CYCLE_MSGS;

	cnav_channel_build_msg(chan);

	/* Determine offset within the 6-second message (0..599 symbols, 0..9 ms) */
	double sec_in_msg = start_tow - (floor(start_tow / 6.0) * 6.0);
	int total_ms = (int)floor(sec_in_msg * 1000.0);
	if (total_ms < 0) total_ms = 0;
	if (total_ms >= 6000) total_ms = 5999;

	chan->symbol_idx = total_ms / 10;
	chan->ms_in_symbol = total_ms % 10;
}

void cnav_channel_step_ms(cnav_channel_t *chan)
{
	chan->ms_in_symbol++;
	if (chan->ms_in_symbol >= 10) {
		chan->ms_in_symbol = 0;
		chan->symbol_idx++;

		if (chan->symbol_idx >= CNAV_MSG_SYMBOLS) {
			chan->symbol_idx = 0;
			chan->tow_count++;
			chan->cycle_idx = (chan->cycle_idx + 1) % CNAV_CYCLE_MSGS;
			cnav_channel_build_msg(chan);
		}
	}
}

int cnav_channel_get_symbol(const cnav_channel_t *chan)
{
	/* Map binary symbol 0 -> +1, 1 -> -1 */
	return (chan->msg_symbols[chan->symbol_idx] == 0) ? 1 : -1;
}

void cnav_channel_get_2msg_symbols(const cnav_channel_t *chan, signed char *symbols_1200)
{
	/* 1. First 600 symbols: current message */
	for (int i = 0; i < CNAV_MSG_SYMBOLS; i++) {
		symbols_1200[i] = (chan->msg_symbols[i] == 0) ? 1 : -1;
	}

	/* 2. Next 600 symbols: next consecutive message (tow_count + 1) */
	int next_bits[CNAV_MSG_BITS];
	int next_symbols[CNAV_MSG_SYMBOLS];
	unsigned int next_tow = chan->tow_count + 1;
	int next_cycle = next_tow % CNAV_CYCLE_MSGS;

	if (next_cycle == 0) {
		cnav_build_type10(&chan->eph, chan->prn, next_tow, next_bits);
	} else if (next_cycle == 1) {
		cnav_build_type11(&chan->eph, chan->prn, next_tow, next_bits);
	} else {
		cnav_build_type30(&chan->eph, &chan->ionoutc, chan->prn, next_tow, next_bits);
	}

	cnav_fec_t next_fec = chan->fec;
	cnav_fec_encode_msg(&next_fec, next_bits, CNAV_MSG_BITS, next_symbols);

	for (int i = 0; i < CNAV_MSG_SYMBOLS; i++) {
		symbols_1200[CNAV_MSG_SYMBOLS + i] = (next_symbols[i] == 0) ? 1 : -1;
	}
}

/* Self-test and verification */
int test_cnav_verification(int verbose)
{
	int errors = 0;

	if (verbose) {
		printf("======================================================================\n");
		printf(" CNAV Message Generation, CRC-24Q & Rate-1/2 FEC Verification Suite\n");
		printf("======================================================================\n\n");
	}

	/* 1. Test CRC-24Q calculation and self-checking remainder property */
	int test_msg[300];
	memset(test_msg, 0, sizeof(test_msg));
	test_msg[0] = 1; test_msg[4] = 1; test_msg[6] = 1; test_msg[7] = 1; /* Preamble 10001011 */
	for (int i = 8; i < 276; i++) {
		test_msg[i] = (i * 37 + 13) % 2;
	}
	unsigned int crc = cnav_crc24q(test_msg, 276);
	pack_uint(test_msg, 276, 24, crc);

	/* Full 300-bit CRC check */
	unsigned int check_crc = cnav_crc24q(test_msg, 300);
	if (check_crc != 0) {
		if (verbose) printf("[FAIL] CRC-24Q check failed! Remainder = 0x%06X (expected 0x0)\n", check_crc);
		errors++;
	} else if (verbose) {
		printf("[PASS] CRC-24Q parity verified: remainder over full 300 bits is 0x000000 (CRC=0x%06X)\n", crc);
	}

	/* 2. Test Message Type 10, 11, 30 packing with sample ephemeris */
	ephem_t eph;
	memset(&eph, 0, sizeof(ephem_t));
	eph.vflg = 1;
	eph.toe.week = 2190;
	eph.toe.sec = 518400.0;
	eph.toc.sec = 518400.0;
	eph.A = 26559710.0 + 1000.0;
	eph.deltan = 4.5e-9;
	eph.m0 = 1.0;
	eph.ecc = 0.005;
	eph.aop = 2.0;
	eph.omg0 = 0.5;
	eph.omgdot = -4.5e-9;
	eph.inc0 = 0.96;
	eph.idot = 1.2e-10;
	eph.cis = 1.0e-7;
	eph.cic = 2.0e-7;
	eph.crs = 10.0;
	eph.crc = 200.0;
	eph.cus = 3.0e-6;
	eph.cuc = 4.0e-6;
	eph.af0 = 1.2e-4;
	eph.af1 = -1.5e-12;
	eph.af2 = 0.0;
	eph.tgd = 2.5e-9;

	ionoutc_t iono;
	memset(&iono, 0, sizeof(ionoutc_t));
	iono.enable = 1;
	iono.alpha0 = 0.1118e-07;
	iono.alpha1 = -0.7451e-08;
	iono.beta0 = 0.1004e+06;

	int bits_mt10[300], bits_mt11[300], bits_mt30[300];
	cnav_build_type10(&eph, 1, 86400, bits_mt10);
	cnav_build_type11(&eph, 1, 86401, bits_mt11);
	cnav_build_type30(&eph, &iono, 1, 86402, bits_mt30);

	/* Verify preamble */
	unsigned int pre10 = 0, pre11 = 0, pre30 = 0;
	for (int i = 0; i < 8; i++) {
		pre10 = (pre10 << 1) | bits_mt10[i];
		pre11 = (pre11 << 1) | bits_mt11[i];
		pre30 = (pre30 << 1) | bits_mt30[i];
	}
	if (pre10 != CNAV_PREAMBLE || pre11 != CNAV_PREAMBLE || pre30 != CNAV_PREAMBLE) {
		if (verbose) printf("[FAIL] Preamble mismatch: 0x%02X, 0x%02X, 0x%02X\n", pre10, pre11, pre30);
		errors++;
	} else if (verbose) {
		printf("[PASS] Preambles verified (0x8B = 10001011)\n");
	}

	/* Verify CRC on all 3 message types */
	if (cnav_crc24q(bits_mt10, 300) != 0 ||
	    cnav_crc24q(bits_mt11, 300) != 0 ||
	    cnav_crc24q(bits_mt30, 300) != 0) {
		if (verbose) printf("[FAIL] Generated CNAV message CRC check failed!\n");
		errors++;
	} else if (verbose) {
		printf("[PASS] Message Type 10, 11, 30 generation and CRC-24Q parity verified!\n");
	}

	/* 3. Test QZSS PRN wrapping (PRN 193 -> 1) */
	int qzss_bits[300];
	cnav_build_type10(&eph, 193, 86400, qzss_bits);
	unsigned int qzss_prn_id = 0;
	for (int i = 8; i < 14; i++) qzss_prn_id = (qzss_prn_id << 1) | qzss_bits[i];
	if (qzss_prn_id != 1 || cnav_crc24q(qzss_bits, 300) != 0) {
		if (verbose) printf("[FAIL] QZSS PRN 193 CNAV packing failed (prn_id=%d)\n", qzss_prn_id);
		errors++;
	} else if (verbose) {
		printf("[PASS] QZSS PRN 193 -> PRN ID %d CNAV packing & CRC verified!\n", qzss_prn_id);
	}

	/* 4. Test Convolutional Encoder Rate 1/2 (K=7, G1=171 oct, G2=133 oct) */
	cnav_fec_t fec;
	cnav_fec_init(&fec);
	int s1, s2;
	cnav_fec_encode_bit(&fec, 1, &s1, &s2);
	/* Initial state 0. Bit 1 in MSB of 7-bit reg: 0b1000000 = 0x40.
	 * G1 = 0b1111001: 0x40 & 0x79 = 0x40 -> parity 1.
	 * G2 = 0b1011011: 0x40 & 0x5B = 0x40 -> parity 1.
	 */
	if (s1 != 1 || s2 != 1) {
		if (verbose) printf("[FAIL] FEC initial bit encode mismatch: s1=%d, s2=%d (expected 1, 1)\n", s1, s2);
		errors++;
	} else if (verbose) {
		printf("[PASS] Rate 1/2 Convolutional FEC encoder impulse response verified.\n");
	}

	/* 5. Test CNAV channel streaming & 10ms symbol alignment */
	cnav_channel_t chan;
	cnav_channel_init(&chan, &eph, &iono, 1, 518400.0);
	int sym_count = 0;
	for (int ms = 0; ms < 6000; ms++) {
		int sym = cnav_channel_get_symbol(&chan);
		if (sym != 1 && sym != -1) {
			errors++;
			break;
		}
		if (chan.ms_in_symbol == 0) sym_count++;
		cnav_channel_step_ms(&chan);
	}
	if (sym_count != 600) {
		if (verbose) printf("[FAIL] Expected 600 symbols in 6000 ms, got %d\n", sym_count);
		errors++;
	} else if (verbose) {
		printf("[PASS] CNAV 100 sps symbol stream timing (600 symbols / 6 sec) verified!\n");
	}

	if (verbose) {
		printf("\n======================================================================\n");
		if (errors == 0)
			printf(" All CNAV verification tests PASSED successfully!\n");
		else
			printf(" CNAV verification tests FAILED with %d errors!\n", errors);
		printf("======================================================================\n\n");
	}

	return errors;
}
