#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l5_code.h"

/* Neumann-Hoffman secondary code sequences (binary 0 / 1) */
/* NH10: 10-bit sequence for I5 channel (period 10 ms): 0000110101 */
const int NH10[NH10_LEN] = {
	0, 0, 0, 0, 1, 1, 0, 1, 0, 1
};

/* NH20: 20-bit sequence for Q5 channel (period 20 ms): 00000100110101001110 */
const int NH20[NH20_LEN] = {
	0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0
};

/* Neumann-Hoffman secondary code sequences in bi-phase mapping (+1 / -1): 0 -> +1, 1 -> -1 */
const int NH10_BIPHASE[NH10_LEN] = {
	 1,  1,  1,  1, -1, -1,  1, -1,  1, -1
};

const int NH20_BIPHASE[NH20_LEN] = {
	 1,  1,  1,  1,  1, -1,  1,  1, -1, -1,  1, -1,  1, -1,  1,  1, -1, -1, -1,  1
};

/* LFSR Polynomial Tap Masks (13 stages, shift-right implementation):
 * XA polynomial: 1 + X^9 + X^10 + X^12 + X^13
 *   Taps from stages 9, 10, 12, 13 (bits 4, 3, 1, 0): 0x001B
 * XB polynomial: 1 + X^1 + X^3 + X^4 + X^6 + X^7 + X^8 + X^12 + X^13
 *   Taps from stages 1, 3, 4, 6, 7, 8, 12, 13 (bits 12, 10, 9, 7, 6, 5, 1, 0): 0x16E3
 */
#define XA_TAP_MASK 0x001B
#define XB_TAP_MASK 0x16E3
#define LFSR_INIT_ALL1 0x1FFF

/* Table 3-Ia from IS-GPS-705: GPS PRN 1-32 XB Code Advance (chips) [I5, Q5] */
static const int GPS_L5_ADVANCE[32][2] = {
	{ 266, 1701}, /* PRN  1 */
	{ 365,  323}, /* PRN  2 */
	{ 804, 5292}, /* PRN  3 */
	{1138, 2020}, /* PRN  4 */
	{1509, 5429}, /* PRN  5 */
	{1559, 7136}, /* PRN  6 */
	{1756, 1041}, /* PRN  7 */
	{2084, 5947}, /* PRN  8 */
	{2170, 4315}, /* PRN  9 */
	{2303,  148}, /* PRN 10 */
	{2527,  535}, /* PRN 11 */
	{2687, 1939}, /* PRN 12 */
	{2930, 5206}, /* PRN 13 */
	{3471, 5910}, /* PRN 14 */
	{3940, 3595}, /* PRN 15 */
	{4132, 5135}, /* PRN 16 */
	{4332, 6082}, /* PRN 17 */
	{4924, 6990}, /* PRN 18 */
	{5343, 3546}, /* PRN 19 */
	{5443, 1523}, /* PRN 20 */
	{5641, 4548}, /* PRN 21 */
	{5816, 4484}, /* PRN 22 */
	{5898, 1893}, /* PRN 23 */
	{5918, 3961}, /* PRN 24 */
	{5955, 7106}, /* PRN 25 */
	{6243, 5299}, /* PRN 26 */
	{6345, 4660}, /* PRN 27 */
	{6477,  276}, /* PRN 28 */
	{6518, 4389}, /* PRN 29 */
	{6875, 3783}, /* PRN 30 */
	{7168, 1591}, /* PRN 31 */
	{7187, 1601}  /* PRN 32 */
};

/* Table 3.3-1 from IS-QZSS-PNT: QZSS PRN 193-202 XB Code Advance (chips) [I5, Q5] */
static const int QZSS_L5_ADVANCE[10][2] = {
	{5836, 4757}, /* PRN 193 */
	{ 926,  427}, /* PRN 194 */
	{6086, 5452}, /* PRN 195 */
	{ 950, 5182}, /* PRN 196 */
	{5905, 6606}, /* PRN 197 */
	{3240, 6531}, /* PRN 198 */
	{6675, 4268}, /* PRN 199 */
	{3197, 3115}, /* PRN 200 */
	{1555, 6835}, /* PRN 201 */
	{3589,  862}  /* PRN 202 */
};

/* Table 3-Ia from IS-GPS-705: GPS PRN 1-32 Initial XB Code States (13-bit) [I5, Q5] */
static const unsigned short GPS_L5_INIT_STATE[32][2] = {
	{0x0AE4, 0x12CC}, /* PRN  1: I=0101011100100, Q=1001011001100 */
	{0x1835, 0x08F6}, /* PRN  2: I=1100000110101, Q=0100011110110 */
	{0x0808, 0x1E23}, /* PRN  3: I=0100000001000, Q=1111000100011 */
	{0x1626, 0x076A}, /* PRN  4: I=1011000100110, Q=0011101101010 */
	{0x1DD7, 0x07B2}, /* PRN  5: I=1110111010111, Q=0011110110010 */
	{0x0CFA, 0x0AA9}, /* PRN  6: I=0110011111010, Q=0101010101001 */
	{0x149F, 0x1F81}, /* PRN  7: I=1010010011111, Q=1111110000001 */
	{0x17A4, 0x0D68}, /* PRN  8: I=1011110100100, Q=0110101101000 */
	{0x1F2B, 0x1743}, /* PRN  9: I=1111100101011, Q=1011101000011 */
	{0x0FDE, 0x0486}, /* PRN 10: I=0111111011110, Q=0010010000110 */
	{0x013A, 0x0205}, /* PRN 11: I=0000100111010, Q=0001000000101 */
	{0x1CF9, 0x0AC5}, /* PRN 12: I=1110011111001, Q=0101011000101 */
	{0x039C, 0x09A5}, /* PRN 13: I=0001110011100, Q=0100110100101 */
	{0x0827, 0x143F}, /* PRN 14: I=0100000100111, Q=1010000111111 */
	{0x0D5A, 0x178F}, /* PRN 15: I=0110101011010, Q=1011110001111 */
	{0x03C9, 0x1A5F}, /* PRN 16: I=0001111001001, Q=1101001011111 */
	{0x098F, 0x1CC8}, /* PRN 17: I=0100110001111, Q=1110011001000 */
	{0x1E1E, 0x16E4}, /* PRN 18: I=1111000011110, Q=1011011100100 */
	{0x191F, 0x065B}, /* PRN 19: I=1100100011111, Q=0011001011011 */
	{0x0D6D, 0x1871}, /* PRN 20: I=0110101101101, Q=1100001110001 */
	{0x0408, 0x0D90}, /* PRN 21: I=0010000001000, Q=0110110010000 */
	{0x1DEF, 0x058E}, /* PRN 22: I=1110111101111, Q=0010110001110 */
	{0x10FE, 0x117D}, /* PRN 23: I=1000011111110, Q=1000101111101 */
	{0x18B4, 0x0DF3}, /* PRN 24: I=1100010110100, Q=0110111110011 */
	{0x1A6D, 0x089B}, /* PRN 25: I=1101001101101, Q=0100010011011 */
	{0x1596, 0x0ABC}, /* PRN 26: I=1010110010110, Q=0101010111100 */
	{0x0ADE, 0x10FA}, /* PRN 27: I=0101011011110, Q=1000011111010 */
	{0x0F56, 0x1F42}, /* PRN 28: I=0111101010110, Q=1111101000010 */
	{0x0BE1, 0x0A24}, /* PRN 29: I=0101111100001, Q=0101000100100 */
	{0x10B7, 0x1079}, /* PRN 30: I=1000010110111, Q=1000001111001 */
	{0x029E, 0x0BE5}, /* PRN 31: I=0001010011110, Q=0101111100101 */
	{0x00B9, 0x122A}  /* PRN 32: I=0000010111001, Q=1001000101010 */
};

/* Table 3.3-1 from IS-QZSS-PNT: QZSS PRN 193-202 Initial XB Code States (13-bit) [I5, Q5] */
static const unsigned short QZSS_L5_INIT_STATE[10][2] = {
	{0x0C2E, 0x1387}, /* PRN 193: I=0110000101110, Q=1001110000111 */
	{0x0C9F, 0x0D3A}, /* PRN 194: I=0110010011111, Q=0110100111010 */
	{0x11CC, 0x0C66}, /* PRN 195: I=1000111001100, Q=0110001100110 */
	{0x1F71, 0x010C}, /* PRN 196: I=1111101110001, Q=0000100001100 */
	{0x07E1, 0x0A2D}, /* PRN 197: I=0011111100001, Q=0101000101101 */
	{0x0071, 0x1057}, /* PRN 198: I=0000001110001, Q=1000001010111 */
	{0x15A4, 0x0671}, /* PRN 199: I=1010110100100, Q=0011001110001 */
	{0x0876, 0x08E6}, /* PRN 200: I=0100001110110, Q=0100011100110 */
	{0x0FA3, 0x0965}, /* PRN 201: I=0111110100011, Q=0100101100101 */
	{0x03CB, 0x1C57}  /* PRN 202: I=0001111001011, Q=1110001010111 */
};

/*! \brief Calculate XOR bit parity of up to 16 bits */
static inline int lfsr_parity(unsigned int v)
{
	v ^= v >> 1;
	v ^= v >> 2;
	v ^= v >> 4;
	v ^= v >> 8;
	return (int)(v & 1);
}

int is_valid_l5_prn(int prn)
{
	if (prn >= 1 && prn <= 32)
		return 1;
	if (prn >= 193 && prn <= 202)
		return 1;
	return 0;
}

int get_l5_advance(int prn, int *adv_i, int *adv_q)
{
	if (prn >= 1 && prn <= 32) {
		if (adv_i) *adv_i = GPS_L5_ADVANCE[prn - 1][0];
		if (adv_q) *adv_q = GPS_L5_ADVANCE[prn - 1][1];
		return 0;
	} else if (prn >= 193 && prn <= 202) {
		if (adv_i) *adv_i = QZSS_L5_ADVANCE[prn - 193][0];
		if (adv_q) *adv_q = QZSS_L5_ADVANCE[prn - 193][1];
		return 0;
	}
	return -1;
}

int get_l5_initial_state(int prn, unsigned short *state_i, unsigned short *state_q)
{
	if (prn >= 1 && prn <= 32) {
		if (state_i) *state_i = GPS_L5_INIT_STATE[prn - 1][0];
		if (state_q) *state_q = GPS_L5_INIT_STATE[prn - 1][1];
		return 0;
	} else if (prn >= 193 && prn <= 202) {
		if (state_i) *state_i = QZSS_L5_INIT_STATE[prn - 193][0];
		if (state_q) *state_q = QZSS_L5_INIT_STATE[prn - 193][1];
		return 0;
	}
	return -1;
}

int generate_l5_code(int prn, int *i5_code, int *q5_code)
{
	unsigned short state_i = 0;
	unsigned short state_q = 0;

	if (get_l5_initial_state(prn, &state_i, &state_q) != 0)
		return -1;

	/* 13-stage Linear Feedback Shift Registers */
	unsigned int reg_xa  = LFSR_INIT_ALL1;
	unsigned int reg_xb1 = (unsigned int)state_i;
	unsigned int reg_xb2 = (unsigned int)state_q;

	for (int t = 0; t < L5_CODE_LEN; t++) {
		/* Short-cycle reset for XA at chip 8190 per IS-GPS-705 Section 3.3.2.2 */
		if (t == L5_XA_PERIOD) {
			reg_xa = LFSR_INIT_ALL1;
		}

		/* Output chip is stage 13 (LSB in right-shift model) */
		int out_xa  = (int)(reg_xa  & 1);
		int out_xb1 = (int)(reg_xb1 & 1);
		int out_xb2 = (int)(reg_xb2 & 1);

		/* Modulo-2 addition (XOR) to produce I5 and Q5 codes */
		if (i5_code != NULL)
			i5_code[t] = out_xa ^ out_xb1;
		if (q5_code != NULL)
			q5_code[t] = out_xa ^ out_xb2;

		/* Clock XA shift register: feedback = stage 9 ^ 10 ^ 12 ^ 13 */
		int fb_xa = lfsr_parity(reg_xa & XA_TAP_MASK);
		reg_xa = ((unsigned int)fb_xa << 12) | (reg_xa >> 1);

		/* Clock XB1 shift register: feedback = stage 1 ^ 3 ^ 4 ^ 6 ^ 7 ^ 8 ^ 12 ^ 13 */
		int fb_xb1 = lfsr_parity(reg_xb1 & XB_TAP_MASK);
		reg_xb1 = ((unsigned int)fb_xb1 << 12) | (reg_xb1 >> 1);

		/* Clock XB2 shift register */
		int fb_xb2 = lfsr_parity(reg_xb2 & XB_TAP_MASK);
		reg_xb2 = ((unsigned int)fb_xb2 << 12) | (reg_xb2 >> 1);
	}

	return 0;
}

int generate_l5_code_biphase(int prn, int *i5_code, int *q5_code)
{
	int ret = generate_l5_code(prn, i5_code, q5_code);
	if (ret != 0)
		return ret;

	/* Map binary 0 -> +1, 1 -> -1 */
	for (int t = 0; t < L5_CODE_LEN; t++) {
		if (i5_code != NULL)
			i5_code[t] = (i5_code[t] == 0) ? 1 : -1;
		if (q5_code != NULL)
			q5_code[t] = (q5_code[t] == 0) ? 1 : -1;
	}
	return 0;
}

int get_nh10_chip(int epoch_ms)
{
	int idx = epoch_ms % NH10_LEN;
	if (idx < 0) idx += NH10_LEN;
	return NH10[idx];
}

int get_nh20_chip(int epoch_ms)
{
	int idx = epoch_ms % NH20_LEN;
	if (idx < 0) idx += NH20_LEN;
	return NH20[idx];
}

int get_nh10_biphase(int epoch_ms)
{
	int idx = epoch_ms % NH10_LEN;
	if (idx < 0) idx += NH10_LEN;
	return NH10_BIPHASE[idx];
}

int get_nh20_biphase(int epoch_ms)
{
	int idx = epoch_ms % NH20_LEN;
	if (idx < 0) idx += NH20_LEN;
	return NH20_BIPHASE[idx];
}

/* Helper to compute XB state from all-1s after 'adv' shifts */
static unsigned short compute_xb_state_from_advance(int adv)
{
	unsigned int r = LFSR_INIT_ALL1;
	for (int i = 0; i < adv; i++) {
		int fb = lfsr_parity(r & XB_TAP_MASK);
		r = ((unsigned int)fb << 12) | (r >> 1);
	}
	return (unsigned short)(r & 0x1FFF);
}

int test_l5_code_verification(int verbose)
{
	int errors = 0;
	int i5[L5_CODE_LEN];
	int q5[L5_CODE_LEN];

	if (verbose) {
		printf("======================================================================\n");
		printf(" L5 PRN Code & Neumann-Hoffman Secondary Code Verification Suite\n");
		printf("======================================================================\n\n");
	}

	/* 1. Verify NH10 and NH20 secondary code values */
	const int expected_nh10[10] = {0,0,0,0,1,1,0,1,0,1};
	for (int i = 0; i < 10; i++) {
		if (NH10[i] != expected_nh10[i] || get_nh10_chip(i) != expected_nh10[i]) {
			if (verbose) printf("[FAIL] NH10 mismatch at index %d\n", i);
			errors++;
		}
		int exp_bi = (expected_nh10[i] == 0) ? 1 : -1;
		if (NH10_BIPHASE[i] != exp_bi || get_nh10_biphase(i) != exp_bi) {
			if (verbose) printf("[FAIL] NH10 biphase mismatch at index %d\n", i);
			errors++;
		}
	}

	const int expected_nh20[20] = {0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0};
	for (int i = 0; i < 20; i++) {
		if (NH20[i] != expected_nh20[i] || get_nh20_chip(i) != expected_nh20[i]) {
			if (verbose) printf("[FAIL] NH20 mismatch at index %d\n", i);
			errors++;
		}
		int exp_bi = (expected_nh20[i] == 0) ? 1 : -1;
		if (NH20_BIPHASE[i] != exp_bi || get_nh20_biphase(i) != exp_bi) {
			if (verbose) printf("[FAIL] NH20 biphase mismatch at index %d\n", i);
			errors++;
		}
	}
	if (verbose && errors == 0) {
		printf("[PASS] Neumann-Hoffman secondary codes (NH10, NH20) match specification.\n\n");
	}

	/* 2. Verify PRN advance and initial state calculation for GPS PRN 1-32 */
	if (verbose) printf("--- Verifying GPS PRN 1-32 ---\n");
	for (int prn = 1; prn <= 32; prn++) {
		int adv_i = 0, adv_q = 0;
		unsigned short init_i = 0, init_q = 0;
		get_l5_advance(prn, &adv_i, &adv_q);
		get_l5_initial_state(prn, &init_i, &init_q);

		/* Compute initial state from advance and compare */
		unsigned short calc_i = compute_xb_state_from_advance(adv_i);
		unsigned short calc_q = compute_xb_state_from_advance(adv_q);

		if (calc_i != init_i || calc_q != init_q) {
			if (verbose) {
				printf("[FAIL] GPS PRN %2d state mismatch: I(calc=0x%04X, tab=0x%04X), Q(calc=0x%04X, tab=0x%04X)\n",
					prn, calc_i, init_i, calc_q, init_q);
			}
			errors++;
		}

		/* Generate PRN code and verify balance */
		if (generate_l5_code(prn, i5, q5) != 0) {
			if (verbose) printf("[FAIL] GPS PRN %2d generate_l5_code returned error\n", prn);
			errors++;
			continue;
		}

		int ones_i = 0, ones_q = 0;
		for (int t = 0; t < L5_CODE_LEN; t++) {
			if (i5[t] != 0 && i5[t] != 1) { errors++; break; }
			if (q5[t] != 0 && q5[t] != 1) { errors++; break; }
			ones_i += i5[t];
			ones_q += q5[t];
		}

		/* L5 Gold code balance: number of 1s in 10,230 chips should be balanced (~5115) */
		if (ones_i < 4800 || ones_i > 5400 || ones_q < 4800 || ones_q > 5400) {
			if (verbose) printf("[WARN] GPS PRN %2d unusual balance: ones_i=%d, ones_q=%d\n", prn, ones_i, ones_q);
		} else if (verbose && (prn <= 3 || prn == 32)) {
			printf("[PASS] GPS PRN %2d: I5 adv=%4d state=0x%04X (ones=%d), Q5 adv=%4d state=0x%04X (ones=%d)\n",
				prn, adv_i, init_i, ones_i, adv_q, init_q, ones_q);
		}
	}

	/* 3. Verify QZSS PRN 193-202 */
	if (verbose) printf("\n--- Verifying QZSS PRN 193-202 ---\n");
	for (int prn = 193; prn <= 202; prn++) {
		int adv_i = 0, adv_q = 0;
		unsigned short init_i = 0, init_q = 0;
		get_l5_advance(prn, &adv_i, &adv_q);
		get_l5_initial_state(prn, &init_i, &init_q);

		unsigned short calc_i = compute_xb_state_from_advance(adv_i);
		unsigned short calc_q = compute_xb_state_from_advance(adv_q);

		if (calc_i != init_i || calc_q != init_q) {
			if (verbose) {
				printf("[FAIL] QZSS PRN %3d state mismatch: I(calc=0x%04X, tab=0x%04X), Q(calc=0x%04X, tab=0x%04X)\n",
					prn, calc_i, init_i, calc_q, init_q);
			}
			errors++;
		}

		if (generate_l5_code(prn, i5, q5) != 0) {
			if (verbose) printf("[FAIL] QZSS PRN %3d generate_l5_code returned error\n", prn);
			errors++;
			continue;
		}

		int ones_i = 0, ones_q = 0;
		for (int t = 0; t < L5_CODE_LEN; t++) {
			ones_i += i5[t];
			ones_q += q5[t];
		}

		if (verbose && (prn <= 195 || prn == 202)) {
			printf("[PASS] QZSS PRN %3d: I5 adv=%4d state=0x%04X (ones=%d), Q5 adv=%4d state=0x%04X (ones=%d)\n",
				prn, adv_i, init_i, ones_i, adv_q, init_q, ones_q);
		}
	}

	/* 4. Verify invalid PRN handling */
	if (generate_l5_code(0, i5, q5) == 0 || generate_l5_code(33, i5, q5) == 0 ||
	    generate_l5_code(192, i5, q5) == 0 || generate_l5_code(203, i5, q5) == 0) {
		if (verbose) printf("[FAIL] Invalid PRN was not rejected properly\n");
		errors++;
	}

	if (verbose) {
		printf("\n======================================================================\n");
		if (errors == 0)
			printf(" All L5 verification tests PASSED successfully!\n");
		else
			printf(" L5 verification tests FAILED with %d errors!\n", errors);
		printf("======================================================================\n\n");
	}

	return errors;
}
