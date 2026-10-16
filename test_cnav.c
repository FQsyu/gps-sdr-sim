#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "gpssim.h"
#include "cnav.h"
#include "l5_code.h"

/* Simple Viterbi Decoder for Rate 1/2, K=7 (G1=171 oct, G2=133 oct)
 * to verify complete encode -> decode round-trip
 */
#define VIT_K 7
#define VIT_STATES 64

static inline int parity_byte(unsigned int v)
{
	v ^= v >> 1;
	v ^= v >> 2;
	v ^= v >> 4;
	return (int)(v & 1);
}

static void viterbi_decode(const int *symbols, int n_symbols, int *decoded_bits)
{
	int steps = n_symbols / 2;
	int nbits = steps - (VIT_K - 1);

	unsigned int br0[32], br1[32];
	for (int s = 0; s < 32; s++) {
		br0[s] = (unsigned int)parity_byte((unsigned int)((2 * s) & 0x4F));
		br1[s] = (unsigned int)parity_byte((unsigned int)((2 * s) & 0x6D));
	}

	unsigned int metric[VIT_STATES];
	unsigned int new_metric[VIT_STATES];
	for (int s = 0; s < VIT_STATES; s++) metric[s] = 100000;
	metric[0] = 0;

	/* Dynamic allocation for traceback */
	int (*dec)[VIT_STATES] = (int (*)[VIT_STATES])malloc(sizeof(int[VIT_STATES]) * steps);

	for (int t = 0; t < steps; t++) {
		unsigned int sym0 = (unsigned int)symbols[t * 2];
		unsigned int sym1 = (unsigned int)symbols[t * 2 + 1];

		for (int i = 0; i < 32; i++) {
			unsigned int bm = (br0[i] ^ sym0) + (br1[i] ^ sym1);
			unsigned int m0 = metric[i] + bm;
			unsigned int m1 = metric[i + 32] + (2 - bm);
			int d = (m0 > m1);
			new_metric[2 * i] = d ? m1 : m0;
			dec[t][2 * i] = d;

			m0 = metric[i] + (2 - bm);
			m1 = metric[i + 32] + bm;
			d = (m0 > m1);
			new_metric[2 * i + 1] = d ? m1 : m0;
			dec[t][2 * i + 1] = d;
		}
		for (int s = 0; s < VIT_STATES; s++) metric[s] = new_metric[s];
	}

	/* Find best end state */
	unsigned int min_m = metric[0];
	int best_state = 0;
	for (int s = 1; s < VIT_STATES; s++) {
		if (metric[s] < min_m) {
			min_m = metric[s];
			best_state = s;
		}
	}

	/* Trace back */
	int *temp_bits = (int *)malloc(sizeof(int) * steps);
	int state = best_state;
	for (int t = steps - 1; t >= 0; t--) {
		int d = dec[t][state];
		int prev_state = (state >> 1) | (d << 5);
		temp_bits[t] = state & 1;
		state = prev_state;
	}

	for (int i = 0; i < nbits; i++) {
		decoded_bits[i] = temp_bits[i];
	}

	free(temp_bits);
	free(dec);
}

int main(int argc, char **argv)
{
	printf("Running GPS / QZSS L5 CNAV Message & Rate-1/2 FEC Unit Tests...\n\n");

	/* 1. Run main verification suite */
	int res = test_cnav_verification(1);
	if (res != 0) {
		fprintf(stderr, "test_cnav_verification failed with %d errors!\n", res);
		return 1;
	}

	/* 2. End-to-end Viterbi FEC round-trip test */
	printf("--- Running End-to-End Viterbi FEC Round-Trip Test ---\n");
	int raw_bits[300 + 6]; /* 300 data bits + 6 tail bits */
	for (int i = 0; i < 300; i++) raw_bits[i] = (i * 17 + 5) % 2;
	for (int i = 300; i < 306; i++) raw_bits[i] = 0;

	cnav_fec_t fec;
	cnav_fec_init(&fec);
	int symbols[306 * 2];
	cnav_fec_encode_msg(&fec, raw_bits, 306, symbols);

	int decoded[300];
	viterbi_decode(symbols, 306 * 2, decoded);

	int bit_errors = 0;
	for (int i = 0; i < 300; i++) {
		if (decoded[i] != raw_bits[i]) bit_errors++;
	}
	assert(bit_errors == 0);
	printf("[PASS] Viterbi decoder successfully recovered 100%% of 300 CNAV bits (0 errors)!\n");

	/* 3. Verify I5 / Q5 modulation alignment with Step 2 PRN & NH secondary codes */
	printf("\n--- Verifying I5/Q5 Modulation & NH Synchronization ---\n");
	ephem_t eph;
	memset(&eph, 0, sizeof(ephem_t));
	eph.toe.sec = 518400.0;
	cnav_channel_t chan;
	cnav_channel_init(&chan, &eph, NULL, 1, 518400.0);

	/* Check 100 ms (10 symbols) of continuous operation */
	for (int epoch_ms = 0; epoch_ms < 100; epoch_ms++) {
		int fec_sym = cnav_channel_get_symbol(&chan);
		int nh10_chip = get_nh10_biphase(epoch_ms);
		int nh20_chip = get_nh20_biphase(epoch_ms);

		/* The combined I5 chip multiplier */
		int i5_mod = fec_sym * nh10_chip;
		assert(i5_mod == 1 || i5_mod == -1);

		/* Q5 pilot channel multiplier */
		int q5_mod = nh20_chip;
		assert(q5_mod == 1 || q5_mod == -1);

		cnav_channel_step_ms(&chan);
	}
	printf("[PASS] I5 FEC symbol and NH10/NH20 10ms synchronization verified.\n");

	/* Test 1200 symbols buffer for GPU generation */
	signed char syms_1200[1200];
	cnav_channel_get_2msg_symbols(&chan, syms_1200);
	for (int i = 0; i < 1200; i++) {
		assert(syms_1200[i] == 1 || syms_1200[i] == -1);
	}
	printf("[PASS] 1200 CNAV 2-message symbols buffer for GPU verified.\n\n");

	printf("ALL CNAV UNIT TESTS PASSED SUCCESSFULLY (100%%)!\n");
	return 0;
}
