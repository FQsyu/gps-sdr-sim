#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "l5_code.h"

int main(int argc, char **argv)
{
	printf("Running GPS / QZSS L5 PRN Code & Neumann-Hoffman Unit Tests...\n\n");

	/* 1. Run complete verification suite */
	int res = test_l5_code_verification(1);
	if (res != 0) {
		fprintf(stderr, "test_l5_code_verification failed with %d errors!\n", res);
		return 1;
	}

	/* 2. Specific check: IS-GPS-705 PRN 1 initial state & first chips */
	unsigned short state_i = 0, state_q = 0;
	int adv_i = 0, adv_q = 0;
	get_l5_advance(1, &adv_i, &adv_q);
	get_l5_initial_state(1, &state_i, &state_q);
	assert(adv_i == 266);
	assert(adv_q == 1701);
	assert(state_i == 0x0AE4); /* 0101011100100 in binary */
	assert(state_q == 0x12CC); /* 1001011001100 in binary */

	/* 3. Specific check: IS-QZSS-PNT PRN 193 (QZS-1) initial state */
	get_l5_advance(193, &adv_i, &adv_q);
	get_l5_initial_state(193, &state_i, &state_q);
	assert(adv_i == 5836);
	assert(adv_q == 4757);
	assert(state_i == 0x0C2E); /* 0110000101110 in binary */
	assert(state_q == 0x1387); /* 1001110000111 in binary */

	/* 4. Check bi-phase mapping */
	int i5_bin[L5_CODE_LEN], q5_bin[L5_CODE_LEN];
	int i5_bi[L5_CODE_LEN], q5_bi[L5_CODE_LEN];
	generate_l5_code(1, i5_bin, q5_bin);
	generate_l5_code_biphase(1, i5_bi, q5_bi);
	for (int i = 0; i < L5_CODE_LEN; i++) {
		assert(i5_bi[i] == ((i5_bin[i] == 0) ? 1 : -1));
		assert(q5_bi[i] == ((q5_bin[i] == 0) ? 1 : -1));
	}
	printf("[PASS] Bi-phase (+1/-1) mapping verified for PRN 1.\n");

	/* 5. Check null pointer safety */
	assert(generate_l5_code(1, i5_bin, NULL) == 0);
	assert(generate_l5_code(1, NULL, q5_bin) == 0);
	assert(generate_l5_code(1, NULL, NULL) == 0);
	printf("[PASS] Null-pointer safety verified.\n");

	/* 6. Verify NH secondary code modulo access */
	for (int epoch = 0; epoch < 100; epoch++) {
		assert(get_nh10_chip(epoch) == NH10[epoch % 10]);
		assert(get_nh20_chip(epoch) == NH20[epoch % 20]);
		assert(get_nh10_biphase(epoch) == NH10_BIPHASE[epoch % 10]);
		assert(get_nh20_biphase(epoch) == NH20_BIPHASE[epoch % 20]);
	}
	printf("[PASS] NH10 and NH20 multi-epoch rollover verified.\n\n");

	printf("ALL UNIT TESTS PASSED SUCCESSFULLY (100%%)!\n");
	return 0;
}
