#ifndef L5_CODE_H
#define L5_CODE_H

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief L5 primary PRN code length in chips (1 ms epoch at 10.23 Mcps) */
#define L5_CODE_LEN (10230)

/*! \brief Length of XA code sequence (chips) before short-cycle reset */
#define L5_XA_PERIOD (8190)

/*! \brief Length of XB code sequence (chips) */
#define L5_XB_PERIOD (8191)

/*! \brief Length of Neumann-Hoffman 10-bit secondary code (I5 channel, 10 ms epoch) */
#define NH10_LEN (10)

/*! \brief Length of Neumann-Hoffman 20-bit secondary code (Q5 channel, 20 ms epoch) */
#define NH20_LEN (20)

/* Secondary code sequences (binary 0 or 1) */
/* NH10: 0000110101 */
extern const int NH10[NH10_LEN];

/* NH20: 00000100110101001110 */
extern const int NH20[NH20_LEN];

/* Secondary code sequences in bi-phase (+1 / -1) mapping (0 -> +1, 1 -> -1) */
extern const int NH10_BIPHASE[NH10_LEN];
extern const int NH20_BIPHASE[NH20_LEN];

/*! \brief Check whether a PRN number is valid for L5 (GPS PRN 1-32, QZSS PRN 193-202).
 *  \param[in] prn PRN number
 *  \return 1 if valid, 0 otherwise.
 */
int is_valid_l5_prn(int prn);

/*! \brief Get the XB code advance values for a given PRN.
 *  \param[in] prn PRN number (1-32 for GPS, 193-202 for QZSS)
 *  \param[out] adv_i Output pointer for I5 XB advance (can be NULL)
 *  \param[out] adv_q Output pointer for Q5 XB advance (can be NULL)
 *  \return 0 on success, -1 if PRN is invalid.
 */
int get_l5_advance(int prn, int *adv_i, int *adv_q);

/*! \brief Get the 13-bit initial shift register state of XB for a given PRN.
 *  \param[in] prn PRN number (1-32 for GPS, 193-202 for QZSS)
 *  \param[out] state_i Output pointer for I5 initial state (13 bits, can be NULL)
 *  \param[out] state_q Output pointer for Q5 initial state (13 bits, can be NULL)
 *  \return 0 on success, -1 if PRN is invalid.
 */
int get_l5_initial_state(int prn, unsigned short *state_i, unsigned short *state_q);

/*! \brief Generate L5 I5 and/or Q5 primary PRN codes (10,230 chips, values 0 or 1).
 *  Uses 13-stage shift registers (XA and XB1/XB2) per IS-GPS-705 and IS-QZSS-PNT.
 *  \param[in] prn PRN number (1-32 for GPS, 193-202 for QZSS)
 *  \param[out] i5_code Output array of L5_CODE_LEN ints (0 or 1), can be NULL
 *  \param[out] q5_code Output array of L5_CODE_LEN ints (0 or 1), can be NULL
 *  \return 0 on success, -1 if PRN is invalid.
 */
int generate_l5_code(int prn, int *i5_code, int *q5_code);

/*! \brief Generate L5 I5 and/or Q5 primary PRN codes in bi-phase format (+1 or -1).
 *  0 maps to +1, 1 maps to -1.
 *  \param[in] prn PRN number (1-32 for GPS, 193-202 for QZSS)
 *  \param[out] i5_code Output array of L5_CODE_LEN ints (+1 or -1), can be NULL
 *  \param[out] q5_code Output array of L5_CODE_LEN ints (+1 or -1), can be NULL
 *  \return 0 on success, -1 if PRN is invalid.
 */
int generate_l5_code_biphase(int prn, int *i5_code, int *q5_code);

/*! \brief Get NH10 code bit for a given millisecond epoch (0..9).
 *  \param[in] epoch_ms 1 ms epoch index
 *  \return 0 or 1
 */
int get_nh10_chip(int epoch_ms);

/*! \brief Get NH20 code bit for a given millisecond epoch (0..19).
 *  \param[in] epoch_ms 1 ms epoch index
 *  \return 0 or 1
 */
int get_nh20_chip(int epoch_ms);

/*! \brief Get NH10 biphasic value (+1 or -1) for a given millisecond epoch.
 *  \param[in] epoch_ms 1 ms epoch index
 *  \return +1 or -1
 */
int get_nh10_biphase(int epoch_ms);

/*! \brief Get NH20 biphasic value (+1 or -1) for a given millisecond epoch.
 *  \param[in] epoch_ms 1 ms epoch index
 *  \return +1 or -1
 */
int get_nh20_biphase(int epoch_ms);

/*! \brief Self-test and verification of L5 code generation and NH secondary codes.
 *  \param[in] verbose If non-zero, print detailed verification logs to stdout.
 *  \return 0 if all tests pass, non-zero on failure.
 */
int test_l5_code_verification(int verbose);

#ifdef __cplusplus
}
#endif

#endif /* L5_CODE_H */
