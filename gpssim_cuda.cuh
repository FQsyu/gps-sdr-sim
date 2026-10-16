#ifndef GPSSIM_CUDA_CUH
#define GPSSIM_CUDA_CUH

#include "gpssim.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int prn;
	int gain;
	unsigned int carr_phase;
	int carr_phasestep;
	double code_phase;
	double f_code;
	int icode;
	int ibit;
	int iword;
	unsigned long dwrd[N_DWRD];
	int ca[CA_SEQ_LEN];
} cuda_channel_t;

typedef struct {
	int prn;
	int sat_idx;                     /*!< 0..MAX_SAT-1 (prn2sat) */
	int gain;
	unsigned int carr_phase;
	int carr_phasestep;
	double code_phase;
	double f_code;
	long long int ms_in_msg;         /*!< 0..5999 */
	signed char cnav_symbols[1200];  /*!< Current (0..599) and next (600..1199) CNAV symbols (+1 or -1) */
} cuda_channel_l5_t;

/*! \brief Initialize CUDA resources (buffers, streams, pinned memory, constant tables) */
int cuda_init_simulator(int iq_buff_size, int data_format, int band);

/*! \brief Free all CUDA resources */
void cuda_cleanup_simulator(void);

/*! \brief Queue a 0.1s block for asynchronous GPU generation */
void cuda_queue_block(
	const channel_t *chan,
	int num_channels,
	const int *gain,
	double delt,
	int iumd,
	int band
);

/*! \brief Synchronize completed block and retrieve pointer to output buffer for writing */
void cuda_sync_and_get_output(void **out_ptr, size_t *out_bytes, int iumd);

#ifdef __cplusplus
}
#endif

#endif // GPSSIM_CUDA_CUH
