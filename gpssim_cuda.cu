#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>
#include "gpssim_cuda.cuh"

// Constant memory for trigonometric tables and secondary codes
__constant__ int c_sinTable512[512];
__constant__ int c_cosTable512[512];
__constant__ int c_nh10[10];
__constant__ int c_nh20[20];

#ifdef __cplusplus
extern "C" {
#endif
extern int sinTable512[512];
extern int cosTable512[512];
#ifdef __cplusplus
}
#endif

static cudaStream_t g_streams[2];
static void *g_d_channels[2] = {NULL, NULL};
static void *g_h_channels[2] = {NULL, NULL};

static short *g_d_iq16[2] = {NULL, NULL};
static void *g_d_out[2] = {NULL, NULL};
static void *g_h_out[2] = {NULL, NULL};

// Primary PRN code memory for L5 (all 42 satellites)
static signed char *g_d_l5i = NULL;
static signed char *g_d_l5q = NULL;

static size_t g_out_bytes = 0;
static int g_iq_buff_size = 0;
static int g_data_format = 0;
static int g_initialized = 0;

__global__ void generate_iq_samples_kernel(
	const cuda_channel_t *channels,
	int num_channels,
	int iq_buff_size,
	double delt,
	short *d_iq16
)
{
	int isamp = blockIdx.x * blockDim.x + threadIdx.x;
	if (isamp >= iq_buff_size)
		return;

	int i_acc = 0;
	int q_acc = 0;

	for (int i = 0; i < num_channels; i++)
	{
		if (channels[i].prn <= 0)
			continue;

		// Carrier phase for this sample
		unsigned int carr_phase = channels[i].carr_phase + (unsigned int)isamp * (unsigned int)channels[i].carr_phasestep;
		int iTable = (carr_phase >> 16) & 0x1ff;

		// Code phase for this sample
		double total_chips = channels[i].code_phase + (double)isamp * (channels[i].f_code * delt);
		long long int total_chip_int = (long long int)floor(total_chips);
		int code_phase_idx = (int)(total_chip_int % CA_SEQ_LEN);
		if (code_phase_idx < 0)
			code_phase_idx += CA_SEQ_LEN;

		int codeCA = channels[i].ca[code_phase_idx] * 2 - 1;

		// Navigation data bit for this sample
		long long int total_codes = (total_chip_int >= 0) ? (total_chip_int / CA_SEQ_LEN) : ((total_chip_int - CA_SEQ_LEN + 1) / CA_SEQ_LEN);
		long long int total_ms = channels[i].icode + total_codes;
		long long int total_bits = total_ms / 20;
		long long int bit_idx = channels[i].ibit + total_bits;
		long long int word_idx = channels[i].iword + bit_idx / 30;
		int bit_in_word = (int)(bit_idx % 30);

		int dataBit = ((channels[i].dwrd[word_idx] >> (29 - bit_in_word)) & 0x1UL) * 2 - 1;

		int ip = dataBit * codeCA * c_cosTable512[iTable] * channels[i].gain;
		int qp = dataBit * codeCA * c_sinTable512[iTable] * channels[i].gain;

		i_acc += ip;
		q_acc += qp;
	}

	// Scaled by 2^7
	i_acc = (i_acc + 64) >> 7;
	q_acc = (q_acc + 64) >> 7;

	d_iq16[isamp * 2]     = (short)i_acc;
	d_iq16[isamp * 2 + 1] = (short)q_acc;
}

__global__ void generate_iq_samples_l5_kernel(
	const cuda_channel_l5_t *channels,
	const signed char *d_l5i_all,
	const signed char *d_l5q_all,
	int num_channels,
	int iq_buff_size,
	double delt,
	short *d_iq16
)
{
	int isamp = blockIdx.x * blockDim.x + threadIdx.x;
	if (isamp >= iq_buff_size)
		return;

	int i_acc = 0;
	int q_acc = 0;

	for (int i = 0; i < num_channels; i++)
	{
		if (channels[i].prn <= 0)
			continue;

		// Carrier phase for this sample (1176.45 MHz carrier Doppler phase)
		unsigned int carr_phase = channels[i].carr_phase + (unsigned int)isamp * (unsigned int)channels[i].carr_phasestep;
		int iTable = (carr_phase >> 16) & 0x1ff;
		int cos_val = c_cosTable512[iTable];
		int sin_val = c_sinTable512[iTable];

		// Code phase for this sample (10.23 Mcps, 10230 chips per 1 ms)
		double total_chips = channels[i].code_phase + (double)isamp * (channels[i].f_code * delt);
		long long int total_chip_int = (long long int)floor(total_chips);
		int code_phase_idx = (int)(total_chip_int % 10230);
		if (code_phase_idx < 0)
			code_phase_idx += 10230;

		// Elapsed 1 ms code epochs within this block
		long long int total_codes = (total_chip_int >= 0) ? (total_chip_int / 10230) : ((total_chip_int - 10229) / 10230);

		// PRN code lookup (+1 or -1)
		int sat_offset = channels[i].sat_idx * 10230;
		int codeL5I = (int)d_l5i_all[sat_offset + code_phase_idx];
		int codeL5Q = (int)d_l5q_all[sat_offset + code_phase_idx];

		// Millisecond within the 6-second CNAV message epoch
		long long int total_ms = channels[i].ms_in_msg + total_codes;

		// CNAV symbol index (10 ms per symbol, 100 sps): 0..1199 in 2-message window
		int sym_idx = (int)(total_ms / 10);
		if (sym_idx < 0) sym_idx = 0;
		else if (sym_idx >= 1200) sym_idx = 1199;
		int cnav_symbol = (int)channels[i].cnav_symbols[sym_idx]; // +1 or -1

		// NH10 secondary code (10 ms period, 10 chips: 1 chip / ms)
		int ms_in_symbol = (int)(total_ms % 10);
		if (ms_in_symbol < 0) ms_in_symbol += 10;
		int nh10 = c_nh10[ms_in_symbol]; // +1 or -1

		// NH20 secondary code (20 ms period, 20 chips: 1 chip / ms)
		int ms_epoch = (int)(total_ms % 20);
		if (ms_epoch < 0) ms_epoch += 20;
		int nh20 = c_nh20[ms_epoch]; // +1 or -1

		// I5 component: CNAV symbol * NH10 * I5 PRN
		int s_i = cnav_symbol * nh10 * codeL5I;

		// Q5 component: NH20 * Q5 PRN (pilot, no data)
		int s_q = nh20 * codeL5Q;

		// QPSK carrier modulation: (s_i * cos - s_q * sin) + j (s_i * sin + s_q * cos)
		int ip = (s_i * cos_val - s_q * sin_val) * channels[i].gain;
		int qp = (s_i * sin_val + s_q * cos_val) * channels[i].gain;

		i_acc += ip;
		q_acc += qp;
	}

	// Scaled by 2^7
	i_acc = (i_acc + 64) >> 7;
	q_acc = (q_acc + 64) >> 7;

	d_iq16[isamp * 2]     = (short)i_acc;
	d_iq16[isamp * 2 + 1] = (short)q_acc;
}

__global__ void quantize_sc08_kernel(const short *d_in, signed char *d_out, int total_elements)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < total_elements)
	{
		d_out[idx] = (signed char)(d_in[idx] >> 4);
	}
}

__global__ void quantize_sc01_kernel(const short *d_in, unsigned char *d_out, int num_bytes)
{
	int byte_idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (byte_idx < num_bytes)
	{
		unsigned char b = 0;
		int base = byte_idx * 8;
		#pragma unroll
		for (int bit = 0; bit < 8; bit++)
		{
			if (d_in[base + bit] > 0)
				b |= (1 << (7 - bit));
		}
		d_out[byte_idx] = b;
	}
}

int cuda_init_simulator(int iq_buff_size, int data_format, int band)
{
	if (g_initialized)
		return 0;

	int device_count = 0;
	cudaError_t err = cudaGetDeviceCount(&device_count);
	if (err != cudaSuccess || device_count == 0)
	{
		fprintf(stderr, "CUDA error: No CUDA-capable devices found (%s)\n", cudaGetErrorString(err));
		return -1;
	}

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	fprintf(stderr, "Using CUDA Device: %s (Compute %d.%d)\n", prop.name, prop.major, prop.minor);

	// Copy trigonometric tables to constant memory
	err = cudaMemcpyToSymbol(c_sinTable512, sinTable512, sizeof(int) * 512);
	if (err != cudaSuccess)
		fprintf(stderr, "CUDA error copying sinTable: %s\n", cudaGetErrorString(err));
	err = cudaMemcpyToSymbol(c_cosTable512, cosTable512, sizeof(int) * 512);
	if (err != cudaSuccess)
		fprintf(stderr, "CUDA error copying cosTable: %s\n", cudaGetErrorString(err));

	g_iq_buff_size = iq_buff_size;
	g_data_format = data_format;

	if (data_format == SC01)
		g_out_bytes = (size_t)iq_buff_size / 4;
	else if (data_format == SC08)
		g_out_bytes = (size_t)iq_buff_size * 2 * sizeof(signed char);
	else
		g_out_bytes = (size_t)iq_buff_size * 2 * sizeof(short);

	size_t iq16_bytes = (size_t)iq_buff_size * 2 * sizeof(short);
	size_t max_single_chan = sizeof(cuda_channel_t) > sizeof(cuda_channel_l5_t) ? sizeof(cuda_channel_t) : sizeof(cuda_channel_l5_t);
	size_t chan_bytes = (size_t)MAX_CHAN * max_single_chan;

	for (int b = 0; b < 2; b++)
	{
		cudaStreamCreateWithFlags(&g_streams[b], cudaStreamNonBlocking);

		cudaMalloc(&g_d_channels[b], chan_bytes);
		cudaHostAlloc(&g_h_channels[b], chan_bytes, cudaHostAllocDefault);

		cudaMalloc(&g_d_iq16[b], iq16_bytes);

		if (data_format == SC16)
		{
			g_d_out[b] = g_d_iq16[b];
		}
		else
		{
			cudaMalloc(&g_d_out[b], g_out_bytes);
		}

		cudaHostAlloc(&g_h_out[b], g_out_bytes, cudaHostAllocDefault);
	}

	if (band == BAND_L5)
	{
		// Copy NH10 and NH20 secondary codes to constant memory
		err = cudaMemcpyToSymbol(c_nh10, NH10_BIPHASE, sizeof(int) * 10);
		if (err != cudaSuccess)
			fprintf(stderr, "CUDA error copying c_nh10: %s\n", cudaGetErrorString(err));
		err = cudaMemcpyToSymbol(c_nh20, NH20_BIPHASE, sizeof(int) * 20);
		if (err != cudaSuccess)
			fprintf(stderr, "CUDA error copying c_nh20: %s\n", cudaGetErrorString(err));

		// Allocate and upload all 42 satellites' primary PRN codes
		size_t l5_codes_bytes = (size_t)MAX_SAT * CA_SEQ_LEN_L5 * sizeof(signed char);
		cudaMalloc(&g_d_l5i, l5_codes_bytes);
		cudaMalloc(&g_d_l5q, l5_codes_bytes);

		signed char *h_l5i = (signed char *)malloc(l5_codes_bytes);
		signed char *h_l5q = (signed char *)malloc(l5_codes_bytes);
		int temp_i[CA_SEQ_LEN_L5], temp_q[CA_SEQ_LEN_L5];

		for (int sv = 0; sv < MAX_SAT; sv++)
		{
			int prn = sat2prn(sv);
			if (is_valid_l5_prn(prn))
			{
				generate_l5_code_biphase(prn, temp_i, temp_q);
				for (int k = 0; k < CA_SEQ_LEN_L5; k++)
				{
					h_l5i[sv * CA_SEQ_LEN_L5 + k] = (signed char)temp_i[k];
					h_l5q[sv * CA_SEQ_LEN_L5 + k] = (signed char)temp_q[k];
				}
			}
			else
			{
				memset(&h_l5i[sv * CA_SEQ_LEN_L5], 0, CA_SEQ_LEN_L5);
				memset(&h_l5q[sv * CA_SEQ_LEN_L5], 0, CA_SEQ_LEN_L5);
			}
		}

		cudaMemcpy(g_d_l5i, h_l5i, l5_codes_bytes, cudaMemcpyHostToDevice);
		cudaMemcpy(g_d_l5q, h_l5q, l5_codes_bytes, cudaMemcpyHostToDevice);

		free(h_l5i);
		free(h_l5q);
	}

	g_initialized = 1;
	return 0;
}

void cuda_cleanup_simulator(void)
{
	if (!g_initialized)
		return;

	if (g_d_l5i) { cudaFree(g_d_l5i); g_d_l5i = NULL; }
	if (g_d_l5q) { cudaFree(g_d_l5q); g_d_l5q = NULL; }

	for (int b = 0; b < 2; b++)
	{
		cudaStreamSynchronize(g_streams[b]);
		cudaStreamDestroy(g_streams[b]);

		if (g_d_channels[b]) cudaFree(g_d_channels[b]);
		if (g_h_channels[b]) cudaFreeHost(g_h_channels[b]);

		if (g_d_iq16[b]) cudaFree(g_d_iq16[b]);
		if (g_data_format != SC16 && g_d_out[b]) cudaFree(g_d_out[b]);

		if (g_h_out[b]) cudaFreeHost(g_h_out[b]);

		g_d_channels[b] = NULL;
		g_h_channels[b] = NULL;
		g_d_iq16[b] = NULL;
		g_d_out[b] = NULL;
		g_h_out[b] = NULL;
	}

	g_initialized = 0;
}

void cuda_queue_block(
	const channel_t *chan,
	int num_channels,
	const int *gain,
	double delt,
	int iumd,
	int band
)
{
	int idx = (iumd - 1) % 2;
	cudaStream_t stream = g_streams[idx];

	int threads = 256;
	int blocks = (g_iq_buff_size + threads - 1) / threads;

	if (band == BAND_L5)
	{
		cuda_channel_l5_t *h_chans_l5 = (cuda_channel_l5_t *)g_h_channels[idx];
		for (int i = 0; i < num_channels && i < MAX_CHAN; i++)
		{
			h_chans_l5[i].prn = chan[i].prn;
			h_chans_l5[i].sat_idx = prn2sat(chan[i].prn);
			h_chans_l5[i].gain = gain[i];
			h_chans_l5[i].carr_phase = chan[i].carr_phase;
			h_chans_l5[i].carr_phasestep = chan[i].carr_phasestep;
			h_chans_l5[i].code_phase = chan[i].code_phase;
			h_chans_l5[i].f_code = chan[i].f_code;
			h_chans_l5[i].ms_in_msg = chan[i].ms_in_msg;
			if (chan[i].prn > 0)
			{
				cnav_channel_get_2msg_symbols(&chan[i].cnav, h_chans_l5[i].cnav_symbols);
			}
			else
			{
				memset(h_chans_l5[i].cnav_symbols, 0, sizeof(h_chans_l5[i].cnav_symbols));
			}
		}

		size_t chan_bytes = (size_t)num_channels * sizeof(cuda_channel_l5_t);
		cudaMemcpyAsync(g_d_channels[idx], h_chans_l5, chan_bytes, cudaMemcpyHostToDevice, stream);

		generate_iq_samples_l5_kernel<<<blocks, threads, 0, stream>>>(
			(const cuda_channel_l5_t *)g_d_channels[idx],
			g_d_l5i,
			g_d_l5q,
			num_channels,
			g_iq_buff_size,
			delt,
			g_d_iq16[idx]
		);
	}
	else
	{
		cuda_channel_t *h_chans = (cuda_channel_t *)g_h_channels[idx];
		for (int i = 0; i < num_channels && i < MAX_CHAN; i++)
		{
			h_chans[i].prn = chan[i].prn;
			h_chans[i].gain = gain[i];
			h_chans[i].carr_phase = chan[i].carr_phase;
			h_chans[i].carr_phasestep = chan[i].carr_phasestep;
			h_chans[i].code_phase = chan[i].code_phase;
			h_chans[i].f_code = chan[i].f_code;
			h_chans[i].icode = chan[i].icode;
			h_chans[i].ibit = chan[i].ibit;
			h_chans[i].iword = chan[i].iword;
			memcpy(h_chans[i].dwrd, chan[i].dwrd, sizeof(unsigned long) * N_DWRD);
			memcpy(h_chans[i].ca, chan[i].ca, sizeof(int) * CA_SEQ_LEN);
		}

		size_t chan_bytes = (size_t)num_channels * sizeof(cuda_channel_t);
		cudaMemcpyAsync(g_d_channels[idx], h_chans, chan_bytes, cudaMemcpyHostToDevice, stream);

		generate_iq_samples_kernel<<<blocks, threads, 0, stream>>>(
			(const cuda_channel_t *)g_d_channels[idx],
			num_channels,
			g_iq_buff_size,
			delt,
			g_d_iq16[idx]
		);
	}

	cudaError_t kerr = cudaGetLastError();
	if (kerr != cudaSuccess)
		fprintf(stderr, "CUDA kernel launch error: %s\n", cudaGetErrorString(kerr));

	if (g_data_format == SC08)
	{
		int total_elem = g_iq_buff_size * 2;
		int q_blocks = (total_elem + threads - 1) / threads;
		quantize_sc08_kernel<<<q_blocks, threads, 0, stream>>>(
			g_d_iq16[idx],
			(signed char *)g_d_out[idx],
			total_elem
		);
	}
	else if (g_data_format == SC01)
	{
		int num_bytes = g_iq_buff_size / 4;
		int q_blocks = (num_bytes + threads - 1) / threads;
		quantize_sc01_kernel<<<q_blocks, threads, 0, stream>>>(
			g_d_iq16[idx],
			(unsigned char *)g_d_out[idx],
			num_bytes
		);
	}

	cudaMemcpyAsync(g_h_out[idx], g_d_out[idx], g_out_bytes, cudaMemcpyDeviceToHost, stream);
}

void cuda_sync_and_get_output(void **out_ptr, size_t *out_bytes, int iumd)
{
	int idx = (iumd - 1) % 2;
	cudaStreamSynchronize(g_streams[idx]);
	*out_ptr = g_h_out[idx];
	*out_bytes = g_out_bytes;
}
