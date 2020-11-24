/*************************************************/
/* Author: Chao-wang Huang                       */
/* Date: Monday, February 04, 2008               */
/* An (2,1,6) Convolutional code is simulated    */
/* Decoding algorithm: SOVA decoding             */
/* Interleaver: S-random Interleaver         	 */
/* Modulation: QPSK with OFDM					 */
/* ITU Vehicular A Channel                       */
/* SISO case									 */
/*************************************************/

#include "stdafx.h"
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <conio.h>
#include <limits.h>
#include <float.h>
#include <cstdlib>
#include <iostream>
using namespace std;
#include <fstream>
#include <cstdio>

void AWGN_noise(float, double, double *);
void JakesFading(double, double, double, int, double *);
void Multipath_fade_pattern(double **, int);
int Error_count(int, int);
void Trellis_diagram(int, int, int, int *);
void CC_Encoder(int *, int *, int, int, int *, int, int);
void Viterbi_dec_SOVA(double *, int *, double *, int, int, int);
void dfour1(double *, unsigned long, int);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
const int n_fft = 1024;
#define Pi 		3.14159265358979
const int n = 2;
const int k = 1;
const int m = 6;
#define K		1024		  			// packet length of information bit stream
#define N		K*(n/k)	 				// packet length of coded bit stream
#define num_packet 10				// number of packets simulated
#define N_state	 64
const int num_state = (int)pow(2.0,(double)m);			// number of states
const int num_in_sym = (int)pow(2.0,(double)k);		// number of input symbols
//const int gp[2] = {171,133};			// Generator polynomial of CC given in Octal
int gp[2] = {79,109};           		// Generator polynomial of CC given in Decimal and in reverse order
const int Mc = 2;						// modulation order (2 for QPSK)
const int CP_length = n_fft/8;			// Length of Cyclic Prefix
const double sample_rate = 11.2e6;		/* Sampling Frequency */
// Power Delay Profile of ITU Vehicular A Channel
const int path_num = 6;
const int path_delay[6] = {0,3,8,12,19,28};	// Delay (sample) {0,310,710,1090,1730,2510} (ns)
const double path_gain[6] = {0.485,0.3853,0.0611,0.0485,0.0153,0.0049};	// Power Gain (linear scale)
const int delay_spread = 28;			// Maximum path delay (samples)
const double vc = 0.0;					/* speed of vehicle in km/hr */
const double C = 3.0e8;					/* light speed */
const double fc = 2.5e9;				/* carrier frequency */
const double OFDM_sym_duration = 102.86e-6;		/* OFDMA Symbol Duration with Cyclic Prefix */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)
const int L = 10000;					// Channel Realizations
double *ts;

struct trel_diag						// Data structure of trellis diagram for each branch
		{
      		int from;					// from_state of trellis diagram
			int to;						// to_state of trellis diagram6
			int in;						// input data bit of trellis diagram
			int out[n];					// output codeword symbol of trellis diagram
		};
struct trel_diag Trellis[N_state][2];	// Trellis[num_state][num_in_sym]

struct surv        			// Data structure of survival path
   		{
			double metric;		// Path metric
            int data_in[K];		// input bit stream, K: packet length, global
            int code_bit[N];	// output code bit stream, N: packet length, global
            int state[K];		// state transition sequence, K: packet length, global
        };
	
	struct surv survival_path, survival_temp[N_state];
   								// Survival_path[N_state], N_state: number of states of CC code, global

int main(void)
{
	time_t  t, start, end;
	int i, j, l, q, p, *data_bit, *coded_bit, *Hk, *delay, err_sum, bler_sum;
	int *S_random, *de_inter, *interleaved, *err_count, *bler_count, BLER;
	double snr, Eb_No, noise_pwr, noise[2], *Yk, /*err_rate, *bler_rate*/ err_rate, bler_rate;
	double Ps, *QPSK_sym, *Sk, *Dk, *LLR, *TX_signal_I, *TX_signal_Q, *RX_signal_I, *RX_signal_Q;
	double **ISI_I, **ISI_Q, **fade, **path_I, **path_Q, *H_f, *gain;
	FILE *ber, *records, *s_inter;
	
	start = time(NULL);
	printf("BLER Performance of coded-OFDM in ITU Vehicular A Channel\n");
	printf("Coding Scheme: (2,1,6) Convolutional Code\n");
	printf("Generator polynomials are {171,133} in Octal\n");
	printf("Minimum free distance = 10\n");
	printf("Decoder: SOVA Decoder\n");
	printf("Interleaver: S-random Interleaver\n");
	printf("Speed of the vehicle = %f (km/h)\n", vc);
	printf("Carrier Frequency = %e (Hz)\n", fc);
	printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
	printf("f_d * t = %f\n", Doppler*OFDM_sym_duration);
	printf("Maximum number of bits of simulation = %d\n", (K-m)*num_packet);
	printf("This program is running. Don't close, please!\n\n");

	fopen_s(&records, "Records_COFDM_QPSK.log", "a");
	fprintf(records, "BLER Performance of coded-OFDM in ITU Vehicular A Channel\n");
	fprintf(records, "Coding Scheme: (2,1,6) Convolutional Code\n");
	fprintf(records, "Generator polynomials are {171,133} in Octal\n");
	fprintf(records, "Minimum free distance = 10\n");
	fprintf(records, "Decoder: SOVA Decoder\n");
	fprintf(records, "Interleaver: S-random Interleaver\n");
	fprintf(records, "Speed of the vehicle = %f (km/h)\n", vc);
	fprintf(records, "Carrier Frequency = %e (Hz)\n", fc);
	fprintf(records, "Maximum Doppler Frequency = %f (Hz)\n", Doppler);
	fprintf(records, "f_d * t = %f\n", Doppler*OFDM_sym_duration);
	fprintf(records, "Maximum number of bits of simulation = %d\n", (K-m)*num_packet);
	fprintf(records, "Eb/No      BER              BLER\n");
	fflush(records);
	
	data_bit = new int[K];
	coded_bit = new int[N];
	Hk = new int[K];              		// Hard Decision
	Sk = new double[N];					// Soft Decision
	Dk = new double[N];					// De-Interleaved
	LLR = new double[N];				// Log-likelihood Ratio
	Yk = new double[2*n_fft+1];			// Received signal
	S_random = new int[N];
	de_inter = new int[N];
	interleaved = new int[N];
	QPSK_sym = new double[2*n_fft+1];
	TX_signal_I = new double[n_fft+CP_length];
	TX_signal_Q = new double[n_fft+CP_length];
	RX_signal_I = new double[n_fft+CP_length];
	RX_signal_Q = new double[n_fft+CP_length];
	path_I = new double*[path_num];
	path_Q = new double*[path_num];
	for(i=0; i<path_num; i++)
	{
		path_I[i] = new double[n_fft+CP_length+delay_spread];
		path_Q[i] = new double[n_fft+CP_length+delay_spread];
	}
	H_f = new double[2*n_fft+1];
	ts = new double[path_num];
	delay = new int[path_num];
	ISI_I = new double*[path_num];
	for(i=0; i<path_num; i++)
		ISI_I[i] = new double[delay_spread];
	ISI_Q = new double*[path_num];
	for(i=0; i<path_num; i++)
		ISI_Q[i] = new double[delay_spread];
	gain = new double[path_num];
	fade = new double*[path_num];
	for(i=0; i<path_num; i++)
		fade[i] = new double[2];
	err_count = new int[L];
	bler_count = new int[L];
	//err_rate = new double[L];
	//bler_rate = new double[L];
	
	// Channel Weighting Gain and Channel path delay
	for(i=0; i<path_num; i++)
	{
		delay[i] = path_delay[i];
		gain[i] = sqrt(path_gain[i]);
	}
	
	srand((unsigned) time(&t));

	Trellis_diagram(n, k, m, &gp[0]);
	
/**************************************************/
/* S-random interleaver and de-interleaver (S=32) */
/**************************************************/
	fopen_s(&s_inter, "s_random.txt", "r");
	for(i=0; i<N; i++)		// Interleaver
   		fscanf_s(s_inter, "%d %d", &de_inter[i], &S_random[i]);
	for(i=0; i<N; i++)		// De-interleaver
   		de_inter[S_random[i]] = i;
	fclose(s_inter);

/************************/
/* main simulation loop */
/************************/ 
	fopen_s(&ber, "ber_COFDM_SOVA.log", "w");
	for(snr=2; snr<=18; snr+=4)
	{
		for(l=0; l<L; l++)
			err_count[l] = bler_count[l] = 0;
   		// noise power calculation
		Eb_No = (double)snr;
		Ps = 1.0;
		noise_pwr = 0.5*Ps/(n_fft*(k/(float)n)*pow(10.0, Eb_No/10.0));	// QPSK, Nyquist filter assumption
		printf("Eb_No = %.1f, ", Eb_No);

		// Initial time of Rayleigh fading pattern
		for(l=0; l<L; l++)
		{
			for(i=0; i<path_num; i++)
				ts[i] = 100.0 + l*100.0 + i*10;
		
			// Multi-path Rayleigh fading pattern
			Multipath_fade_pattern(fade, path_num);

			// Initialize ISI buffer
			for(i=1; i<path_num; i++)		// path index
				for(j=0; j<delay[i]; j++)	// path delay
				{
					ISI_I[i][j] = 0.0;
					ISI_Q[i][j] = 0.0;
				}
			
			p = 0;
			do
			{
				// Generate random information bit stream
				for(i=0; i<K-m; i++)
					if(rand()/(float)RAND_MAX>=0.5)
						data_bit[i] = 1;
					else
						data_bit[i] = 0;

				for(i=K-m; i<K; i++)
					data_bit[i] = 0;

				// Convolutional Encoder 
				CC_Encoder(data_bit, coded_bit, n, m, &gp[0], K, num_state);

  				// Interleaving
				for(i=0; i<N; i++)
           			interleaved[S_random[i]] = coded_bit[i];
    
/*****************************************************/
/* QPSK mapping, IFFT and ITU Vehicular A Channel    */
/*****************************************************/
	   			// QPSK mapping
				for(q=0; q<n_fft; q++)
   				{
					QPSK_sym[2*q+1] = (2*interleaved[2*q]-1)/sqrt(2.0);
					QPSK_sym[2*q+2] = (2*interleaved[2*q+1]-1)/sqrt(2.0);
				}

				// IFFT
				dfour1(&QPSK_sym[0], 1024, 1);  // IFFT for transmitted signal must be multiplied by n_fft^(-1)

				for(j=0; j<n_fft; j++)
				{
         			QPSK_sym[2*j+1] = (QPSK_sym[2*j+1])/(double)n_fft;
					QPSK_sym[2*j+2] = (QPSK_sym[2*j+2])/(double)n_fft;
				}

				// Cyclic Prefix
				for(i=0; i<CP_length; i++)
				{
					TX_signal_I[i] = QPSK_sym[2*(i+n_fft-CP_length)+1];
					TX_signal_Q[i] = QPSK_sym[2*(i+n_fft-CP_length)+2];
				}

				for(i=CP_length; i<n_fft+CP_length; i++)
				{
					TX_signal_I[i] = QPSK_sym[2*(i-CP_length)+1];
					TX_signal_Q[i] = QPSK_sym[2*(i-CP_length)+2];
				}

				
				// Multi-path Rayleigh fading pattern
				//Multipath_fade_pattern(fade, path_num);
				
				// Multipath Signal
				for(i=0; i<path_num; i++)
					for(j=delay[i]; j<n_fft+CP_length+delay[i]; j++)
					{
						path_I[i][j] = gain[i] * (TX_signal_I[j-delay[i]]*fade[i][0] - TX_signal_Q[j-delay[i]]*fade[i][1]);
						path_Q[i][j] = gain[i] * (TX_signal_I[j-delay[i]]*fade[i][1] + TX_signal_Q[j-delay[i]]*fade[i][0]);
					}
			
				// ISI Records
				for(i=1; i<path_num; i++)
					for(j=0; j<delay[i]; j++)
					{
						path_I[i][j] = ISI_I[i][j];
						path_Q[i][j] = ISI_Q[i][j];
						ISI_I[i][j] = path_I[i][j+n_fft+CP_length];
						ISI_Q[i][j] = path_Q[i][j+n_fft+CP_length];
					}

				for(j=0; j<n_fft+CP_length; j++)
				{
					RX_signal_I[j] = 0.0;
					RX_signal_Q[j] = 0.0;

					// Multipath Channel
					for(i=0; i<path_num; i++)
					{
						RX_signal_I[j] = RX_signal_I[j] + path_I[i][j];
						RX_signal_Q[j] = RX_signal_Q[j] + path_Q[i][j];
					}
				
					// AWGN noise
					AWGN_noise(0.0, noise_pwr, &noise[0]);
		  	   		RX_signal_I[j] += noise[0];
			 		RX_signal_Q[j] += noise[1];
				}
				
/******************************************************************************/
/*                FFT, Equalization and QPSK De-mapping                       */
/******************************************************************************/
				// Remove Cyclic Prefix
				for(i=0; i<n_fft; i++)
				{
					Yk[2*i+1] = RX_signal_I[i+CP_length];
					Yk[2*i+2] = RX_signal_Q[i+CP_length];
				}
			
				// FFT
				dfour1(&Yk[0], 1024, -1);

				/* Channel Frequency Response */
				for(i=0; i<=2*n_fft; i++)
					H_f[i] = 0.0;

				// Time Domain Impulse Response
				for(j=0; j<path_num; j++)	
				{
					H_f[2*delay[j]+1] = gain[j] * fade[j][0];
					H_f[2*delay[j]+2] = gain[j] * fade[j][1];
				}

				// FFT (Frequency Response)
				dfour1(&H_f[0], 1024, -1);

				// Frequency Domain Equalization (Matched Filter)
				for(i=0; i<n_fft; i++)
				{
					Sk[2*i] = (Yk[2*i+1]*H_f[2*i+1] + Yk[2*i+2]*H_f[2*i+2]);
			    	Sk[2*i+1] = (Yk[2*i+2]*H_f[2*i+1] - Yk[2*i+1]*H_f[2*i+2]);
				}
/*
			// Frequency Domain Equalization (ZF)
			for(i=0; i<n_fft; i++)
			{
				if((pow(H_f[2*i+1],2.0)+pow(H_f[2*i+2],2.0)) == 0.0)
					Sk[2*i] = Sk[2*i+1] = 0.0;
				else
				{
					Sk[2*i] = (Yk[2*i+1]*H_f[2*i+1]+Yk[2*i+2]*H_f[2*i+2]) / (pow(H_f[2*i+1],2.0)+pow(H_f[2*i+2],2.0));
            		Sk[2*i+1] = (Yk[2*i+2]*H_f[2*i+1]-Yk[2*i+1]*H_f[2*i+2]) / (pow(H_f[2*i+1],2.0)+pow(H_f[2*i+2],2.0));
				}
			}
*/
				// De-Interleaving
				for(i=0; i<N; i++)
					Dk[de_inter[i]] = Sk[i];

				// Soft out Viterbi decoding
	            Viterbi_dec_SOVA(Dk, Hk, LLR, K, num_state, num_in_sym);
				
				// Error Count
				BLER = 0;
				for(i=0; i<K-m; i++)	// Bit error count
	            {
					err_count[l] += Error_count(data_bit[i], Hk[i]);
				
					if(Error_count(data_bit[i], Hk[i]) == 1) // Codeword error count
						BLER = 1;
			    }

				if(BLER == 1)
					bler_count[l]++;
/******************************************************************************/
				p++;
			} while(/*err_count[l] <= 1000 && */p < num_packet);

			// Statistics and records
			//printf("Error Rate:\n");
			//err_rate[l] = err_count[l] / (double)((K-m)*p);
      		//bler_rate[l] = bler_count[l] / (double)p;
			//printf("%e, %e\n", err_rate[l], bler_rate[l]);
		}
		
		// Statistics and records
		printf("Error Rate:\n");
		err_sum = bler_sum = 0;
		for(l=0; l<L; l++)
		{
			err_sum += err_count[l];
			bler_sum += bler_count[l];
		}
		
		err_rate = err_sum / (double)((K-m)*p*L);
   		bler_rate = bler_sum / (double)(p*L);

		printf("%e, %e\n", err_rate, bler_rate);

		fprintf(ber, "%.1f ", Eb_No);
		fprintf(records, "%.1f ", Eb_No);
		fprintf(ber, "%e %e ", err_rate, bler_rate);
		fprintf(records, "    %e    %e", err_rate, bler_rate);
	/*
		for(l=0; l<L; l++)
		{
			fprintf(ber, "%e %e ", err_rate[l], bler_rate[l]);
			fprintf(records, "    %e    %e", err_rate[l], bler_rate[l]);
		}
*/
		fprintf(ber, "\n");
		fprintf(records, "\n");
		fflush(records);
		fflush(ber);
	}

	delete data_bit;
	delete coded_bit;
	delete Hk;
	delete Yk;
	delete Sk;
	delete Dk;
	delete LLR;
	delete S_random;
	delete de_inter;
	delete interleaved;
	delete QPSK_sym;
	delete ISI_I;
	delete ISI_Q;
	delete TX_signal_I;
	delete TX_signal_Q;
	delete RX_signal_I;
	delete RX_signal_Q;
	delete H_f;
	delete ts;
	delete delay;
	for(i=0; i<path_num; i++)
	{
		delete path_I[i];
		delete path_Q[i];
	}
	delete path_I;
	delete path_Q;
	for(i=0; i<path_num; i++)
		delete ISI_I[i];
	delete ISI_I;
	for(i=0; i<path_num; i++)
		delete ISI_Q[i];
	delete ISI_Q;
	delete gain;
	for(i=0; i<path_num; i++)
		delete fade[i];
	delete fade;
	delete err_count;
	delete bler_count;
	//delete err_rate;
	//delete bler_rate;
	
	end = time(NULL);
	printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
	fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
	fclose(ber);
	fclose(records);
	printf("This program is ended. Press 'Enter' to continue.\n");
	getchar();

	return 0;
}

void Trellis_diagram(int n, int k, int m, int *gp)
{
/**********************************************************************************/
/* Generate Trellis Diagram for (n=2, k=1, m=6) {171,133} CC code */
/**********************************************************************************/
	int i, j, input, out_bit, out_sym, from_state, to_state, tmp_state;
    int num_state = (int)pow(2.0,(double)m), num_in_sym = (int)pow(2.0,(double)k);

   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input of trellis diagram for (n, k, m)
      {
      	tmp_state = from_state;
         out_sym = 0;                      // output codeword symbol of trellis diagram
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         for(i=0; i<n; i++)
         {
         	out_bit = 0;						// output bit of trellis diagram
            for(j=m; j>=0; j--)
            	out_bit ^= ((tmp_state & gp[i]) >> j) & 1;  	// Calculate output bit

            out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
         }
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram

         Trellis[from_state][input].from = from_state;
         Trellis[from_state][input].to = to_state;
         Trellis[from_state][input].in = input;
         Trellis[from_state][input].out[0] = ((out_sym>>1)&1);
         Trellis[from_state][input].out[1] = out_sym&1;
      }
   }
}

void CC_Encoder(int *Data_in, int *Code_bit, int n, int m, int *gp, int Packet_length, int Num_state)
{
/******************************************************************************/
/* Convolutional Encoder (n=2, k=1, m=6) Generator polynomial: {171, 133}     */
/******************************************************************************/
   int i, j, l, from_state, tmp_state, out_bit;

   from_state = 0;
   for(i=0; i<Packet_length; i++)
   {
   	tmp_state = from_state;
//    out_sym = 0;                      // output codeword symbol of trellis diagram
		tmp_state = (tmp_state << 1) ^ (Data_in[i] & 0x01);  // read input bit
      for(j=0; j<n; j++)
      {
      	out_bit = 0;						// output bit of trellis diagram
         for(l=m; l>=0; l--)
         	out_bit ^= ((tmp_state & gp[j]) >> l) & 1;  	// Calculate output bit

         Code_bit[2*i+j] = out_bit;
//       out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
		}
      from_state = tmp_state & (Num_state-1); 					// to_state of trellis diagram
	}
}
/******************************************************************************/

void Viterbi_dec_SOVA(double *Data_in, int *Data_out, double *Soft_out, int Packet_length, int Num_state, int Num_in)
{
/*========================================================================================*/
/* Soft Output Viterbi Algorithm decoder (SOVA)                                           */
/* Convolutional Decoder (n=2, k=1, m=6) Generator polynomial: {171, 133}                 */
/*========================================================================================*/
	int i, j, l, q, pre_state;
   double **mju_f, *mju, *survival_metric, metric, ***branch_metric, **mju_b, mju_tmp;

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[K];	// input bit stream, K: packet length, global
            int code_bit[N];	// output code bit stream, N: packet length, global
            int state[K];		// state transition sequence, K: packet length, global
         };
   //struct surv survival_path[N_state], survival_temp[N_state];
   									// Survival_path[N_state], N_state: number of states of CC code, global

   struct surv *survival_path, *survival_temp;
   									// Survival_path[N_state], N_state: number of states of CC code, global
   survival_path = new struct surv [N_state];
   survival_temp = new struct surv [N_state];

   survival_metric = new double[Num_state];
   mju = new double[Num_in];							// minimum path metric
   mju_f = new double*[Packet_length+1];			// forward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_f[i] = new double[Num_state];
   mju_b = new double*[Packet_length+1];       	// backward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_b[i] = new double[Num_state];
   branch_metric = new double**[Packet_length];	// branch[time index][state][input]
   for(i=0; i<Packet_length; i++)
      branch_metric[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         branch_metric[i][j] = new double[Num_in];

   // Initialize survival path
   for(i=0; i<Num_state; i++)
   {
// 	survival_path[i].metric = DBL_MAX;	// Initial maximum value for Euclidean distance
		survival_path[i].metric = -DBL_MAX;	// Initial minimum value for cross-correlation
//    mju_f[0][i] = DBL_MAX;					// Initial maximum value for Euclidean distance
//    mju_b[K][i] = DBL_MAX;					// Initial maximum value for Euclidean distance
		mju_f[0][i] = -DBL_MAX;					// Initial minimum value for cross-correlation
      mju_b[K][i] = -DBL_MAX;					// Initial minimum value for cross-correlation
	}
   survival_path[0].metric = 0.0;
   mju_f[0][0] = 0.0;
   mju_b[K][0] = 0.0;

/*********************/
/* Forward Recursion */
/*********************/
	for(i=0; i<Packet_length; i++)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = -DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	// branch metric, Euclidean Distance
/*          branch_metric[i][j][l] = 0.0;
				branch_metric[i][j][l] += pow(Data_in[2*i]-(2*Trellis[j][l].out[0]-1),2);		// brahch metric
            branch_metric[i][j][l] += pow(Data_in[2*i+1]-(2*Trellis[j][l].out[1]-1),2);	// branch metric
            metric = survival_path[j].metric + branch_metric[i][j][l];
*/
				branch_metric[i][j][l] = 0.0;		// branch metric, Cross-correlation
            branch_metric[i][j][l] += (Data_in[2*i] * (2*Trellis[j][l].out[0]-1));		// brahch metric
            branch_metric[i][j][l] += (Data_in[2*i+1] * (2*Trellis[j][l].out[1]-1));	// branch metric
            metric = survival_path[j].metric + branch_metric[i][j][l];

            // find the survival path metric
//          if(metric < survival_metric[Trellis[j][l].to])	//	Euclidean distance (Minimize)
				if(metric > survival_metric[Trellis[j][l].to])	// Cross-correlation (Maximize)
            {
            	survival_metric[Trellis[j][l].to] = metric;

               // Record and refresh the survival path
               /*for(q=0; q<i; q++)
               {
               	survival_temp[Trellis[j][l].to].data_in[q] = survival_path[j].data_in[q];
                  survival_temp[Trellis[j][l].to].code_bit[2*q] = survival_path[j].code_bit[2*q];
                  survival_temp[Trellis[j][l].to].code_bit[2*q+1] = survival_path[j].code_bit[2*q+1];
                  survival_temp[Trellis[j][l].to].state[q] = survival_path[j].state[q];
               } */
               survival_temp[Trellis[j][l].to].data_in[i] = l;
               survival_temp[Trellis[j][l].to].code_bit[2*i] = Trellis[j][l].out[0];
               survival_temp[Trellis[j][l].to].code_bit[2*i+1] = Trellis[j][l].out[1];
               survival_temp[Trellis[j][l].to].state[i] = Trellis[j][l].from;
            }
			}

		// Record and refresh the survival path
      for(j=0; j<Num_state; j++)		// to_state index
      {
      	survival_path[j].metric = survival_metric[j];
         mju_f[i+1][j] = survival_metric[j];
         /*for(q=0; q<=i; q++)
         {
         	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
            survival_path[j].code_bit[2*q] = survival_temp[j].code_bit[2*q];
            survival_path[j].code_bit[2*q+1] = survival_temp[j].code_bit[2*q+1];
            survival_path[j].state[q] = survival_temp[j].state[q];
         }*/
      }
	}

   for(j=0; j<Num_state; j++)		// to_state index   
   {
      survival_path[j].data_in[Packet_length-1] = survival_temp[j].data_in[(Packet_length-1)];
      survival_path[j].code_bit[2*(Packet_length-1)] = survival_temp[j].code_bit[2*(Packet_length-1)];
      survival_path[j].code_bit[2*(Packet_length-1)+1] = survival_temp[j].code_bit[2*(Packet_length-1)+1];
  //    survival_path[j].state[Packet_length-1] =
	//	          survival_temp[j].state[Packet_length-1];
   }

   for(j=0; j<Num_state; j++)		// to_state index
   {
   	pre_state = survival_temp[j].state[Packet_length-1];  // from state

   	for( q=Packet_length-2; q>=0; q--)
   	{
      	survival_path[j].data_in[q] = survival_temp[pre_state].data_in[q];
	      survival_path[j].code_bit[2*q] = survival_temp[pre_state].code_bit[2*q];
   	   survival_path[j].code_bit[2*q+1] = survival_temp[pre_state].code_bit[2*q+1];
      	// survival_path[j].state[q] = survival_temp[pre_state].state[q];
	      pre_state = survival_temp[pre_state].state[q];  // from state
   	}
  	}

/****************************************/
/* Backward Recursion and Soft Decision */
/****************************************/
	for(i=Packet_length-1; i>=0; i--)
   {
   	for(j=0; j<Num_state; j++)				// Initialize the survival path metric
//    	survival_metric[j] = DBL_MAX;		// Initial maximum value for Euclidean distance
			survival_metric[j] = -DBL_MAX;	// Initial minimum value for cross-correlation

		for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)		// input bit
         {
         	metric = mju_b[i+1][Trellis[j][l].to] + branch_metric[i][j][l];

            // find the survival path metric
//          if(metric < survival_metric[j])	//	Euclidean distance (Minimize)
				if(metric > survival_metric[j])	// Cross-correlation (Maximize)
            	survival_metric[j] = metric;
			}

		// Record the survival path metric
      for(j=0; j<Num_state; j++)		// from_state index
      	mju_b[i][j] = survival_metric[j];
/*
      // LLR Calculation for the information bit
      mju[survival_path[0].data_in[i]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].data_in[i]+1)%2] = DBL_MAX;		//	Euclidean distance (Minimize)
		mju[(survival_path[0].data_in[i]+1)%2] = -DBL_MAX;    // Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      {
      	mju_tmp = mju_f[i][j] + branch_metric[i][j][(survival_path[0].data_in[i]+1)%2]
         			 + mju_b[i+1][Trellis[j][(survival_path[0].data_in[i]+1)%2].to];

//       if(mju_tmp < mju[(survival_path[0].data_in[i]+1)%2])	//	Euclidean distance (Minimize)
			if(mju_tmp > mju[(survival_path[0].data_in[i]+1)%2])	// Cross-correlation (Maximize)
         	mju[(survival_path[0].data_in[i]+1)%2] = mju_tmp;
		}

//    Soft_out[i] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[i] = mju[1] - mju[0];		// Cross-correlation
*/
      // LLR Calculation (for the 1st code bit)
      mju[survival_path[0].code_bit[2*i]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].code_bit[2*i]+1)%2] = DBL_MAX;			//	Euclidean distance (Minimize)
		mju[(survival_path[0].code_bit[2*i]+1)%2] = -DBL_MAX;    	// Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)
         	if(Trellis[j][l].out[0] != survival_path[0].code_bit[2*i])
            {
            	mju_tmp = mju_f[i][j] + branch_metric[i][j][l] + mju_b[i+1][Trellis[j][l].to];
//             if(mju_tmp < mju[(survival_path[0].code_bit[2*i]+1)%2])	//	Euclidean distance (Minimize)
					if(mju_tmp > mju[(survival_path[0].code_bit[2*i]+1)%2])	// Cross-correlation (Maximize)
               	mju[(survival_path[0].code_bit[2*i]+1)%2] = mju_tmp;
				}

//    Soft_out[2*i] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[2*i] = mju[1] - mju[0];		// Cross-correlation

		// LLR Calculation (for the 2nd code bit)
      mju[survival_path[0].code_bit[2*i+1]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].code_bit[2*i+1]+1)%2] = DBL_MAX;		//	Euclidean distance (Minimize)
		mju[(survival_path[0].code_bit[2*i+1]+1)%2] = -DBL_MAX;		// Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)
         	if(Trellis[j][l].out[1] != survival_path[0].code_bit[2*i+1])
            {
            	mju_tmp = mju_f[i][j] + branch_metric[i][j][l] + mju_b[i+1][Trellis[j][l].to];
//   	         if(mju_tmp < mju[(survival_path[0].code_bit[2*i+1]+1)%2])	//	Euclidean distance (Minimize)
					if(mju_tmp > mju[(survival_path[0].code_bit[2*i+1]+1)%2])	// Cross-correlation (Maximize)
               	mju[(survival_path[0].code_bit[2*i+1]+1)%2] = mju_tmp;
				}

//    Soft_out[2*i+1] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[2*i+1] = mju[1] - mju[0];		// Cross-correlation

      Data_out[i] = survival_path[0].data_in[i];
   }

   delete survival_path;
   delete survival_temp;

	delete[] survival_metric;
   delete[] mju;
   for(i=0; i<=Packet_length; i++)
      delete[] mju_f[i];
   delete[] mju_f;
   for(i=0; i<=Packet_length; i++)
       delete[] mju_b[i];
   delete[] mju_b;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete[] branch_metric[i][j];
   for(i=0; i<Packet_length; i++)
      delete[] branch_metric[i];
   delete[] branch_metric;
}
/******************************************************************************/

void AWGN_noise(float mu, double variance, double *noise)
{
//	const  float Pi = 3.14159265358979;
   double u1, u2;
   do
   {
   	u1 = (double)rand()/(double)RAND_MAX;
      u2 = (double)rand()/(double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrt(-2.0*log(u1))*cos(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
   *(noise+1) = (sqrt(-2.0*log(u1))*sin(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
}

void JakesFading(double f_c/*Hz*/, double v/*m/s*/, double t/*s*/, int type, double *fade)
{
	//const double C = 3.0e8;     // (m/s)
   //const float Pi = 3.14159265358979;
   int n, Np, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2, theta_n;

   lamda = C/f_c;     // wave length (meter)
   w_m = 2.0*Pi*v/lamda;    // maximum Doppler frequency
   Np = 2*(2*N_o+1);

   switch(type)
   {
   	case 1:
   		alpha = 0.0;
         T_c2 = (double)N_o;
         T_s2 = (double)N_o + 1.0;
         break;
      case 2:
      	alpha = 0.0;
         T_c2 = (double)N_o + 1.0;
         T_s2 = (double)N_o;
         break;
      case 3:
      	alpha = Pi/4.0;
         T_c2 = (double)N_o + 0.5;
         T_s2 = (double)N_o + 0.5;
         break;
      default:
      	printf("\nInvalid type selection for Jake's fading channel model.\n");
         break;
   }

	if(v < 0.0)
	{
		printf("Warning!! The vehicle speed is invalid.\n");
		cout << "Vehicle speed: " << v << "km/h" << endl;
		getchar();
   		//*(fade+0) = 1.0;
		//*(fade+1) = 0.0;
   }
   else
   {
   	*(fade+0) = sqrt(1.0/T_c2)*cos(alpha)*cos(w_m*t);
      *(fade+1) = sqrt(1.0/T_s2)*sin(alpha)*cos(w_m*t);

      for(n = 1; n <= N_o; n++)
      {
      	switch(type)
         {
         	case 1:
            	beta_n = (double)n*Pi/((double)N_o+1.0);
               break;
            case 2:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
            case 3:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
         	default:
            	break;
         }
         w_n = w_m*cos(2.0*Pi*(double)n/(double)Np);
            theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			//theta_n = 0.0;
         *(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n);
         *(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n);
		}
	}
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

void dfour1(double data[],unsigned long nn,int isign)
//double data[];
//unsigned long nn;
//int isign;
{
	unsigned long n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2)
   {
   	if (j > i)
      {
      	SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m)
      {
      	j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax)
   {
   	istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2)
      {
      	for (i=m;i<=n;i+=istep)
         {
         	j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
   }
}
#undef SWAP

void Multipath_fade_pattern(double **fade, int path_num)
{
	int i;
	double gain[2];

	for(i=0; i<path_num; i++)
	{
		JakesFading(fc, vc*1000/3600.0, ts[i], 2, &gain[0]);
		ts[i] += (OFDM_sym_duration * num_packet);
		fade[i][0] = gain[0];
		fade[i][1] = gain[1];
	}
}
