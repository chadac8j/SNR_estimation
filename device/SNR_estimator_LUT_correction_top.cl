/******************************************************************************
*  @file    snr_estimator_top.cl
*  @author  Chad Cole
*  @date    6/3/2021
*  @version 1.0
*
*  @brief Device side OpenCL implemenation of SNR Estimator.
*                                                             
*
*  @section DESCRIPTION
*
*  This is the top level device side interface for a FPGA implemetation 
*  of a SNR Estimator that will be used in the DVB-S2 and DVB-S2X waveforms. It uses channels to 
*  stitch together three OpenCL kernels.  One for reading data from memory, 
*  one for SNR estimation, and one for writing the SNR estimate
*  back to memory.  
*
*  Inputs:  char - I channel input
*           char - Q channel input
*           uint - data length
*
*  Outputs: Unsigned Long - Two values: The numerator and denominator of the SNR Estimate for current data frame.
*										These values can be divided and turned into a floating point number in main.cpp
*
*  Assumptions:  Data coming in is at a rate of one sample per symbol
*
*******************************************************************************/

#define SNR_DWELL_LENGTH  1024

#pragma OPENCL EXTENSION cl_intel_channels : enable


//data format for data sent to coarse frequency detector, reuse struct for SNR estimator
typedef struct {
	char2 data;         
	char sof;           
	unsigned short plHeaderCnt;
	int pilotCumCnt;
	char pilotsActive;
	char pilotCnt;
	char pls_lock;
	char sof_lock;
	char tracking_active;
	char pls_active;
}__freqDetIn;

// Channel declarations
channel __freqDetIn SNR_DET_DIN_LUT           __attribute__((depth(8)));
channel int SNR_DOUT           __attribute__((depth(8)));

__kernel 
void data_in(     __global char* dataIn_I, 
                  __global char* dataIn_Q,
                  unsigned int dataInLen) 
{
	printf("In read kernel, dataInLen= %d \n", dataInLen);
	for(uint i=0; i< dataInLen; i++){ 
		__freqDetIn din;
		din.data.x = dataIn_I[i];
		din.data.y = dataIn_Q[i];

		// for testing...
		//  <maybe_a_counter_for_this?>
		// set SOF on last of 90 SOF/PLS header symbols
		if ((i%SNR_DWELL_LENGTH) == 0){
			printf("In read kernel, sof==True, input ind= %d \n", i);
			din.sof = 1;
		}else{
			din.sof = 0;
		}

		write_channel_intel(SNR_DET_DIN_LUT, din);
	} //end for
}


__kernel 
void data_out(	__global short* snr_est_out,
				int num_output_frames) 
{
// give output kernel time to look for SOF indicator - wait, the read_channel_intel() is blocking so no wait necessary?
		int output_ind = 0;
		for(uint i=0; i< num_output_frames; i++){ 
			short temp_snr_est = 0;
			temp_snr_est = read_channel_intel(SNR_DOUT);
			snr_est_out[output_ind] = temp_snr_est;
	        printf("In write kernel, snr_est_out= %d, output_ind=%d \n", snr_est_out[output_ind], output_ind);
			output_ind += 1;
	        //printf("In write kernel, denominator=%lu \n", snr_est_out[1]);
		}
}

// Include the datapath kernels
#include "SNR_estimator_LUT_correction.cl"

//void data_out(	__global int* restrict snr_est_out,
