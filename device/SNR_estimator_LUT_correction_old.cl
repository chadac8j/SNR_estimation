
#include "cordic.h"

//#define SNR_SYMBOL_LENGTH  1040  //Use more samples than num_samp_to_average
#define SNR_SYMBOL_LENGTH  1024  //Use more samples than num_samp_to_average
//#define SNR_SYMBOL_LENGTH  90  //Use QPSK frame preamble to estimate SNR


__kernel 
//__attribute__((task))
void snr_est_LUT_correction(	unsigned int slotLen  // do I need SOF indicator here or in output kernel?
			) 
{
	
	char2 shift_reg[SNR_SYMBOL_LENGTH+1];
	#pragma unroll
	for(uint i=0; i<=SNR_SYMBOL_LENGTH; i++){
		shift_reg[i].x = 0;
		shift_reg[i].y = 0;
	}
	unsigned long numerator=0;
	unsigned long noiseVarSum=0;
	unsigned long noiseVarSum_final=0;
	int tmpI = 0;
	int tmpQ = 0;
	int cordic_abs = 0;
	unsigned long abs_energy_sum = 0;
	int temp_snr_est = 0;

	int num_samp_to_average = 1024;  // so division is shift by 10 bits
	int bits_to_shift = 10;
	int while_loop_cntr = 0;	
	int carry = 0;
	
	while(1){
		__freqDetIn freqDetIn;
		
		freqDetIn=read_channel_intel(SNR_DET_DIN_LUT);  // need to create this new Channel

		if (freqDetIn.sof==1) {
			// reset the accumulators of numerator/denominator and loop counter at sof
			// Have a running sum of two main components of SNR estimate
			while_loop_cntr = 0;		
		}
		while_loop_cntr += 1;

		abs_energy_sum = 0;
//		int samp_idx = 0;
		#pragma unroll
		for(uint j=SNR_SYMBOL_LENGTH; j>0; j--){
			shift_reg[j].x = shift_reg[j-1].x;
			shift_reg[j].y = shift_reg[j-1].y;
			/* if it's more efficient, you can add newest data point and subtract oldest one
			   instead of doing the sum each loop iteration  */
//			if (samp_idx < num_samp_to_average) {
				// ********************************
				//   Magnitude Calculation
				// We multiply cordic inputs by 2^8 to ensure numerical precision.  Must remember to divide by 2^8 before final SNR estimate.
				//  We also don't normalize average, but keep the sum to avoid the loss of precision due to division.  This also must be accounted
				//    for in final SNR calculation
				// ********************************	
				tmpI = (int)(shift_reg[SNR_SYMBOL_LENGTH - j].x)<<8;
				tmpQ = (int)(shift_reg[SNR_SYMBOL_LENGTH - j].y)<<8;
				cordic_abs = mag_cordic(tmpI, tmpQ);
				abs_energy_sum += cordic_abs;
//			}
//			samp_idx += 1;
		}
		shift_reg[0].x = freqDetIn.data.x;
		shift_reg[0].y = freqDetIn.data.y;

		noiseVarSum = 0;
		#pragma unroll
		for(uint j=0; j<while_loop_cntr; j++)
		{
			tmpI = (int)(shift_reg[j].x)<<8;
			tmpQ = (int)(shift_reg[j].y)<<8;
			cordic_abs = mag_cordic(tmpI, tmpQ);
			noiseVarSum += ((cordic_abs<<bits_to_shift) - abs_energy_sum)*((cordic_abs<<bits_to_shift) - abs_energy_sum); //at this point we are 8 + bits_to_shift bits total shifted
		}
		

		// dividing twice and rounding down adds downward bias, but not dividing each operand leads to saturation beyond 32 integer bits
		carry = (1&(noiseVarSum>>(15-1)));  //the division by 2^15 here is to account for the <<8 done prior to cordic and then that value is squared.  There is also a multiplication by 2 that brings it from 16 down to 15 bit shifts
		noiseVarSum_final=(noiseVarSum>>15) + carry;
		noiseVarSum_final = noiseVarSum_final<<1; //multiply by 2

		numerator = ((abs_energy_sum*abs_energy_sum)<<(bits_to_shift))>>(2*8);  // the final division by 2^16 here is to account for the <<8 done prior to cordic and then that value is squared. 

		/* SNR estimate will be calculated in main.cpp, but will look like this:
		dout_estimate = (double)numerator/((double)(noiseVarSum_final));
		*/

	/* Here is the simplest version of the Matlab SNR estimate equation from the simulation:
	Es_No_est = 0.5*(num_trials*sum_mag^2)/((sum((num_trials*abs(y_I_vec + 1i*y_Q_vec) - sum_mag).^2)))
	*/

		
/*		printf("estimator Kernel loop ind=%d. abs_energy_avg=%d, abs_energy_sum=%lu \n", while_idx, abs_energy_avg, abs_energy_sum);
		printf("noiseVarSum=%lu \n", noiseVarSum);
		printf("noiseVarSum_final=%lu \n", noiseVarSum_final);
		printf("numerator=%lu \n", numerator);
		printf(" dout_estimate =  %f \n", dout_estimate);
*/

		if (while_loop_cntr==num_samp_to_average) { // return SNR estimate after num_samp_to_average samples
	        printf("In SNR kernel SOF_detect==True, numerator= %lu \n", numerator);
	        printf("In SNR kernel SOF_detect==True, denominator=%lu \n", noiseVarSum_final);
			temp_snr_est = round(100*log10((float)(numerator)/((float)(noiseVarSum_final))));
			write_channel_intel(SNR_DOUT, temp_snr_est);
		}
		
	}  // end main while loop
}



