/******************************************************************************
*  @file    cordic.
*  @author  Jason DuChez
*  @date    4/18/2019
*  @version 1.0
*
*  @brief This is a CORDIC implementation.     	*                                                    
*	@Rev 1 - Initial support is for 16bit inputs and outputs.
*     The cordic defaults to using 24 iterations and supports
*     arctan, and vector magnitude functionality.  Through the
*     following function calls:
*         unsigned short mag_cordic(short x, short y)
*         short arctan_cordic(short x, short y)
*     Additional CORDIC functions can easily be added.
*
*
*******************************************************************************/

#ifndef CORDIC_H_
#define CORDIC_H_

#include "round.h"


#define SHIFT 16
#define CORDIC_GAIN (0x9b75)  //  2^16*0.60725293 

typedef struct vector3 {
	int x;
	int y;
	int z;
} vector3;

typedef struct cos_sin {
	int cos;
	int sin;
} cos_sin;

typedef enum { ROTATION, VECTOR } cordic_mode_t;


//ATANTABLESZ sets the number of iterations the cordic
// will use for calcuation.  if this changes the 
// atanTable should be regenerated.  The following MATALB
// code can be used for re-generating the table
//
//	ATANTABLESZ = 24;
//	atanTable=zeros(1,ATANTABLESZ);
//	
//	for n=0:ATANTABLESZ-1
//	    atanTable(n+1)=(round((2^ATANTABLESZ)*atan(2^-n)/pi))
//	end
//
#define ATANTABLESZ 24
__constant int atanTable[] = {
	4194304,
	2476042,
	1308273,
	664100,
	333339,
	166832,
	83436,
	41721,
	20861,
	10430,
	5215,
	2608,
	1304,
	652,
	326,
	163,
	81,
	41,
	20,
	10,
	5,
	3,
	1,
	1
};



/*
* This is to reduce our amount of code for one more if statement per
* loop. It checks if (z < 0) and if (y >= 0). Depending on the mode, it
* will return
*/
int rot_decision(cordic_mode_t mode, int val)
{

	int result;
	if (val < 0)
	{
		result = 1;
	}
	else {
		result = 0;
	}

	/* The function will return the result if we are in rotation mode,
	* and the opposite if we are in vector mode. This is because rotation
	* mode runs if( z < 0) whereas vector mode runs if (y >= 0)
	*/
	if (mode == ROTATION)
	{
		return result;
	}
	else {
		return !result;
	}
}


//~ /**
//~ * Cordic in rotation mode. 
//~ * x: The x coordinate in the vector.
//~ * y: The y coordinate in the vector.
//~ * z: The angle accumulator. Reaches 0 after n iterations  
//~ *
//~ */
//~ vector3 cordic(long x, long y, long z, cordic_mode_t mode)
//~ {
	//~ int i;
	//~ long x_temp;
	//~ long * val;

	//~ // If the mode is rotation, our check for rotation direction is on the z value, but
	//~ // if it is in vector, it is on the y value.
	//~ if (mode == ROTATION)
	//~ {
		//~ val = &z;
	//~ }
	//~ else {
		//~ val = &y;
	//~ }
	
	//~ //x = x << SHIFT;
	//~ //y = y << SHIFT;
	//~ //z = z << SHIFT;

	//~ for (i = 0; i < ATANTABLESZ; i++)
	//~ {
		//~ x_temp = x;
		//~ //if (rot_decision(mode, *val))
		//~ if (rot_decision(mode, *val))
		//~ {
			//~ x = x + (y >> i);
			//~ y = y - (x_temp >> i);
			//~ z = z + atanTable[i];
		//~ }
		//~ else {
			//~ x = x - (y >> i);
			//~ y = y + (x_temp >> i);
			//~ z = z - atanTable[i];
		//~ }
	//~ }
	//~ vector3 result;
	//~ result.x = x;
	//~ result.y = y;
	//~ result.z = z;
	//~ return result;
//~ }



/**
* Pipelined Cordic in rotation mode. 
* x: The x coordinate in the vector.
* y: The y coordinate in the vector.
* z: The angle accumulator. Reaches 0 after n iterations  
*
*/
vector3 cordic(long x, long y, long z, cordic_mode_t mode)
{

	long * val;

	long __attribute__((register)) x_vec[ATANTABLESZ+1]; 
	long __attribute__((register)) y_vec[ATANTABLESZ+1];
	long __attribute__((register)) z_vec[ATANTABLESZ+1];

	x_vec[0]=x;
	y_vec[0]=y;
	z_vec[0]=z;	

	#pragma unroll
	for (int i = 0; i < ATANTABLESZ; i++)
	{
		// If the mode is rotation, our check for rotation direction is on the z value, but
		// if it is in vector, it is on the y value.
		if (mode == ROTATION)
		{
			val = &z_vec[i];
		}
		else {
			val = &y_vec[i];
		}

		if (rot_decision(mode, *val))
		{
			x_vec[i+1] = x_vec[i] + (y_vec[i] >> i);
			y_vec[i+1] = y_vec[i] - (x_vec[i] >> i);
			z_vec[i+1] = z_vec[i] + atanTable[i];
		}
		else {
			x_vec[i+1] = x_vec[i] - (y_vec[i] >> i);
			y_vec[i+1] = y_vec[i] + (x_vec[i] >> i);
			z_vec[i+1] = z_vec[i] - atanTable[i];
		}
	}
	vector3 result;
	result.x = x_vec[ATANTABLESZ];
	result.y = y_vec[ATANTABLESZ];
	result.z = z_vec[ATANTABLESZ];
	return result;
}


//detemines the quadrant and makes appropriate shifts
// this allows to cordic to opeate from 0 to 2pi
struct vector3 checkQuadrant(short x, short y, short z){
	int quadrant;
	if (x>0 & y>=0){
		quadrant=0;
	}else if (x<=0 & y>0){
		quadrant=1;
	}else if (x<0 & y<=0){
		quadrant=2;
	}else {
		quadrant=3;  
	}
	
	if (x>=0 & y>=0){
		x=x; y=y; 
	}else if (x<0 & y>=0){
		int xtmp=-x; x=y; y=xtmp; 
	}else if (x<0 & y<0){
		x=-x; y=-y; 
	}else {
		int ytmp=-y; y=x; x=ytmp;  
	}
	
	struct vector3 dout;
	dout.x=x;
	dout.y=y;
	dout.z=quadrant;
	
	return dout;
}

//this function returns the sin and cos for a given angle (theta).
// This is intended to be used for 24bit input values.
//theta = 2^24*(angle_in_radians)/(2*pi),  where angle_in_radians
// ranges from 0 to 2*pi.  This value must be no larger than 24 bits.
// 0x0 represents 0 and 0xffffff represents 2*pi
struct cos_sin sin_cos_cordic_24b(int theta){
	
	const int cordic_gain_inv_24b=0x9B6F23;  //(1/1.647 * 2^24)
	int theta_corrected;
	char invert_cos;
	char invert_sin;
	
	char bit22 = (theta>>22)&0x1;
	char bit23 = (theta>>23)&0x1;
	int theta_22b = theta & 0x003fffff;
		
	if ( (bit23==0) & (bit22==0)){
		theta_corrected=theta_22b;
		invert_cos=0;
		invert_sin=0;
	}
	if ( (bit23==0) & (bit22==1)){
		theta_corrected=0x003fffff-theta_22b;
		invert_cos=1;
		invert_sin=0;
	}
	if ( (bit23==1) & (bit22==0)){
		theta_corrected=theta_22b;
		invert_cos=1;
		invert_sin=1;
	}
	if ( (bit23==1) & (bit22==1)){
		theta_corrected=0x003fffff-theta_22b;
		invert_cos=0;
		invert_sin=1;
	}
	
	//printf("%d %d %x %x %x %d\n", bit23, bit22, theta, theta_22b, theta_corrected, invert_cos);
		
	struct vector3 output;	
	output = cordic( cordic_gain_inv_24b, 0, theta_corrected<<1, ROTATION);
	
	struct cos_sin cos_sin_s;
	
	if (invert_cos==1)
		cos_sin_s.cos=-output.x;
	else
		cos_sin_s.cos=output.x;
	
	if (invert_sin==1)	
		cos_sin_s.sin=-output.y;
	else
		cos_sin_s.sin=output.y;
	
	return cos_sin_s;
	
}

// int sin_cordic_24b(int theta){
	// return cordic( CORDIC_GAIN, 0, theta, ROTATION).y;
// }


//the arctan functtion function returns a value replresenting
// an angle from 0 to 2*pi.  The hexadecimal MATALB equivalent is
// dec2hex(round(2^16*wrapTo2Pi(angle(x+yi))/(2*pi)))
unsigned short arctan_cordic(short x, short y){
	
	unsigned short atan;
	
	int xi = x << SHIFT;
	int yi = y << SHIFT;
	
	//check for conditions when the sample is on an axis
	//this can lead to the angle being in the wrong
	//quadrant due to rounding issues
	if (x==0){
		if (y>=0) atan=16384;
		else      atan=49152;
	}else if (y==0){
		if (x>=0) atan=0;
		else      atan=32768;
	}else{
	
		struct vector3 in = checkQuadrant( x, y, 0);
		char quadrant=in.z;
		
		unsigned int atan_i = cordic( in.x, in.y, 0, VECTOR).z;
		unsigned short atan_1q = round_i(atan_i, 9);
		//printf("%d %d \n", quadrant,  ((atan_1q)&0x3fff) );
		atan = (quadrant<<14) | ((atan_1q)&0x3fff);
	}
	
	return atan;
}

//the arctan functtion function returns a value replresenting
// an angle from 0 to 2*pi.  The hexadecimal MATALB equivalent is
// dec2hex(round(2^24*wrapTo2Pi(angle(x+yi))/(2*pi)))
unsigned int arctan_cordic_24b(short x, short y){
	
	unsigned int atan;
	
	int xi = x << SHIFT;
	int yi = y << SHIFT;
	
	//check for conditions when the sample is on an axis
	//this can lead to the angle being in the wrong
	//quadrant due to rounding issues
	if (x==0){
		if (y>=0) atan=0x400000;
		else      atan=0xc00000;
	}else if (y==0){
		if (x>=0) atan=0;
		else      atan=0x800000;
	}else{
	
		struct vector3 in = checkQuadrant( x, y, 0);
		char quadrant=in.z;
		
		unsigned int atan_i = cordic( in.x, in.y, 0, VECTOR).z;
		unsigned int atan_1q = round_i(atan_i, 1);
		
		//atan = (quadrant<<22) | ((atan_1q)&0x3fffff);
		atan = ((quadrant<<22) + (atan_1q)) & 0xffffff;
	}
	
	return atan;
}


// in vector mode, the magnitude equals x[n]
// x[n]=An*sqrt(x[0]^2+y[0]^2)
// where An is the CORDIC_GAIN
// The hexadecimal MATALB equivalent is
// dec2hex(round(abs(x+yi)))
unsigned int mag_cordic(int x, int y){
	
	//works in the 1st quadrant
	if (x<0) x=-x;
	if (y<0) y=-y;
	
		
	unsigned int magScaled = cordic( x, y, 0, VECTOR).x;
	unsigned long long mag = (unsigned long)magScaled * (unsigned long)CORDIC_GAIN;
	long magScaled_l = round_l(mag, 16);
	unsigned int magScaled_rnd = (unsigned int) magScaled_l;
	
	return magScaled_rnd;
}


// in vector mode, the magnitude equals x[n]
// x[n]=An*sqrt(x[0]^2+y[0]^2)
// where An is the CORDIC_GAIN
// The hexadecimal MATALB equivalent is
// dec2hex(round(abs(x+yi)))
unsigned long mag_cordic_l(int x, int y){
	
	//works in the 1st quadrant
	if (x<0) x=-x;
	if (y<0) y=-y;
	
		
	unsigned int magScaled = cordic( x, y, 0, VECTOR).x;
	unsigned long mag = (unsigned long)magScaled * (unsigned long)CORDIC_GAIN;
		
	return mag;
}



#endif
