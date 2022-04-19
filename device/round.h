#ifndef ROUND_H_
#define ROUND_H_


short round_s(short din, char bits)
{
	short dout_us;
	if (((din>>(bits-1))&0x1)==0x1)
		dout_us=(din>>bits)+1;
	else
		dout_us=(din>>bits);
	
	return dout_us;
}

int round_i(int din, char bits)
{
	int dout_us;
	if (((din>>(bits-1))&0x1)==0x1)
		dout_us=(int)(din>>bits)+1;
	else
		dout_us=(int)(din>>bits);
	
	//signed dout=(signed)dout_us;
	return dout_us;
}

long round_l(long din, char bits)
{
	long dout_us;
	if (((din>>(bits-1))&0x1)==0x1)
		dout_us=(din>>bits)+1;
	else
		dout_us=(din>>bits);
	
	return dout_us;
}

#endif