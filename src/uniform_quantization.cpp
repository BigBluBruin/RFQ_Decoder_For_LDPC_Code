#include "uniform_quantization.h"
int Quantize( double& in, double& min, double& max, double& interval, int& total)
{
	//This function realized quantization
	//Input low: lowest num; high: highest num
	//total maximu_num
	if (in > max)
	{
		return total - 1;
	}
	else if (in<min)
	{
		return  0;
	}
	else
	{
		return  (int)((in - min) / interval);
	}
}