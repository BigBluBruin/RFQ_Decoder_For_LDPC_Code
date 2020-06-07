#include "channel_quantizer.h"

channel_quantizer::channel_quantizer(std::string Quantizer_file)
{
	quantizer_file = Quantizer_file;
}

void channel_quantizer::read_channel_quantizer()
{
	std::ifstream myfile(quantizer_file);
	if (myfile.is_open())
	{
        std::cout<<"------------------"<<std::endl;
        std::cout<<quantizer_file<<std::endl;
		myfile>>Min_value>>Max_value>>Partition_num;
		quantizer = D1i(Partition_num);
		for (int ii = 0; ii < Partition_num; ii++)
		{
			myfile >> quantizer[ii];
		}
		std::cout << "Info: finish message loading ..." << std::endl;
	}
	else
	{
		std::cout << "Wrong Info: Can't find quantization file..." << std::endl;
	}
}