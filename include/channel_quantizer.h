#pragma once
#include "allocate_pointer.h"
#include <fstream>
#include <string>
#include <iostream>

class channel_quantizer
{
public:
	std::string quantizer_file;
	int* quantizer;
	double Max_value;
	double Min_value;
	int Partition_num;

public:
	channel_quantizer(std::string Quantizer_file=" ");
	void read_channel_quantizer();
};