#include "full_BP.h"

double boxplus(double argument_1, double argument_2)
{
    /*double sum1, max1, max2, output1;
	sum1 = argument_1 + argument_2;
	if(sum1 >= 0)
		max1 = sum1;
	else
		max1 = 0;
	
	if(argument_1 >= argument_2)
		max2 = argument_1;
	else
		max2 = argument_2;
	
	output1 = max1 - max2 + log(1 + exp(-abs(sum1))) - log(1+exp(-abs(argument_1 - argument_2)));
	return output1;*/
    double result = 2 * atanh(tanh(argument_1 / 2) * tanh(argument_2 / 2));
    if (result==INFINITY)
    {
        result = 3.0*pow(10,7);
    }
    else if (result==-INFINITY)
    {
        result = - 3.0*pow(10,7);

    } 
    return result;
}

double check_node_operation(std::vector<double> input)
{
    double output, fir, sec;
    fir = input[0];
    for (unsigned index = 0; index < input.size() - 1; index++)
    {
        sec = input[index + 1];
        fir = boxplus(fir, sec);
    }
    output = fir;
    return output;
}
double check_node_operation(std::vector<double> input, int l, int r)
{
    double output, fir, sec;
    fir = input[0];
    for (unsigned index = 0; index < input.size() - 1; index++)
    {
        sec = input[index + 1];
        fir = boxplus(fir, sec);
        LP(fir,l,r);
    }
    output = fir;
    return output;
}

double vari_node_operation(std::vector<double> input)
{
    double sum;
    sum = std::accumulate(input.begin(), input.end(), 0.0);
    return sum;
}
double vari_node_operation(std::vector<double> input, int l, int r)
{
    double output, fir, sec;
    fir = input[0];
    for (unsigned index = 0; index < input.size() - 1; index++)
    {
        sec = input[index + 1];
        fir = fir+sec;
        LP(fir,l,r);
    }
    output = fir;
    return output;
}

unsigned min_sum(unsigned &input1, unsigned &input2, const unsigned &quantization_size)
{
    if ((input1 + 1) > quantization_size / 2 && (input2 + 1) > quantization_size / 2)
    {
        return std::min(input1, input2);
    }
    else if ((input1 + 1) <= quantization_size / 2 && (input2 + 1) <= quantization_size / 2)
    {
        return std::min(quantization_size - 1 - input1, quantization_size - 1 - input2);
    }
    else
    {
        if (input1 > input2)
        {
            return std::max(input2, quantization_size - 1 - input1);
        }
        else
        {
            return std::max(input1, quantization_size - 1 - input2);
        }
    }
}


double min_sum(double & input1, double & input2)
{
    return sgn(input1)*sgn(input2)*std::min(input1,input2);
}

unsigned check_node_operation_min(std::vector<unsigned> &input, const unsigned &quan_size)
{
    unsigned fir = input[0];
    unsigned sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = min_sum(fir, sec, quan_size);
    }
    return fir;
}


double sgn(double input)
{
    if (input >= 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}
double check_min_double(double &input_1, double &input_2)
{
    return sgn(input_1) * sgn(input_2) * std::min(abs(input_2), abs(input_1));
}

double check_node_operation_min_double(std::vector<double> &input)
{
    double fir = input[0];
    double sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = check_min_double(fir, sec);
    }
    return fir;
}

double check_node_operation_att_min_double(std::vector<double> &input, const double att)
{
    double fir = input[0];
    double sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = check_min_double(fir, sec);
    }
    return att * fir;
}



void LP(double &in, int l, int r)
{

    double rpow, lpow, minpow, maxpow;
    rpow = pow(2, r);
    lpow = pow(2, (l - 1));
    minpow = (-1) * lpow;
    maxpow = lpow - (double)1.0 / (double)rpow;
    if (in > maxpow)
    {
        (in) = maxpow;
    }
    else if (in < minpow)
    {
        (in) = minpow;
    }
    else
    {
        (in) = (double)floor((in)*rpow + 0.5) / (double)rpow;
    }
}

void quan_two_vec(std::vector<std::vector<double>> &aa, int l, int r)
{
    for (auto &inner : aa)
    {
        for (auto &bb : inner)
        {
            LP(bb, l, r);
        }
    }
}

void quan_one_vec(std::vector<double> &aa, int l, int r)
{
    for (auto &bb : aa)
    {
        LP(bb, l, r);
    }
}

double minstar(double &input_1, double &input_2)
{
    double part1= -sgn(input_1)*sgn(input_2)*std::min(abs(input_1), abs(input_2));
    double part2= log(1+exp(-abs(input_2-input_1)));
    double part3= log(1+exp(-abs(input_1+input_2)));
    return sgn(input_1)*sgn(input_2)*-1*(part1+part2-part3);
}

double check_node_operation_minstar(std::vector<double> & input)
{
    double fir = input[0];
    double sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = minstar(fir, sec);
    }
    return fir;
}

void display(std::vector<int> input)
{
    for (const auto aa:input)
    {
        std::cout<<aa<<"  ";
    }
    std::cout<<"  "<<std::endl;

}
void display(std::vector<std::vector<int>> input)
{
    for (const auto aa:input)
    {
        display(aa);
    }
}