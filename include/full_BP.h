#include <vector>
#include <numeric>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iostream>

//--------------------------------------------------
double sgn(double  input);
double vari_node_operation(std::vector<double> input);
double vari_node_operation(std::vector<double> input, int l, int r); 
//---------------------------------------------------
double boxplus(double argument_1, double argument_2);
double check_min_double(double & input_1, double & input_2);
double minstar(double & input_1, double & input_2);
unsigned min_sum(unsigned & input1, unsigned & input2, const unsigned & quantization_size);
double min_sum(double & input1, double & input2);
//---------------------------------------------------
double check_node_operation(std::vector<double> input);
double check_node_operation(std::vector<double> input, int l, int r);
double check_node_operation_min_double(std::vector<double> & input);
double check_node_operation_minstar(std::vector<double> & input);
double check_node_operation_att_min_double(std::vector<double> & input, const double att);
unsigned check_node_operation_min(std::vector<unsigned> & input, const unsigned & quan_size);
std::vector <double> check_node_operation_fast_minsum(std::vector<double> & input);
//--------------------------------------------------
void LP(double &in, int l, int r);
void quan_two_vec(std::vector<std::vector<double>>& aa, int l, int r);
void quan_one_vec(std::vector<double> &aa, int l, int r);
//--------------------------------------------------
void display(std::vector<int> input);
void display(std::vector<std::vector<int>> input);
//--------------------------------------------------
std::vector <double> check_sepcail_quan (int quan_size);
std::vector <double> vari_sepcail_recons (int quan_size);