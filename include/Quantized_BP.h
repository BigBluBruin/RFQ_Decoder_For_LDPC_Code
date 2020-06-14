#pragma  once
#include "uniform_quantization.h"
#include "Parity_Check_Matrix_Info.h"
#include "channel_quantizer.h"
#include "full_BP.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <string.h>

class Quantized_BP{

public:
std::string h_filename;
std::string channel_recon_file, recon_file, threshold_file;
std::vector<double> parameters; //noise
Parity_Check_Matrix_Info h_ins;
channel_quantizer cq_ins;
int target_error;
int cardi_msg;
int max_iter;
std::vector<std::vector<double>> check_recons;
std::vector<std::vector<double>> vari_recons;
std::vector<double> channel_recons;
std::vector<std::vector<double>> check_threshold;
std::vector<std::vector<double>> vari_threshold;
public:
std::vector<int> total_frames;
std::vector<int> total_iterations;

public:


Quantized_BP(std::string H_filename, std::string Quantizer_name, std::vector<double> Parameters,
             int Target_error, std::string Channel_recon_file, std::string Recon_file,
             std::string Threshold_file);
void noise_generator(std::vector<double> & cwds,double parameter);
void fill_in();

//----
bool decoder(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits);

//-----Overloads of quantized-min sum decoder-----
bool decoder_min(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit);
bool decoder_min(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int l, int r);
//-----main simulation results-----
void main_simulation(const char ind[], const char suffix[],std::string filename, int l, int r);
void main_simulation(const char ind[], const char suffix[], int l, int r);
void main_simulation(const char ind[], const char suffix[],std::string filename);
void main_simulation(const char ind[],const char suffix[]);


//-----main simulation vertical layer decoder-----
void main_simulation_vertical_layered(const char ind[], const char suffix[],std::string filename, int l, int r, int layer_size);
bool decoder_min_vertical_layered(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int layer_size);
bool decoder_min_vertical_layered(std::vector<double> cwds ,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int l, int r, int layer_size);

//-----main simulation Horizontal  layer decoder-----
void main_simulation_horizontal_layered(const char ind[], const char suffix[],std::string filename, int l, int r, int layer_size);
bool decoder_min_horizontal_layered(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int layer_size);
bool decoder_min_horizontal_layered(std::vector<double> cwds ,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int l, int r, int layer_size);

//-----some tools-----
unsigned threshold_quantization(std::vector<double> & threshold, double value);
bool iscwds(std::vector<int> final_bits);
void generate_whole_bits(std::vector<int> & inputbit);
void generate_codeword(std::vector<int> &inputbit, std::vector<double> & codewords);

//-----what happened to wrong code ? -----
void decoder_min_track(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int ind);
void decoder_min_track(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int l, int r, int ind);
};