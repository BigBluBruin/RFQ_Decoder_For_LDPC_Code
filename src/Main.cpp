
#include "Quantized_BP.h"
#include "string.h"


int main(int argc,char* argv[])
{


    //-----Input Commands Instruction-----------------------------------------------------------------------------------------------------------
    //  -Hfile      parity check file name 
    //  -Gfile      G file name,  set "all_zeros" means all zero codewords are simulated
    //  -suf        suffix of RFQ related file
    //  -ebno       eb/no you want to calculated, for parallel compuation, only one ebno is allowed
    //  -tarerr     target number error
    //  -ind        suffix of output txtfile, for parallel computation
    //  -layersize  number of variable nodes in each group------- set "0" means flooding schedule
    //  -l          the number of bits allocated to integer part 
    //  -r          the number of bits allocated to fraction part
    //  -t          type: v-vertical h-horizontal, it only works when layered size is larger than 0
    // Example: ./runner -Hfile 80211_irr_648_1296.txt -Gfile 80211_G.txt  -suf 092 -ebno 1.40 -tarerr 10 -ind 0 -layersize 54 -l 5  -r 5 -t H
    //----------------------------------------------------------------------------------------------------------------------------------------- 


    int i ;
    int num_of_in=(argc-1)/2;
    std::string g_filename;
    std::string ind;
    int layer_size;
    int l,r,target_error;
    std::vector<double> parameters;
    std::string suffix;
    std::string parity_file;
    std::string layer_type;
    if (num_of_in != 10)
    {
        std::cout<<"Info: Number of Input is not equal to 10, please check again."<<std::endl;

    }
    for (i = 0; i < num_of_in; ++i)
    {
        if (strcmp(argv[2 * i + 1], "-Hfile") == 0)
        {
            parity_file = argv[2 * i + 2];
        }
        else if (strcmp(argv[2 * i + 1], "-Gfile") == 0)
        {
            g_filename = argv[2 * i + 2];
        }
        else if (strcmp(argv[2 * i + 1], "-suf") == 0)
        {
            suffix = argv[2 * i + 2];
        }
        else if (strcmp(argv[2 * i + 1], "-ebno") == 0)
        {
            parameters = {std::stod(argv[2 * i + 2])} ;
        }
        else if (strcmp(argv[2 * i + 1], "-tarerr") == 0)
        {
            target_error = std::stoi(argv[2 * i + 2]) ;
        }
        else if (strcmp(argv[2 * i + 1], "-ind") == 0)
        {
            ind = argv[2 * i + 2] ;
        }
        else if (strcmp(argv[2 * i + 1], "-layersize") == 0)
        {
            layer_size = std::stoi(argv[2 * i + 2]) ;
        }
        else if (strcmp(argv[2 * i + 1], "-l") == 0)
        {
            l = std::stoi(argv[2 * i + 2]) ;
        }
        else if (strcmp(argv[2 * i + 1], "-r") == 0)
        {
            r = std::stoi(argv[2 * i + 2]) ;
        }
        else if (strcmp(argv[2 * i + 1], "-t") == 0)
        {
            layer_type = argv[2 * i + 2] ;
        }
        else
        {
            std::cout<<argv[2 * i + 1]<<" is not recgonized ... plz check again"<<std::endl;
            return 0;
        }
    }


    std::string quantizer_file = "channel_quantizer_" + suffix + ".txt";
    std::string channel_recons_file = "channel_reconstruction_" + suffix + ".txt";
    std::string recons_file = "reconstruction_" + suffix + ".txt";
    std::string threshold_file = "threshold_" + suffix + ".txt";
    Quantized_BP this_quanbp(parity_file, quantizer_file, parameters, target_error, channel_recons_file, recons_file, threshold_file);
    this_quanbp.fill_in();

    if (layer_size == 0)
    {
        if (strcmp(g_filename.data(), "all_zeros") == 0)
        {
            // all-zero codeword simulation
            if (l == 0 || r == 0)
            {
                // no internal quantization
                this_quanbp.main_simulation(ind.data(), suffix.data());
            }
            else
            {
                // with internal quantization
                this_quanbp.main_simulation(ind.data(), suffix.data(), l, r);
            }
        }
        else
        {
            if (l == 0 || r == 0)
            {
                // no internal quantization{}
                this_quanbp.main_simulation(ind.data(), suffix.data(), g_filename);
            }
            else
            {
                // with internal quantization
                this_quanbp.main_simulation(ind.data(), suffix.data(), g_filename.data(), l, r);
            }
        }
    }
    else
    {
        if (strcmp(layer_type.data(),"v")==0||strcmp(layer_type.data(),"V")==0)
        {
            this_quanbp.main_simulation_vertical_layered(ind.data(),suffix.data(),g_filename.data(),l,r,layer_size);
        }
        else if (strcmp(layer_type.data(),"h")==0||strcmp(layer_type.data(),"H")==0)
        {
            this_quanbp.main_simulation_horizontal_layered(ind.data(),suffix.data(),g_filename.data(),l,r,layer_size);
        }
        else
        {
            std::cout<<"Info: there is no type: ||"<<layer_type<<"||, please check again ..."<<std::endl;
            return 0;
        }
        
        
    }
}

//-------------TRASH----------------------------
//-----This part is used for tracking trapping set using our decoding scheme---------------------------
// int total_num=35;
// int iteration =0;
// for (int cur_ind = 0; cur_ind<total_num; cur_ind++)
// {
//     //cur_ind=11;
//     std::string cur_initial_bit_filename="initial_bits_"+std::to_string(cur_ind)+".txt";
//     std::string cur_cwd_filename = "wrong_cwd_"+std::to_string(cur_ind)+".txt";
//     std::ifstream filehandle;
//     std::vector<int> inital_bits(this_quanbp.h_ins.vari_num,0);
//     std::vector<int> final_bits(this_quanbp.h_ins.vari_num,0);
//     std::vector<double> cwds(this_quanbp.h_ins.vari_num,0);
//     filehandle.open(cur_initial_bit_filename);
//     if(filehandle.is_open())
//     {
//         for(int ii=0; ii<this_quanbp.h_ins.vari_num;ii++)
//         {
//             filehandle>>inital_bits[ii];
//         }
//     }
//     else
//     {
//         std::cout<<"Info: Could not open INITIAL BIT file: "<<cur_initial_bit_filename<<". please check your input"<<std::endl;
//         return 0;
//     }
//     filehandle.close();
//     filehandle.open(cur_cwd_filename);
//     if (filehandle.is_open())
//     {
//         for (int ii = 0; ii < this_quanbp.h_ins.vari_num; ii++)
//         {
//             filehandle >> cwds[ii];
//         }
//     }
//     else
//     {
//         std::cout << "Info: Could not open WRONG CWDs file: " << cur_cwd_filename << ". please check your input" << std::endl;
//         return 0;
//     }
//     //this_quanbp.decoder_min_track(cwds,iteration,final_bits,inital_bits,cur_ind);
//     //this_quanbp.decoder_min_track(cwds,iteration,final_bits,inital_bits,5,5,cur_ind);
//     this_quanbp.decoder_min_track(cwds,iteration,final_bits,inital_bits,5,3,cur_ind);
// }