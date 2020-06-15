#include "Quantized_BP.h"

Quantized_BP::Quantized_BP(std::string H_filename, std::string Quantizer_name, std::vector<double> Parameters, int Target_error, std::string Channel_recon_file, std::string Recon_file, std::string Threshold_file)
{
    h_filename = H_filename;
    parameters = Parameters;
    target_error = Target_error;
    total_frames.assign(parameters.size(), -1);
    total_iterations.assign(parameters.size(), -1);
    cq_ins.quantizer_file = Quantizer_name;
    channel_recon_file = Channel_recon_file;
    recon_file = Recon_file;
    threshold_file = Threshold_file;
}

void Quantized_BP::fill_in()
{
    //--------parity check matrix info---------------------
    h_ins.filename = h_filename;
    if (h_ins.Read_Parity_Check_Matrix())
    {
        std::cout << "Parity Matrix Success"
                  << std::endl;
    }
    else
    {
        std::cout << "Parity Check Matrix Fails"
                  << std::endl;
    }
    //------------------------------------------------------

    //---------quantization info----------------------------
    cq_ins.read_channel_quantizer();
    //-------------------------------------------------------

    //---------channel reconstruction------------------------
    std::ifstream handle_channel_rec(channel_recon_file);
    unsigned int quan_size;
    if (handle_channel_rec.is_open())
    {
        handle_channel_rec >> quan_size;
        channel_recons.assign(quan_size, -1.0);
        for (unsigned index = 0; index < quan_size; index++)
            handle_channel_rec >> channel_recons[index];
        std::cout << "Channel reconstruction file Successfull..."
                  << std::endl;
    }
    else
    {
        std::cout << "Channel reconstruction file NOT found"
                  << std::endl;
    }
    //-------------------------------------------------------

    //--------reconstruction file-----------------------------
    std::ifstream handel_reconstruction(recon_file);
    if (handel_reconstruction.is_open())
    {
        handel_reconstruction >> quan_size >> max_iter;
        check_recons.resize(max_iter);
        vari_recons.resize(max_iter);
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            check_recons[index1].assign(quan_size, -1);
            for (unsigned index2 = 0; index2 < quan_size; index2++)
            {
                handel_reconstruction >> check_recons[index1][index2];
            }
        }
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            vari_recons[index1].assign(quan_size, -1);
            for (unsigned index2 = 0; index2 < quan_size; index2++)
            {
                handel_reconstruction >> vari_recons[index1][index2];
            }
        }
        std::cout << "Reconstruction file Successfull..."
                  << std::endl;
    }
    else
    {
        std::cout << "Reconstruction file NOT found"
                  << std::endl;
    }

    //--------------------------------------------------------

    //--------threshold files---------------------------------
    std::ifstream handle_quan_threshold(threshold_file);
    if (handle_quan_threshold.is_open())
    {
        handle_quan_threshold >> quan_size >> max_iter;
        check_threshold.assign(max_iter, std::vector<double>(1, -1));
        vari_threshold.assign(max_iter, std::vector<double>(1, -1));
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            check_threshold[index1].assign(quan_size / 2 - 1, -1);
            for (unsigned index2 = 0; index2 < quan_size / 2 - 1; index2++)
            {
                handle_quan_threshold >> check_threshold[index1][index2];
            }
        }
        for (int index1 = 0; index1 < max_iter; index1++)
        {
            vari_threshold[index1].assign(quan_size / 2 - 1, -1);
            for (unsigned index2 = 0; index2 < quan_size / 2 - 1; index2++)
            {
                handle_quan_threshold >> vari_threshold[index1][index2];
            }
        }
        std::cout << "Threshold file Successfull..."
                  << std::endl;
    }
    else
    {
        std::cout << "Threshold file NOT found"
                  << std::endl;
    }
    //--------------------------------------------------------
}

unsigned Quantized_BP::threshold_quantization(std::vector<double> &threshold, double value)
{
    unsigned ind = 10000;
    unsigned quan_size = threshold.size() + 1;
    bool negtive = true;
    if (value > 0)
    {
        negtive = false;
        value = -1 * value;
    }
    if (value <= threshold[0])
    {
        ind = 0;
    }
    else if (value > threshold[quan_size - 2])
    {
        ind = quan_size - 1;
    }
    else
    {
        for (unsigned index = 0; index < quan_size - 2; index++)
        {
            if (value > threshold[index] && value <= threshold[index + 1])
            {
                ind = index + 1;
            }
        }
    }
    if (!negtive)
    {
        ind = 2 * quan_size - 1 - ind;
    }
    return ind;
}

void Quantized_BP::noise_generator(std::vector<double> &cwds, double parameter)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> distribution(0, sqrt(parameter));
    for (unsigned inc = 0; inc < cwds.size(); inc++)
    {
        cwds[inc] += distribution(generator);
    }
}

bool Quantized_BP::iscwds(std::vector<int> final_bits)
{
    std::vector<int> edge_bits(h_ins.edge_num, -1);
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            edge_bits[h_ins.edge_v[ii][jj]] = final_bits[ii];
        }
    }

    for (int ii = 0; ii < h_ins.check_num; ii++)
    {
        int checks = 0;
        for (int jj = 0; jj < h_ins.check_degreetable[ii]; jj++)
        {
            checks = (checks + edge_bits[h_ins.edge_c[ii][jj]]) % 2;
        }
        if (checks > 0)
        {
            return false;
        }
    }

    return true;
}

bool Quantized_BP::decoder(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits)
{
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<double> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quantized_val;
    double full_return;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double final_codewords;
    std::vector<double> message;
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = rx[ii];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj)
                    {
                        message.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if (message.size() != cur_dc - 1)
                {
                    std::cout << "not collecto enough data from variable nodes" << std::endl;
                    std::cout << message.size() << " " << cur_dc - 1 << std::endl;
                }
                full_return = check_node_operation(message);
                quantized_val = threshold_quantization(check_threshold[cur_iter], full_return);
                if (quantized_val != 10000)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]] = check_recons[cur_iter][quantized_val];
                    //std::cout<<msg_c2v[h_ins.edge_c[ii][jj]]<<std::endl;
                }
                else
                {
                    std::cout << "variable node quantization failed..." << std::endl;
                }
            }
        }

        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message.clear();
                message.push_back(rx[ii]);
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj)
                    {
                        message.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                if (message.size() != cur_dv)
                {
                    std::cout << "not collect enough data from chekc nodes" << std::endl;
                }
                full_return = vari_node_operation(message);
                quantized_val = threshold_quantization(vari_threshold[cur_iter], full_return);
                if (quantized_val != 10000)
                {
                    msg_v2c[h_ins.edge_v[ii][jj]] = vari_recons[cur_iter][quantized_val];
                }
                else
                {
                    std::cout << "variable node quantization failed..." << std::endl;
                }
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {

            int cur_dv = h_ins.vari_degreetable[ii];
            message.clear();
            message.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message);
            //std::cout<<final_codewords<<"  ";
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        //std::cout<<std::endl;
        //return false;

        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            int sum = 0;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                sum = sum + final_bits[ii];
            }
            if (sum == 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    iteration = iteration + max_iter;
    return false;
}

void Quantized_BP::main_simulation(const char ind[], const char suffix[])
{
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
    //-----Basic Information Output------------
    std::cout << "Info: All zeros Codewords..." << std::endl;
    std::cout << "Quantization: Full_Precision" << std::endl;
    //----------------------------------------------
    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii])) / (2.0 * h_ins.rate);
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            initial_bits.assign(h_ins.vari_num, 0);
            //add noise
            noise_generator(codewords, cur_para);
            //bool result = decoder(codewords, iteration, final_bits);
            bool result = decoder_min(codewords, iteration, final_bits, initial_bits);
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}

void Quantized_BP::main_simulation(const char ind[], const char suffix[], int l, int r)
{
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
    //-----Basic Information Output------------
    std::cout << "Info: All zeros Codewords..." << std::endl;
    std::cout << "Quantization: Integer--" << l << "--decimal--" << r << std::endl;
    quan_two_vec(vari_recons, l, r);
    quan_two_vec(vari_threshold, l, r);
    quan_two_vec(check_recons, l, r);
    quan_two_vec(check_threshold, l, r);
    quan_one_vec(channel_recons, l, r);
    //----------------------------------------------
    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii]) / (2.0 * h_ins.rate));
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            initial_bits.assign(h_ins.vari_num, 0);
            //add noise
            noise_generator(codewords, cur_para);
            bool result = decoder_min(codewords, iteration, final_bits, initial_bits, l, r);
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
                //you can also output failed codewords here....
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}

void Quantized_BP::main_simulation(const char ind[], const char suffix[], std::string filename)
{
    // random codeword generation, without any quantization
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<> d(1, 0.5);

    //-----Basic Information Output-----------------
    std::cout << "Info: Random Codewords..." << std::endl;
    std::cout << "Quantization: Full_Precision" << std::endl;
    //----------------------------------------------
    h_ins.Read_G_Matrix(filename);
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii])) / (2.0 * h_ins.rate);
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            initial_bits.assign(h_ins.vari_num, 0);
            for (unsigned ii = 0; ii < (h_ins.vari_num - h_ins.check_num); ii++)
            {
                initial_bits[ii] = d(gen);
            }
            generate_whole_bits(initial_bits);
            generate_codeword(initial_bits, codewords);
            //add noise
            noise_generator(codewords, cur_para);
            bool result = decoder_min(codewords, iteration, final_bits, initial_bits);
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
                //you can also output failed codewords here....
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}

void Quantized_BP::main_simulation(const char ind[], const char suffix[], std::string filename, int l, int r)
{
    // random codeword generation, with quantization
    std::random_device rd;
    std::mt19937 gen(rd());
    // perform 4 trials, each succeeds 1 in 2 times
    std::binomial_distribution<> d(1, 0.5);
    h_ins.Read_G_Matrix(filename);
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
    //-----Basic Information Output-----
    std::cout << "Info: Random Codewords..." << std::endl;
    std::cout << "Quantization: l-" << l << " bit, r-" << r << " bit..." << std::endl;
    //--------------------------
    quan_two_vec(vari_recons, l, r);
    quan_two_vec(vari_threshold, l, r);
    quan_two_vec(check_recons, l, r);
    quan_two_vec(check_threshold, l, r);
    quan_one_vec(channel_recons, l, r);
    //---------------------------
    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii])) / (2.0 * h_ins.rate);
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            initial_bits.assign(h_ins.vari_num, 0);
            for (unsigned ii = 0; ii < (h_ins.vari_num - h_ins.check_num); ii++)
            {
                initial_bits[ii] = d(gen);
            }
            generate_whole_bits(initial_bits);
            generate_codeword(initial_bits, codewords);
            //add noise
            noise_generator(codewords, cur_para);
            bool result = decoder_min(codewords, iteration, final_bits, initial_bits, l, r);
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
                //you can also output failed codewords here....
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}

bool Quantized_BP::decoder_min(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bit)
{

    //randomly generated codeword, without quantization
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<unsigned> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quan_size = 16;
    double full_return_v;
    unsigned full_return_c;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double final_codewords;
    std::vector<double> message_v;
    std::vector<unsigned> message_c;
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message_c.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj && msg_v2c[h_ins.edge_c[ii][kk]] != -1)
                    {
                        message_c.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if (message_c.size() != cur_dc - 1)
                {
                    std::cout << "not collecto enough data from variable nodes" << std::endl;
                    std::cout << message_c.size() << " " << cur_dc - 1 << std::endl;
                }
                full_return_c = check_node_operation_min(message_c, quan_size);
                msg_c2v[h_ins.edge_c[ii][jj]] = check_recons[cur_iter][full_return_c];
            }
        }

        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.clear();
                message_v.push_back(rx[ii]);
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj && msg_c2v[h_ins.edge_v[ii][kk]] != -1)
                    {
                        message_v.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                full_return_v = vari_node_operation(message_v);
                msg_v2c[h_ins.edge_v[ii][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii];
            message_v.clear();
            message_v.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message_v);
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
            std::cout<<final_codewords<<"  ";
        }
        std::cout<<std::endl;
        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bit[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    iteration = iteration + max_iter;
    return false;
}

bool Quantized_BP::decoder_min(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bit, int l, int r)
{
    //randomly generated codeword, with quantization
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<unsigned> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quan_size = 16;
    double full_return_v;
    unsigned full_return_c;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double final_codewords;
    std::vector<double> message_v;
    std::vector<unsigned> message_c;
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message_c.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj && msg_v2c[h_ins.edge_c[ii][kk]] != -1)
                    {
                        message_c.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if (message_c.size() != cur_dc - 1)
                {
                    std::cout << "not collecto enough data from variable nodes" << std::endl;
                    std::cout << message_c.size() << " " << cur_dc - 1 << std::endl;
                }
                full_return_c = check_node_operation_min(message_c, quan_size);
                msg_c2v[h_ins.edge_c[ii][jj]] = check_recons[cur_iter][full_return_c];
            }
        }
        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.clear();
                message_v.push_back(rx[ii]);
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj && msg_c2v[h_ins.edge_v[ii][kk]] != -1)
                    {
                        message_v.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                full_return_v = vari_node_operation(message_v, l, r);
                msg_v2c[h_ins.edge_v[ii][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii];
            message_v.clear();
            message_v.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message_v, l, r);
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bit[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    iteration = iteration + max_iter;
    return false;
}

void Quantized_BP::generate_whole_bits(std::vector<int> &inputbit)
{
    for (unsigned ii = 0; ii < h_ins.G_info[0].size(); ii++)
    {
        inputbit[h_ins.G_info[1][ii]] += inputbit[h_ins.G_info[0][ii]];
    }
    for (unsigned ii = 0; ii < inputbit.size(); ii++)
    {
        inputbit[ii] = inputbit[ii] % 2;
    }
}

void Quantized_BP::generate_codeword(std::vector<int> &inputbit, std::vector<double> &codewords)
{
    for (unsigned ii = 0; ii < inputbit.size(); ii++)
    {
        codewords[ii] = (inputbit[ii] == 1) ? -1.0 : 1.0;
    }
}

void Quantized_BP::decoder_min_track(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bit, int ind)
{
    //randomly generated codeword, without quantization
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<unsigned> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quan_size = 16;
    double full_return_v;
    unsigned full_return_c;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    std::vector<double> final_codewords(h_ins.vari_num, -1);
    std::vector<double> message_v;
    std::vector<unsigned> message_c;
    std::string filename = "trapping_set_decoder_min_no_interquan_" + std::to_string(ind) + ".txt";
    std::string llr_filename = "llr_decoder_min_no_interquan_" + std::to_string(ind) + ".txt";
    std::string tracking_filename = "trapping_set_info_" + std::to_string(ind) + ".txt";
    //code 0
    //std::vector<int> tracked_vari{101, 154, 215, 297, 302, 807, 879, 923, 928, 933, 982, 987, 1047, 1062, 1153, 1271};
    //std::vector <int> tracked_vari{297,982,928,302,987,933,879,215,923,101,807,1271,154,1153};
    //std::vector <int> tracked_vari{297, 923,  215,  879,  933,  987,  302,  928,  982};
    //std::vector <int> tracked_vari{297, 923,  215,  1153,  154,  1271,  101};
    //std::vector <int> tracked_vari{101, 297, 982,  928, 302,807};

    //code 11
    std::vector <int> tracked_vari{732,786,840,176,992,668,1261,405,833,349,1114};

    //code 13
    //std::vector <int> tracked_vari{59,183,309,314,940,994,999,1175,1229};
    
    std::ofstream track_handle(tracking_filename);
    std::ofstream ts_plot("ts_plot.txt");
    std::ofstream myfile(filename);
    std::ofstream llr_file(llr_filename);
    std::vector<int> cur_ts;
    std::vector<double> cur_llr;
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }


    //track information
    track_handle << "Input Varible Node LLR Information: " << std::endl<<std::endl;
    track_handle <<"|vn ind|";
    for (auto aa : tracked_vari)
    {
        track_handle << "   " << aa << "     |";
    }
    track_handle << std::endl<<"|";
    for (unsigned ii = 0; ii < tracked_vari.size()+1; ii++)
        track_handle << ":----:|";
    track_handle << std::endl;
    track_handle <<"|llr val|";
    for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
        track_handle << " " << rx[tracked_vari[ii]] << "  |";
    track_handle << std::endl;
    track_handle << "|llrsign|";
    for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
    {
        if(rx[tracked_vari[ii]]>0)
        {
            track_handle << " +  |";
        }
        else
        {
           track_handle << " -  |";
        }
    }
    track_handle << std::endl;
    track_handle << "|true sign|";
    for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
    {
        if(initial_bit[tracked_vari[ii]]==0)
        {
            track_handle << " +  |";
        }
        else
        {
           track_handle << " -  |";
        }
    }
    for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
    {
        if((initial_bit[tracked_vari[ii]]==0&&rx[tracked_vari[ii]]>0)||(initial_bit[tracked_vari[ii]]==1&&rx[tracked_vari[ii]]<0))
        {
            ts_plot<<"1 ";
        }
        else
        {
           ts_plot<<"0 ";
        }
    }
    ts_plot<<std::endl;
    track_handle << std::endl<<std::endl;


    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message_c.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj && msg_v2c[h_ins.edge_c[ii][kk]] != -1)
                    {
                        message_c.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if (message_c.size() != cur_dc - 1)
                {
                    std::cout << "not collecto enough data from variable nodes" << std::endl;
                    std::cout << message_c.size() << " " << cur_dc - 1 << std::endl;
                }
                full_return_c = check_node_operation_min(message_c, quan_size);
                msg_c2v[h_ins.edge_c[ii][jj]] = check_recons[cur_iter][full_return_c];
            }
        }

        track_handle << "Iteration: " << cur_iter << "" << std::endl<<std::endl;

        for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
        {
            track_handle << "Incomming messages into variable node: " << tracked_vari[ii]  << "  ";
            if (initial_bit[tracked_vari[ii]]==0)
            {
                track_handle<<"+"<<std::endl;
            }
            else
            {
                track_handle<<"-"<<std::endl;
            }
            
            int cur_dv = h_ins.vari_degreetable[tracked_vari[ii]];
            track_handle<<std::endl <<"|";
            for (int jj = 0; jj < cur_dv; jj++)
            {
                track_handle << "  " << h_ins.edge_relation[1][h_ins.edge_v[tracked_vari[ii]][jj]] << "  |";
            }
            track_handle <<std::endl<<"|";
            for (int jj = 0; jj < cur_dv; jj++)
            {
                track_handle << ":----:|";
            }
            track_handle << std::endl;
            track_handle <<"|";
            for (int jj = 0; jj < cur_dv; jj++)
            {
                track_handle << "  " << msg_c2v[h_ins.edge_v[tracked_vari[ii]][jj]] << "  |";
            }
            track_handle << std::endl
                         << std::endl;
        }


        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.clear();
                message_v.push_back(rx[ii]);
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj && msg_c2v[h_ins.edge_v[ii][kk]] != -1)
                    {
                        message_v.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                full_return_v = vari_node_operation(message_v);
                msg_v2c[h_ins.edge_v[ii][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii];
            message_v.clear();
            message_v.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords[ii] = vari_node_operation(message_v);
            if (final_codewords[ii] > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }

        //track information
        track_handle << "Iteration: " << cur_iter << "- Decision Making : " << std::endl
                     << std::endl;
        track_handle << "|vn ind|";
        for (auto aa : tracked_vari)
            track_handle << "  " << aa << "   |";
        track_handle << std::endl
                     << "|";
        for (unsigned ii = 0; ii < tracked_vari.size() + 1; ii++)
            track_handle << ":----:|";
        track_handle << std::endl
                     << "|llr val|";
        for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
            track_handle << " " << final_codewords[tracked_vari[ii]] << " |";
        track_handle << std::endl;
        track_handle << "|llr sign|";
        for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
        {
            if (final_codewords[tracked_vari[ii]] > 0)
            {
                track_handle << " +  |";
            }
            else
            {
                track_handle << " -  |";
            }
        }
        track_handle << std::endl;
        track_handle << "|true sign|";
        for (unsigned ii = 0; ii < tracked_vari.size(); ii++)
        {
            if (initial_bit[tracked_vari[ii]] == 0)
            {
                track_handle << " +  |";
            }
            else
            {
                track_handle << " -  |";
            }
            if (initial_bit[tracked_vari[ii]] != final_bits[tracked_vari[ii]])
            {
                ts_plot << "0 ";
            }
            else
            {
                ts_plot << "1 ";
            }
        }
        ts_plot << std::endl;
        track_handle << std::endl
                     << std::endl
                     << std::endl
                     << std::endl;


        //check sum
        iteration = iteration + cur_iter + 1;
        cur_llr.clear();
        cur_ts.clear();
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            if (final_bits[ii] != initial_bit[ii])
            {
                cur_ts.push_back(ii);
                cur_llr.push_back(final_codewords[ii]);
            }
        }
        for (auto aa : cur_ts)
        {
            myfile << aa << "  ";
        }
        for (auto aa : cur_llr)
        {
            llr_file << aa << "  ";
        }
        myfile << std::endl;
        llr_file << std::endl;
    }
    myfile.close();
    llr_file.close();
    track_handle.close();
    ts_plot.close();
    iteration = iteration + max_iter;
}

void Quantized_BP::decoder_min_track(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bit, int l, int r, int ind)
{
    quan_two_vec(vari_recons, l, r);
    quan_two_vec(vari_threshold, l, r);
    quan_two_vec(check_recons, l, r);
    quan_two_vec(check_threshold, l, r);
    quan_one_vec(channel_recons, l, r);
    //randomly generated codeword, with quantization
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<unsigned> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    std::vector<double> final_codewords(h_ins.vari_num, -1);
    unsigned quan_size = 16;
    double full_return_v;
    unsigned full_return_c;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    std::vector<double> message_v;
    std::vector<unsigned> message_c;
    std::string filename = "trapping_set_decoder_min_l_" + std::to_string(l) + "_r_" + std::to_string(r) + "_" + std::to_string(ind) + ".txt";
    std::string llr_filename = "llr_decoder_min_l_" + std::to_string(l) + "_r_" + std::to_string(r) + "_" + std::to_string(ind) + ".txt";
    std::ofstream myfile(filename);
    std::ofstream llr_file(llr_filename);
    std::vector<int> cur_ts;
    std::vector<double> cur_llr;
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message_c.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj && msg_v2c[h_ins.edge_c[ii][kk]] != -1)
                    {
                        message_c.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if (message_c.size() != cur_dc - 1)
                {
                    std::cout << "not collecto enough data from variable nodes" << std::endl;
                    std::cout << message_c.size() << " " << cur_dc - 1 << std::endl;
                }
                full_return_c = check_node_operation_min(message_c, quan_size);
                msg_c2v[h_ins.edge_c[ii][jj]] = check_recons[cur_iter][full_return_c];
            }
        }
        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.clear();
                message_v.push_back(rx[ii]);
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj && msg_c2v[h_ins.edge_v[ii][kk]] != -1)
                    {
                        message_v.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                full_return_v = vari_node_operation(message_v, l, r);
                msg_v2c[h_ins.edge_v[ii][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii];
            message_v.clear();
            message_v.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords[ii] = vari_node_operation(message_v, l, r);
            if (final_codewords[ii] > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        iteration = iteration + cur_iter + 1;
        //check sum
        iteration = iteration + cur_iter + 1;
        cur_llr.clear();
        cur_ts.clear();
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            if (final_bits[ii] != initial_bit[ii])
            {
                cur_ts.push_back(ii);
                cur_llr.push_back(final_codewords[ii]);
            }
        }
        for (auto aa : cur_ts)
        {
            myfile << aa << "  ";
        }
        for (auto aa : cur_llr)
        {
            llr_file << aa << "  ";
        }
        myfile << std::endl;
        llr_file << std::endl;
    }
    myfile.close();
    llr_file.close();
    iteration = iteration + max_iter;
}

void Quantized_BP::main_simulation_vertical_layered(const char ind[], const char suffix[], std::string filename, int l, int r, int layer_size)
{
    // random codeword generation, with quantization
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<> d(1, 0.5);
    if(strcmp(filename.data(),"all_zeros")!=0)
    {
        h_ins.Read_G_Matrix(filename);
    }
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    bool result;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
    //-----Basic Information Output-----
    std::cout << "Info: Random Codewords..." << std::endl;
    std::cout << "Quantization: l-" << l << " bit, r-" << r << " bit..." << std::endl;
    if (r!=0 || l!=0)
    {
        quan_two_vec(vari_recons, l, r);
        quan_two_vec(vari_threshold, l, r);
        quan_two_vec(check_recons, l, r);
        quan_two_vec(check_threshold, l, r);
        quan_one_vec(channel_recons, l, r);
    }

    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii])) / (2.0 * h_ins.rate);
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            initial_bits.assign(h_ins.vari_num, 0);
            if (strcmp(filename.data(),"all_zeros")!=0)
            {
                // if g matrix is loaded, then we do random encoding
                for (unsigned ii = 0; ii < (h_ins.vari_num - h_ins.check_num); ii++)
                {
                    initial_bits[ii] = d(gen);
                }
                generate_whole_bits(initial_bits);
            }
            generate_codeword(initial_bits, codewords);
            noise_generator(codewords, cur_para);
            if ( l!=0 || r!=0 )
            {
                // make internal message finite 
                result = decoder_min_vertical_layered(codewords ,iteration,  final_bits,  initial_bits, l, r, layer_size);
            }
            else
            {
                // make internal message full precision
                result = decoder_min_vertical_layered(codewords ,iteration,  final_bits,  initial_bits,layer_size);
            }            
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
                //you can also output failed codewords here....
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}

bool Quantized_BP::decoder_min_vertical_layered(std::vector<double> cwds ,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int l, int r, int layer_size)
{
    //randomly generated codeword, with quantization
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<unsigned> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quan_size = 16;
    double full_return_v;
    unsigned full_return_c;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double final_codewords;
    int dv,dc;
    int cur_check_node,cur_edge;
    std::vector<double> message_v;
    std::vector<unsigned> message_c;
    int total_layer_num;
    if (h_ins.vari_num%layer_size!=0)
    {
        std::cout<<"Info: variable number ("<<h_ins.vari_num<<") is invisible by layer_size ("<<layer_size<<"). please check again."<<std::endl;
        return false;
    }
    else
    {
        total_layer_num = h_ins.vari_num/layer_size;
    }
    
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }

 
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)];
        }
    }

    //start iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        // go through each layer 
        for (int cur_layer = 0 ; cur_layer<total_layer_num; cur_layer++)
        {

            int ini_vari_ind = cur_layer*layer_size;
            int end_vari_ind = (cur_layer+1)*layer_size-1;
            //update check node message
            for (int cur_vari = ini_vari_ind; cur_vari<=end_vari_ind; cur_vari++)
            {
                
                //std::cout<<cur_vari<<std::endl;
                dv=h_ins.vari_degreetable[cur_vari]; // get degree
                for (int check_ind=0; check_ind<dv; check_ind++)
                {
                    // go through all messages connected to this variable node
                    // and update message
                    cur_edge = h_ins.edge_v[cur_vari][check_ind];
                    cur_check_node = h_ins.edge_relation[1][cur_edge];
                    message_c.clear();
                    dc = h_ins.check_degreetable[cur_check_node];
                    for (int kk = 0; kk < dc; kk++)
                    {
                        if (h_ins.edge_c[cur_check_node][kk]!=cur_edge && msg_v2c[h_ins.edge_c[cur_check_node][kk]] != -1)
                        {
                            message_c.push_back(msg_v2c[h_ins.edge_c[cur_check_node][kk]]);
                        }
                    }
                    if (message_c.size() != dc - 1)
                    {
                        std::cout << "not collecto enough data from variable nodes";
                        std::cout << message_c.size() << " " << dc - 1 << std::endl;
                    }
                    full_return_c = check_node_operation_min(message_c, quan_size);
                    msg_c2v[cur_edge] = check_recons[cur_iter][full_return_c];
                }
            }
            // // update variabe node messages
            for (int cur_vari = ini_vari_ind; cur_vari <= end_vari_ind; cur_vari++)
            {
                int cur_dv = h_ins.vari_degreetable[cur_vari]; //find current dv
                for (int jj = 0; jj < cur_dv; jj++)
                {
                    message_v.clear();
                    message_v.push_back(rx[cur_vari]);
                    //std::cout<<"have been 2.5"<<std::endl;
                    for (int kk = 0; kk < cur_dv; kk++)
                    {
                        //collect data
                        if (kk != jj && msg_c2v[h_ins.edge_v[cur_vari][kk]] != -1)
                        {
                            message_v.push_back(msg_c2v[h_ins.edge_v[cur_vari][kk]]);
                        }
                    }
                    full_return_v = vari_node_operation(message_v, l, r);
                    msg_v2c[h_ins.edge_v[cur_vari][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
                }
            }
            
        }
        


        // go through each layer 
        // for (int cur_layer = 0 ; cur_layer<total_layer_num; cur_layer++)
        // {

        //     int ini_vari_ind = cur_layer*layer_size;
        //     int end_vari_ind = (cur_layer+1)*layer_size-1;

        //     // update variabe node messages
        //     for (int cur_vari = ini_vari_ind; cur_vari <= end_vari_ind; cur_vari++)
        //     {
        //         int cur_dv = h_ins.vari_degreetable[cur_vari]; //find current dv
        //         for (int jj = 0; jj < cur_dv; jj++)
        //         {
        //             message_v.clear();
        //             message_v.push_back(rx[cur_vari]);
        //             //std::cout<<"have been 2.5"<<std::endl;
        //             for (int kk = 0; kk < cur_dv; kk++)
        //             {
        //                 //collect data
        //                 if (kk != jj && msg_c2v[h_ins.edge_v[cur_vari][kk]] != -1)
        //                 {
        //                     message_v.push_back(msg_c2v[h_ins.edge_v[cur_vari][kk]]);
        //                 }
        //             }
        //             full_return_v = vari_node_operation(message_v, l, r);
        //             msg_v2c[h_ins.edge_v[cur_vari][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
        //         }
        //     }

        // }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii];
            message_v.clear();
            message_v.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message_v, l, r);
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bit[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    iteration = iteration + max_iter;
    return false;
}

bool Quantized_BP::decoder_min_vertical_layered(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bit, int layer_size)
{
    //randomly generated codeword, with quantization
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<unsigned> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    unsigned quan_size = 16;
    double full_return_v;
    unsigned full_return_c;
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double final_codewords;
    int dv, dc;
    int cur_check_node, cur_edge;
    std::vector<double> message_v;
    std::vector<unsigned> message_c;
    int total_layer_num;
    if (h_ins.vari_num % layer_size != 0)
    {
        std::cout << "Info: variable number (" << h_ins.vari_num << ") is invisible by layer_size (" << layer_size << "). please check again." << std::endl;
        return false;
    }
    else
    {
        total_layer_num = h_ins.vari_num / layer_size;
    }

    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }
    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)];
        }
    }

    //start iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {

        // go through each layer
        for (int cur_layer = 0; cur_layer < total_layer_num; cur_layer++)
        {

            int ini_vari_ind = cur_layer * layer_size;
            int end_vari_ind = (cur_layer + 1) * layer_size - 1;

            // update check node message
            for (int cur_vari = ini_vari_ind; cur_vari < end_vari_ind; cur_vari++)
            {
                dv = h_ins.vari_degreetable[cur_vari]; // get degree
                for (int check_ind = 0; check_ind < dv; check_ind++)
                {
                    // go through all messages connected to this variable node
                    // and update message
                    cur_edge = h_ins.edge_v[cur_vari][check_ind];
                    cur_check_node = h_ins.edge_relation[cur_edge][1];
                    message_c.clear();
                    dc = h_ins.check_degreetable[cur_check_node];
                    for (int kk = 0; kk < dc; kk++)
                    {
                        if (h_ins.edge_c[cur_check_node][kk] != cur_edge && msg_v2c[h_ins.edge_c[cur_check_node][kk]] != -1)
                        {
                            message_c.push_back(msg_v2c[h_ins.edge_c[cur_check_node][kk]]);
                        }
                    }
                    if (message_c.size() != dc - 1)
                    {
                        std::cout << "not collecto enough data from variable nodes" << std::endl;
                        std::cout << message_c.size() << " " << dc - 1 << std::endl;
                    }
                    full_return_c = check_node_operation_min(message_c, quan_size);
                    msg_c2v[cur_edge] = check_recons[cur_iter][full_return_c];
                }
            }

            // update variabe node messages
            for (int cur_vari = ini_vari_ind; cur_vari < end_vari_ind; cur_vari++)
            {
                int cur_dv = h_ins.vari_degreetable[cur_vari]; //find current dv
                for (int jj = 0; jj < cur_dv; jj++)
                {
                    message_v.clear();
                    message_v.push_back(rx[cur_vari]);
                    //std::cout<<"have been 2.5"<<std::endl;
                    for (int kk = 0; kk < cur_dv; kk++)
                    {
                        //collect data
                        if (kk != jj && msg_c2v[h_ins.edge_v[cur_vari][kk]] != -1)
                        {
                            message_v.push_back(msg_c2v[h_ins.edge_v[cur_vari][kk]]);
                        }
                    }
                    full_return_v = vari_node_operation(message_v);
                    msg_v2c[h_ins.edge_v[cur_vari][jj]] = threshold_quantization(vari_threshold[cur_iter], full_return_v);
                }
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii];
            message_v.clear();
            message_v.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message_v.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message_v);
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bit[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    iteration = iteration + max_iter;
    return false;
}

void Quantized_BP::main_simulation_horizontal_layered(const char ind[], const char suffix[],std::string filename, int l, int r, int layer_size)
{
    // random codeword generation, with quantization
    //-----------basic information output------------
    std::cout << "Layer Type: Horizontal" << std::endl;
    std::cout << "Layer Size: " << layer_size << std::endl;
    std::cout << "Codewords from: " << filename << std::endl;
    std::cout << "Quantization: l-" << l << " bit, r-" << r << " bit..." << std::endl;
    //-----------------------------------------------
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<> d(1, 0.5);
    if(strcmp(filename.data(),"all_zeros")!=0)
    {
        h_ins.Read_G_Matrix(filename);
    }
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    bool result;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
    //-----Basic Information Output-----
    
    if (r!=0 || l!=0)
    {
        quan_two_vec(vari_recons, l, r);
        quan_two_vec(vari_threshold, l, r);
        quan_two_vec(check_recons, l, r);
        quan_two_vec(check_threshold, l, r);
        quan_one_vec(channel_recons, l, r);
    }

    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii])) / (2.0 * h_ins.rate);
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_den_" + suffix + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            initial_bits.assign(h_ins.vari_num, 0);
            if (strcmp(filename.data(),"all_zeros")!=0)
            {
                // if g matrix is loaded, then we do random encoding
                for (unsigned ii = 0; ii < (h_ins.vari_num - h_ins.check_num); ii++)
                {
                    initial_bits[ii] = d(gen);
                }
                generate_whole_bits(initial_bits);
            }
            generate_codeword(initial_bits, codewords);
            noise_generator(codewords, cur_para);
            if ( l!=0 || r!=0 )
            {
                // make internal message finite 
                result = decoder_min_horizontal_layered(codewords ,iteration,  final_bits,  initial_bits, l, r, layer_size);
            }
            else
            {
                // make internal message full precision
                result = decoder_min_horizontal_layered(codewords ,iteration,  final_bits,  initial_bits,layer_size);
            }            
            if (!result)
            {
                /* fail */
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
                //you can also output failed codewords here....
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}

bool Quantized_BP::decoder_min_horizontal_layered(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bit, int layer_size)
{
    //randomly generated codeword, with quantization
    //Definition area
    unsigned quan_size = 16;
    std::vector<double> msg_c2v(h_ins.edge_num, 0);
    std::vector<double> msg_v2c, updated_c2v;
    std::vector<double> rx(h_ins.vari_num, -1);
    std::vector<double> rec_v = vari_sepcail_recons(quan_size); //  handmade reconstrcution
    std::vector<double> qua_c = check_sepcail_quan(quan_size);  //   handmade quantization
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    int starting_check, end_check;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double temp1;
    int temp2_int;
    int dc;
    int cur_check_node,cur_vari;
    int total_layer_num;
    //std::ofstream myfile ("llr.txt");
    if (h_ins.vari_num % layer_size != 0)
    {
        std::cout << "Info: variable number (" << h_ins.vari_num << ") is invisible by layer_size (" << layer_size << "). please check again." << std::endl;
        return false;
    }
    else
    {
        total_layer_num = h_ins.check_num / layer_size;
    }

    // Initialize Q value 
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }


    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        // go through each layer
        for (int cur_layer = 0; cur_layer < total_layer_num; cur_layer++)
        {

            starting_check = cur_layer * layer_size;
            end_check = (cur_layer + 1) * layer_size-1;
            for (cur_check_node = starting_check; cur_check_node <= end_check; cur_check_node++)
            {
                dc = h_ins.check_degreetable[cur_check_node];
                msg_v2c.assign(dc, -1);
                for (int ii = 0; ii < dc; ii++)
                {
                    //update rx                  
                    cur_vari = h_ins.edge_relation[0][h_ins.edge_c[cur_check_node][ii]];
                    rx[cur_vari] -= msg_c2v[h_ins.edge_c[cur_check_node][ii]];
                    // get c->v msg
                    temp1 = rx[cur_vari];
                    // quantization
                    if (cur_iter == 0)
                    {
                        temp2_int = cq_ins.quantizer[Quantize(temp1, quan_min, quan_max, interval, cq_ins.Partition_num)];
                    }
                    else
                    {
                        temp2_int = (int)threshold_quantization(vari_threshold[cur_iter-1], temp1);
                    }
                    //reconstruction
                    msg_v2c[ii] = rec_v[temp2_int];
                }
                //get updated c2v msgs
                updated_c2v = check_node_operation_fast_minsum(msg_v2c);
                for (int ii = 0; ii < dc; ii++)
                {
                    cur_vari = h_ins.edge_relation[0][h_ins.edge_c[cur_check_node][ii]];
                    // quantiation
                    temp2_int = (int)threshold_quantization(qua_c, updated_c2v[ii]);
                    // reconstruction & update c2v msges
                    msg_c2v[h_ins.edge_c[cur_check_node][ii]] = check_recons[cur_iter][temp2_int];
                    // update posterior
                    rx[cur_vari] += msg_c2v[h_ins.edge_c[cur_check_node][ii]];
                }
            }
           
        }
        int summ =0;
        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            if (rx[ii] > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
                summ++;
            }
        }

        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bit[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    //myfile.close();
    iteration = iteration + max_iter;
    return false;
}

bool Quantized_BP::decoder_min_horizontal_layered(std::vector<double> cwds ,int &iteration, std::vector<int> & final_bits, std::vector<int> & initial_bit, int l, int r, int layer_size)
{
    //randomly generated codeword, with quantization
    //Definition area
    double rpow, lpow, minpow, maxpow;
    rpow = pow(2, r);
    lpow = pow(2, (l - 1));
    minpow = (-1) * lpow;
    maxpow = lpow - (double)1.0 / (double)rpow;
    unsigned quan_size = 16;
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<double> msg_v2c, updated_c2v;
    std::vector<double> rx(h_ins.vari_num, -1);
    std::vector<double> rec_v = vari_sepcail_recons(quan_size); //  handmade reconstrcution
    std::vector<double> qua_c = check_sepcail_quan(quan_size);  //   handmade quantization
    double Max_val = cq_ins.Max_value;
    double Min_val = cq_ins.Min_value;
    int Partition = cq_ins.Partition_num;
    int starting_check, end_check;
    double quan_max = Max_val + (double)Max_val / Partition;
    double quan_min = Min_val - (double)Max_val / Partition;
    double interval = (double)(quan_max - quan_min) / Partition;
    double temp1;
    int temp2_int;
    int dc;
    int cur_check_node,cur_vari;
    int total_layer_num;
    if (h_ins.vari_num % layer_size != 0)
    {
        std::cout << "Info: variable number (" << h_ins.vari_num << ") is invisible by layer_size (" << layer_size << "). please check again." << std::endl;
        return false;
    }
    else
    {
        total_layer_num = h_ins.check_num / layer_size;
    }

    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = channel_recons[cq_ins.quantizer[Quantize(cwds[ii], quan_min, quan_max, interval, cq_ins.Partition_num)]];
    }

    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        // go through each layer
        for (int cur_layer = 0; cur_layer < total_layer_num; cur_layer++)
        {
            starting_check = cur_layer * layer_size;
            end_check = (cur_layer + 1) * layer_size;
            for (cur_check_node = starting_check; cur_check_node <= end_check; cur_check_node++)
            {
                dc = h_ins.check_degreetable[cur_check_node];
                msg_v2c.assign(dc, -1);
                for (int ii = 0; ii < dc; ii++)
                {
                    //update rx
                    cur_vari = h_ins.edge_relation[0][h_ins.edge_c[cur_check_node][ii]];
                    rx[cur_vari] -= msg_c2v[h_ins.edge_c[cur_check_node][ii]];
                    // get c->v msg
                    temp1 = rx[cur_vari];
                    // quantization
                    if (cur_iter == 0)
                    {
                        temp2_int = cq_ins.quantizer[Quantize(cwds[h_ins.edge_relation[0][h_ins.edge_c[cur_check_node][ii]]], quan_min, quan_max, interval, cq_ins.Partition_num)];
                    }
                    else
                    {
                        temp2_int = (int)threshold_quantization(vari_threshold[cur_iter - 1], temp1);
                    }
                    //reconstruction
                    msg_v2c[ii] = rec_v[temp2_int];
                }
                //gte updated c2v msgs
                updated_c2v = check_node_operation_fast_minsum(msg_v2c);

                for (int ii = 0; ii < dc; ii++)
                {
                    cur_vari = h_ins.edge_relation[0][h_ins.edge_c[cur_check_node][ii]];
                    // quantiation
                    temp2_int = (int)threshold_quantization(qua_c, updated_c2v[ii]);
                    // reconstruction & update c2v msges
                    msg_c2v[h_ins.edge_c[cur_check_node][ii]] = check_recons[cur_iter][temp2_int];
                    // update posterior
                    rx[cur_vari] += msg_c2v[h_ins.edge_c[cur_check_node][ii]];
                    if (rx[cur_vari]>maxpow)
                        rx[cur_vari] = maxpow;
                    if (rx[cur_vari]<minpow)
                        rx[cur_vari] = minpow;
                }
            }
        }
        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            if (rx[ii] > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bit[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    iteration = iteration + max_iter;
    return false;
}