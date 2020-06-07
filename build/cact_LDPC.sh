#!/bin/sh
echo -n "enter design eb/no> "
read eb_no
for ii in  2.60
do
   cat Result_n_1296_k_648_den_$eb_no\_Para_$ii\0000_ind\_*.txt > Result_n_2193_k_1161_SNR_$ii\_Rate_0.50_den_$eb_no\.txt
done
