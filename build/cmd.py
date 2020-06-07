import os
import subprocess

H_file = "80211_irr_648_1296.txt"
total_target_num=100
G_file = "80211_G.txt"
l = 5
r = 5
lut_suffix="092"
ebno=1.40
layer_size=54
target_wrong_number=100
ind = 0 

"""para_cmd="./runner -Hfile "+  H_file +" -Gfile  " + G_file + " -suf " + lut_suffix + " -ebno  "+str(ebno)+" -tarerr " + str(total_target_num) +" -ind " + str(ind) +" -layersize "+ str(layer_size)  +" -l " + str(l) + " -r " + str(r)"""
para_cmd="parallel -j 15 ./runner -Hfile "+  H_file +" -Gfile  " + G_file + " -suf " + lut_suffix + " -ebno  {1} -tarerr  5  -ind {2} -layersize "+ str(layer_size)  +" -l " + str(l) + " -r " + str(r)+" ::: 2.4 2.5 2.6 ::: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19" 
subprocess.call(para_cmd,shell=True)
