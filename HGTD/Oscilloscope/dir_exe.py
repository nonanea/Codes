import os
import sys

dir = sys.argv[1]
plot_dir = dir
plot_dir = plot_dir.replace("Data/","Plots/")

if not os.path.exists(plot_dir):
    # os.mkdir(plot_dir)
    os.makedirs(plot_dir)

subdir_list = os.listdir(dir)

# print(subdir)
for subdir in subdir_list:
    # if "735_" not in subdir:    continue
    path = os.path.join(dir,subdir,"")
    # command = "python3 oscilloscope_raw.py -dir "+path
    # os.system(command)
    command = "root -l -b -q file_read.c'(\""+ path.replace("Data/","") + "\",20,50,-0.7)'"
    # print(command)
    os.system(command)