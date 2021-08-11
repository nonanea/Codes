import os

full_set=[410472,364250,364230,345708,410218,410219,410220,410157,413008,346345,341471,341488,342284,342285,364100,364101,364102,364103,364104,364105,364106,364107,364108,364109,364110,364111,364112,364113,364114,364115,364116,364117,364118,364119,364120,364121,364122,364123,364124,364125,364126,364127,450578]

# full_set = ["k"]

list = os.listdir(os.getcwd())
dir = []

prefix = "tree"

for i in range(len(list)):
    if prefix in list[i]:
        dir.append(list[i])

for i in range(len(full_set)):
    flag = 0
    for j in range(len(dir)):
        if str(full_set[i]) in dir[j]:
            flag = 1
    
    if flag == 0:
        print(full_set[i])

# print(dir)