import os

file = "test.root"
command = "root -l -b -q tree_name.c"
command = command + "'" + "(\"" + file + "\")" + "'"
print(command)
os.system(command)
# os.system('root -l -b -q tree_name.c''(' + file+')\'')