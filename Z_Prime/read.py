import os

f = "sample_list.txt"

with open(f,"r") as text:
    content = text.readlines()
    # content = text.read().splitlines()

for line in content:
    print line
