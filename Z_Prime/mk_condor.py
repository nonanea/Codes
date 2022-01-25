import os
import sys

f_name = ["mc16.988603.MGPy8EG_30nlo_Leptophilicmutau_2vZp011",
"mc16.988604.MGPy8EG_30nlo_Leptophilicmutau_2vZp013",
"mc16.988605.MGPy8EG_30nlo_Leptophilicmutau_2vZp015",
"mc16.988606.MGPy8EG_30nlo_Leptophilicmutau_2vZp017",
"mc16.988607.MGPy8EG_30nlo_Leptophilicmutau_2vZp019",
"mc16.988608.MGPy8EG_30nlo_Leptophilicmutau_2vZp023",
"mc16.988609.MGPy8EG_30nlo_Leptophilicmutau_2vZp027",
"mc16.988610.MGPy8EG_30nlo_Leptophilicmutau_2vZp031",
"mc16.988611.MGPy8EG_30nlo_Leptophilicmutau_2vZp035",
"mc16.988612.MGPy8EG_30nlo_Leptophilicmutau_2vZp039",
"mc16.988613.MGPy8EG_30nlo_Leptophilicmutau_2vZp042",
"mc16.988614.MGPy8EG_30nlo_Leptophilicmutau_2vZp045",
"mc16.988615.MGPy8EG_30nlo_Leptophilicmutau_2vZp048",
"mc16.988616.MGPy8EG_30nlo_Leptophilicmutau_2vZp051",
"mc16.988617.MGPy8EG_30nlo_Leptophilicmutau_2vZp054",
"mc16.988618.MGPy8EG_30nlo_Leptophilicmutau_2vZp057",
"mc16.988619.MGPy8EG_30nlo_Leptophilicmutau_2vZp060",
"mc16.988620.MGPy8EG_30nlo_Leptophilicmutau_2vZp063",
"mc16.988621.MGPy8EG_30nlo_Leptophilicmutau_2vZp066",
"mc16.988622.MGPy8EG_30nlo_Leptophilicmutau_2vZp069",
"mc16.988623.MGPy8EG_30nlo_Leptophilicmutau_2vZp072",
"mc16.988624.MGPy8EG_30nlo_Leptophilicmutau_2vZp075",
"mc16.988690.MGPy8EG_30nlo_Leptophilicmutau_2vZp010"]

for name in f_name:
    f = open("sub_list/"+name+".sub","a")
    f.seek(0)
    f.truncate()

    f.write("Universe = vanilla\nExecutable = execute.sh\n")

    # Output = the.out
    # Error = the.err
    f.write("Output = logs/"+ name + ".out\n")
    f.write("Error = logs/"+ name + ".err\n")
    f.write("Queue")

    os.system("condor_submit "+name+".sub")
    # print("condor_submit "+name+".sub")