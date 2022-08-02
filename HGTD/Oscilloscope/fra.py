import os

indir = "Plots/"

# f_list = os.listdir(indir)

f_out = open("fra.csv","a")
f_out.seek(0)
f_out.truncate()

Pos = [[],[]]
Var = [[],[],[],[]]

for i in range(4):
    f = indir + "amplitude_" + str(i+1) + ".csv"
    print f
    with open(f,"r") as text:
        content = text.read().splitlines()

        for var in content:
            var = var.split(",")
            Var[i].append(float(var[2]))
            if i == 0:
                Pos[0].append(var[0])
                Pos[1].append(var[1])

for i,var in enumerate(Var[0]):
    # print Pos[0][i],Pos[1][i],Var[0][i]
    fra1 = (Var[0][i]-Var[1][i])/(Var[0][i]+Var[1][i])
    # fra2 = (Var[2][i]-Var[3][i])/(Var[2][i]+Var[3][i])
    fra2 = (Var[0][i]+Var[2][i]-Var[1][i]-Var[3][i])/(Var[0][i]+Var[2][i]+Var[1][i]+Var[3][i])
    f_out.write( Pos[0][i] + "," + Pos[1][i] + "," + str(fra1) + "," + str(fra2) + "\n" )

f_out.close()
