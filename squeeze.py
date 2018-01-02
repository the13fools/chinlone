

fpath = "build/ex/"
fdir = "cube/"
fpath += fdir

with open(fpath + "0.txt") as f0:
    c0 = f0.readlines()
with open(fpath + "1.txt") as f1:
    c1 = f1.readlines()
with open(fpath + "2.txt") as f2:
    c2 = f2.readlines()


c0 = [x.strip() for x in open(fpath + "0.txt")] 
c1 = [x.strip() for x in open(fpath + "1.txt")] 
c2 = [x.strip() for x in open(fpath + "2.txt")] 

for i in range(len(c0)):
    with open(fpath + 'merged.txt', 'a') as the_file:
        the_file.write(c0[i] + " " + c1[i] + " " + c2[i] + "\n" )
