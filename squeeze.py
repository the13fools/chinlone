

fpath = "test_data/"
fdir = "tet_30/"
fpath += fdir


c0 = [x.strip() for x in open(fpath + "0.txt")] 
c1 = [x.strip() for x in open(fpath + "1.txt")] 
c2 = [x.strip() for x in open(fpath + "2.txt")] 
#c3 = [x.strip() for x in open(fpath + "3.txt")] 

for i in range(len(c0)):
    with open(fpath + 'merged.txt', 'a') as the_file:
     #   the_file.write(c0[i] + " " + c1[i] + " " + c2[i] + " " + c3[i] + "\n" )
        the_file.write(c0[i] + " " + c1[i] + " " + c2[i] +  "\n" )
