import sys
import random

if len(sys.argv) != 4:
    print("Usage: python match_size_random.py [sv_target] [sv_background] [out]")
    sys.exit()

sv1 = sys.argv[1]
sv2 = sys.argv[2]
outfile = sys.argv[3]

fin1=open(sv1, "r")
fout=open(outfile, "w")
fout.write("#chr\tstart\tend\tid\ttotal_available\n")

for line1 in fin1:
    if not line1.startswith("#"):
        SV1 = line1.split("\t")
        size1 = float(SV1[2])-float(SV1[1])
        list=[]
        fin2=open(sv2, "r")
        for line2 in fin2:
            if not line2.startswith("#"):
                SV2 = line2.split("\t")
                size2 = float(SV2[2])-float(SV2[1])
                if abs(size2/size1 -1) < 0.2:
                    list.append(SV2)
        if (len(list) > 0):
            match = random.sample(list, 1)
            stripped_match = [s.rstrip() for s in match[0]]
            fout.write("%s\t%s\n" % ("\t".join(stripped_match), len(list)))
        else:
            print("%s has no match with a size of %s " % (SV1, size1))
