import sys

if len(sys.argv) != 4:
    print("Usage: python overlap_hotspot_sv_region.py [sv] [hotspot] [outfile]  > log.txt 2>&1 ")
    sys.exit()

sv = sys.argv[1]
hotspot = sys.argv[2]
outfile = sys.argv[3]

fin1=open(sv, "r")
fin2=open(hotspot, "r")
num1=len(fin1.readlines())
num2=len(fin2.readlines())-1
num3=num1*num2+1
print("There should be %s*%s+1=%s lines in the log file" % (num1,num2,num3))

fout=open(outfile, "w")
fout.write("chr\tstart\tend\tnum_overlap\tregion_overlap\n")
fin1=open(sv, "r")
for line1 in fin1:
    if not line1.startswith("#"):
        SV = line1.split("\t")[:3]
        chr1 = SV[0]
        s1 = float(SV[1])
        e1 = float(SV[2])
        list=[]
        fin2=open(hotspot, "r")
        for line2 in fin2:
            if not line2.startswith("#"):
                Hotspot = line2.split("\t")
                chr2 = Hotspot[0]
                s2 = float(Hotspot[1])
                e2 = float(Hotspot[2])
                if (chr1 == chr2):
                    if e1<s2 or s1>e2:
                        #print("not overlap")
                        pass
                    elif (s1>s2 and e1<e2) or (s1<s2 and e1>s2) or (s1<e2 and e1>e2) or (s1<s2 and e1>e2):
                        #print("overlap")
                        list.append(line2)
                else:
                    pass
                    #print(woifferent chr")
        stripped_SV = [s.rstrip() for s in SV]
        if len(list) == 0:
            fout.write("%s\t%s\t%s\n" % ("\t".join(stripped_SV), 0, 0))
        else:
            fout.write("%s\t%s\t%s\n" % ("\t".join(stripped_SV), len(list), list))





