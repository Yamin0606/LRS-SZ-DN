#!/usr/bin/python
import argparse
import os
current_path = os.path.abspath(".")
print(current_path)

#usage: python TEDenovo.py --geno M425-0_straglr_TR.txt  --parent M425-1_straglr_TR.txt  --o1 xx -o2 xx -o3 xx -o4 xx -o5 xx

def reverse(s):
  str = ""
  for i in s:
    str = i + str
  return str

def geno_separate(Geno, depth):
    """
    556.9(7);350.9(10);185.0(3)
    """
    genos = Geno.split(";")
    Lengths = [float(g.split("(")[0]) for g in genos]
    reads = [(g.split("(")[1]) for g in genos]
    Reads = [float(r.split(")")[0]) for r in reads]
    loc = [Reads.index(v) for v in Reads if v < float(depth)] #filter TR with reads <5
    loc.reverse()
    [Lengths.pop(v) for v in loc]
    return Lengths


def region_TR_geno(in_file):
    group = in_file.split("/")[-1].split("-")[1].split("_")[0]
    in_h = open(in_file, "r")
    Regions = set()
    RegionGeno = {}
    RegionLength = {}
    for line in in_h:
        line = line.strip()
        if line.startswith("#"):
            continue
        else:
            lines = line.split("\t")
            Chr, Start, End, Repeat, Geno = lines[:5]
            region = (Chr, Start, End)
            if group=="0":
                Lengths = geno_separate(Geno,5)
            elif group=="1" or group=="2":
                Lengths = geno_separate(Geno,3)
            if region not in Regions and len(Lengths) > 0:
                Regions.add(region)
                RegionGeno[region] = [Chr, Start, End, Repeat, Geno]
                RegionLength[region] = Lengths
    in_h.close()
    return RegionLength, RegionGeno

def denovo_TE(proband_file, parent_file, out_file1,out_file2,out_file3,out_file4,lengthRatio):
    id = proband_file.split("/")[-1].split("-")[0]
    out_h1 = open(out_file1, "w")
    out_h1.write("#Chr\tStart\tEnd\tRepeatUnit\tProbandGeno\tid\tParent1Geno\tParent2Geno\tproMax\tparMax\n")
    out_h2 = open(out_file2, "w")
    out_h2.write("#Chr\tStart\tEnd\tRepeatUnit\tProbandGeno\tid\tparmotif0\tparmotif1\n")
    out_h3 = open(out_file3, "w")
    out_h3.write("#Chr\tStart\tEnd\tRepeatUnit\tProbandGeno\tid\tgeno\tproMax\tparMax\n")
    out_h4 = open(out_file4, "w")
    out_h4.write("#Chr\tStart\tEnd\tRepeatUnit\tProbandGeno\tid\tproMax\n")

    lengthRatio = float(lengthRatio)

    probandLength, probandGeno = region_TR_geno(proband_file)

    Parents = parent_file.split(",")
    ParentNum = len(Parents)

    if ParentNum == 2:
        parentLength0, parentGeno0 = region_TR_geno(Parents[0])
        parentLength1, parentGeno1 = region_TR_geno(Parents[1])
        for region in probandLength:
            proLengths = probandLength[region]
            promotif = probandGeno[region][3]
            proMax = max(proLengths)
            if proMax > 100:
                probandRecord = probandGeno[region]
                if region in parentLength0 and region in parentLength1:
                    parLengths0 = parentLength0[region]
                    parLengths1 = parentLength1[region]
                    parmotif0 = parentGeno0[region][3]
                    parmotif1 = parentGeno1[region][3]
                    num_par0 = float(promotif == parmotif0)+float(promotif == reverse(parmotif0))
                    num_par1 = float(promotif == parmotif1)+float(promotif == reverse(parmotif1))
                    if num_par0 >0 and num_par1 >0:
                        parMax = max(parLengths0 + parLengths1)
                        if abs(proMax/parMax-1) > lengthRatio and parMax > 100:
                            geno0 = parentGeno0[region][-1]
                            geno1 = parentGeno1[region][-1]
                            out_h1.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("\t".join(probandRecord), id, geno0, geno1, proMax, parMax))
                    elif num_par0 >0:
                        parMax = max(parLengths0)
                        if abs(proMax/parMax-1) > lengthRatio and parMax > 100:
                            geno0 = parentGeno0[region][-1]
                            geno1 = "NA"
                            out_h1.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("\t".join(probandRecord), id, geno0, geno1, proMax, parMax))
                    elif num_par1 >0:
                        parMax = max(parLengths1)
                        if abs(proMax/parMax-1) > lengthRatio and parMax > 100:
                            geno0 = "NA"
                            geno1 = parentGeno1[region][-1]
                            out_h1.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("\t".join(probandRecord), id, geno0, geno1, proMax, parMax))
                    else:
                        out_h2.write("%s\t%s\t%s\t%s\n" % ("\t".join(probandRecord), id, parmotif0, parmotif1))
                elif region in parentLength0 or region in parentLength1:
                    if region in parentLength0:
                        parLengths = parentLength0[region]
                        parmotif0 = parentGeno0[region][3]
                        id_par = Parents[0].split("/")[-1].split("_")[0]
                        if parmotif0 == promotif or parmotif0 == reverse(promotif):
                            geno = parentGeno0[region][-1]
                            parMax = max(parLengths)
                            if abs(proMax/parMax-1) > lengthRatio and parMax > 100:
                                out_h3.write("%s\t%s\t%s\t%s\t%s\n" % ("\t".join(probandRecord), id_par, geno, proMax, parMax))
                    elif region in parentLength1:
                        parLengths = parentLength1[region]
                        parmotif1 = parentGeno1[region][3]
                        id_par = Parents[1].split("/")[-1].split("_")[0]
                        if parmotif1 == promotif or parmotif1 == reverse(promotif):
                            geno = parentGeno1[region][-1]
                            parMax = max(parLengths)
                            if abs(proMax/parMax-1) > lengthRatio and parMax > 100:
                                out_h3.write("%s\t%s\t%s\t%s\t%s\n" % ("\t".join(probandRecord), id_par, geno, proMax, parMax))
                else:
                    out_h4.write("%s\t%s\t%s\n" % ("\t".join(probandRecord), id, proMax))
            else:
                print("%s,is smaller than 100 bp with the max allele = %s\n" % (probandGeno[region],proMax) )
    elif ParentNum == 1:
        print("Two parents are needed to call de novo TR")
    out_h1.close()
    out_h2.close()
    out_h3.close()
    out_h4.close()

def overlap(out_h4,parent_file,outfile):
    Parents = parent_file.split(",")
    fout=open(outfile, "w")
    fout.write("#chr\tstart\tend\tnum_olp\tregion_olp\n")
    fin1 = open(out_h4,"r")
    for line1 in fin1:
        if not line1.startswith("#"):
            chr1, s1, e1 = line1.split("\t")[:3]
            list=[]
            Regions = set()
            for bed2 in Parents:
                fin2=open(bed2, "r")
                for line2 in fin2:
                    if not line2.startswith("#"):
                        chr2, s2, e2 = line2.split("\t")[:3]
                        region = (chr2,s2,e2)
                        if region not in Regions:
                            Regions.add(region)
                            if (chr1 == chr2):
                                if e1<s2 or s1>e2:
                                    pass
                                elif (s1>s2 and e1<e2) or (s1<s2 and e1>s2) or (s1<e2 and e1>e2) or (s1<s2 and e1>e2):
                                    list.append(line2)
            if len(list) != 0:
                fout.write("%s\t%s\t%s\n" % ("\t".join((chr1, s1, e1)), len(list), list))



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--geno", help="The input file of TR genotype.")
    parser.add_argument("-p", "--parent", help="The input parent file of TR genotype.")
    parser.add_argument("-o1", "--out1", help="The output file of de novo TR shared with both parents.")
    parser.add_argument("-o2", "--out2", help="The output file of TR in the same region but with different motif.")
    parser.add_argument("-o3", "--out3", help="The output file of de novo TR shared only with  one of the parents.")
    parser.add_argument("-o4", "--out4", help="The output file of de novo TR that do not exist in parents.")
    parser.add_argument("-r", "--lengthRagio", default=0.5, help="The threshold of length ratio.")
    parser.add_argument("-o5", "--out5", help="The output file of overlapping.")
    args = parser.parse_args()
    denovo_TE(args.geno, args.parent, args.out1, args.out2, args.out3, args.out4, args.lengthRagio)
    overlap(args.out4, args.parent, args.out5)

if __name__ == "__main__":
    main()

