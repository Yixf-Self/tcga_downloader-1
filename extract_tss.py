#Extract TSS from a GENCODE GFF/GTF3 file
import sys

def main():
    if len(sys.argv) < 3:
        usage()
        sys.exit(2)

    input_name = sys.argv[1]
    output_name = sys.argv[2]

    #name -> (chr, pos)
    tss = {}
    for line in (open(input_name, "r").readlines()):
        if line.startswith("#"):
            pass
        else:
            s = line.strip().split()
            chrom = s[0]
            strand = s[6]
            pos = int(s[3]) if strand == "+" else int(s[4])
            name = [x[len("gene_id")+1:] for x in s[8].split(";") if x.startswith("gene_id=")][0].partition(".")[0]
            if name not in tss:
                tss[name] = (chrom, pos)
            elif strand == "+":
                if pos < tss[name][1]:
                    tss[name] = (chrom, pos)
            else:
                if pos > tss[name][1]:
                    tss[name] = (chrom, pos)

    outfile = open(output_name, "w")
    for k in tss.keys():
        outfile.write("\t".join([k,tss[k][0],str(tss[k][1])]) + "\n")
    outfile.close()


def usage():
    print("extract_tss input.gff output.txt")


if __name__ == "__main__":

    main()