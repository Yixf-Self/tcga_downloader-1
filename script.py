import json
import sys
import getopt
import os
import requests


metadata_file = "./metadata.json"
tss_file = "./tss_hg19_from_gencode.v19.tsv"

def main():

    tumor_type = None
    command = None


    #parser command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "t:c:", ["help", "tumor=", "command="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
        elif o in ("-t", "--tumor"):
            tumor_type = "TCGA-"+a
        elif o in ("-c", "--command"):
            if a.lower() in ("normal", "tumoral", "normalassociated"):
                command = a
            else:
                assert False, "unhandled value for command"
        else:
            assert False, "unhandled option"

    #output directory and check it not exists
    outdirectory = "./" + tumor_type + "-" + command
    downloaddirectory = outdirectory + "/download"
    assert not os.path.exists(outdirectory), "Directory " + outdirectory + " already exists"

    #load metadata json file
    with open(metadata_file) as data_file:
        jmeta = json.load(data_file)


    #filter by tumor type
    jtumor_all = [x for x in jmeta if x["cases"][0]["project"]["project_id"]==tumor_type]

    #distinct cases
    distinct_cases = list(set([x["cases"][0]["case_id"] for x in jtumor_all]))

    #primary tumors and patients
    primary_tumors = [x for x in jtumor_all
                      if int(x["cases"][0]["samples"][0]["sample_type_id"]) == 1]
    pt_patients = list(set([x["cases"][0]["case_id"] for x in primary_tumors]))

    #blood derived normals and patients
    blood_derived_normals = [x for x in jtumor_all
                             if int(x["cases"][0]["samples"][0]["sample_type_id"]) == 10]
    bdn_patients = list(set([x["cases"][0]["case_id"] for x in blood_derived_normals]))

    #solid tissue normals and patients
    solid_tissue_normals = [x for x in jtumor_all
                             if int(x["cases"][0]["samples"][0]["sample_type_id"]) == 11]
    stn_patients = list(set([x["cases"][0]["case_id"] for x in solid_tissue_normals]))

    #tumoral_associated with normal
    normal_ids = bdn_patients + stn_patients
    normal_associated_tumors = [x for x in primary_tumors if x["cases"][0]["case_id"] in normal_ids]

    #print recap
    print("There are", str(len(jtumor_all)), "files associated with", tumor_type,
          "associated with", str(len(distinct_cases)), "patients;", "\n",
          len(primary_tumors), "samples correspond to primary tumor, (",len(pt_patients),"patients)" "\n",
          len(blood_derived_normals), "samples correspond to blood derived normal, (",len(bdn_patients),"patients)", "\n",
          len(solid_tissue_normals), "samples correspond to solid tissue normal, (",len(stn_patients),"patients)")
    print("Normal samples are associated with", len(normal_associated_tumors), "tumoral samples.")

    #select files to be downloaded
    ids2download = []
    if command == "normal":
        ids2download = [x["file_id"] for x in solid_tissue_normals+blood_derived_normals]
    elif command == "tumoral":
        ids2download = [x["file_id"] for x in primary_tumors]
    else:
        ids2download = [x["file_id"] for x in normal_associated_tumors]

    #create folder
    os.makedirs(outdirectory)
    os.makedirs(downloaddirectory)

    #download files
    for i in range(0, len(ids2download)):
        f = ids2download[i]
        sys.stdout.write("\rDownloading file " + str(i) + " of " + str(len(ids2download)))
        sys.stdout.flush()
        fileout = downloaddirectory + "/" + f + ".fpkm.gz"
        while True:
            response = requests.get('https://gdc-api.nci.nih.gov/data/' + f)
            if response.status_code == 200:
                with open(fileout, 'wb') as o:
                    o.write(response.content)
                break

    print("File downloaded, unzipping...")
    os.system("(cd " + downloaddirectory + "; gunzip *)")

    print("Extract expressions...")
    for filename in os.listdir(downloaddirectory):
        os.system("(cd " + downloaddirectory +
                  "; grep  -v __ " + filename +
                  "| cut -f 2 > " + filename + ".expr)")
    filename = [filename for filename in os.listdir(downloaddirectory) if filename.endswith(".fpkm")][0]
    os.system("(cd " + downloaddirectory +
              "; grep -v __ "+ filename +" | cut -f 1  > genes.ids)")
    os.system("(cd " + downloaddirectory +
              "; rm *.fpkm)")

    print("Building the expression matrix...")
    expr_list = []
    expr_list.append([x.strip().partition(".")[0] for x in
                      open(downloaddirectory + "/genes.ids" , 'r').readlines()
                      ])
    for filename in filter(lambda x : x.endswith(".expr"),os.listdir(downloaddirectory)):
        expr_list.append(
            [x.strip() for x in
             open(downloaddirectory + "/" + filename, 'r').readlines()
             ]
        )

    print("Reading TSSs...")
    tss_regions = {}
    for line in open(tss_file, 'r').readlines():
        s = line.strip().split()
        name = s[0]
        chrom = s[1]
        pos = s[2]
        tss_regions[name] = (chrom, pos)

    print("Writing dataset...")
    outfile = open(outdirectory + "/" + tumor_type + "-expression-" + command + ".matrix", 'w')
    for x in zip(*expr_list):
        gene = x[0]
        if gene in tss_regions:
            outfile.write("\t".join([gene, tss_regions[gene][0], tss_regions[gene][1]]+list(x[1:])) + "\n")
    outfile.close()

    print("Removing download directory...")
    os.system("rm -rf " + downloaddirectory)

    print("Done.")


def usage():
    print("script.py -t TUMOR_TAG -c [normal|tumoral|normalAssociated]")


if __name__ == "__main__":

    main()
