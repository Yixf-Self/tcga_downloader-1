import json
import sys
import getopt


metadata_file = "./metadata.json"
tumor_type = "TCGA-BRCA"

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:v", ["help", "output="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            outputt = a
        else:
            assert False, "unhandled option"

    with open(metadata_file) as data_file:
        jmeta = json.load(data_file)

    print(len(jmeta))

    #filter by tumor type
    jtumor_all = [x for x in jmeta if x["cases"][0]["project"]["project_id"]==tumor_type]

    #distinct cases
    distinct_cases = list(set([x["cases"][0]["case_id"] for x in jtumor_all]))

    primary_tumors = [x for x in jtumor_all
                      if int(x["cases"][0]["samples"][0]["sample_type_id"]) == 1]
    pt_patients = list(set([x["cases"][0]["case_id"] for x in primary_tumors]))

    blood_derived_normals = [x for x in jtumor_all
                             if int(x["cases"][0]["samples"][0]["sample_type_id"]) == 10]
    bdn_patients = list(set([x["cases"][0]["case_id"] for x in blood_derived_normals]))

    solid_tissue_normals = [x for x in jtumor_all
                             if int(x["cases"][0]["samples"][0]["sample_type_id"]) == 11]
    stn_patients = list(set([x["cases"][0]["case_id"] for x in solid_tissue_normals]))

    normal_ids = bdn_patients + stn_patients
    normal_associated_tumors = [x for x in primary_tumors if x["cases"][0]["case_id"] in normal_ids]

    print("There are", str(len(jtumor_all)), "file associated with", tumor_type,
          "associated with", str(len(distinct_cases)), "patients;", "\n",
          len(primary_tumors), "samples correspond to primary tumor, (",len(pt_patients),"patients)" "\n",
          len(blood_derived_normals), "samples correspond to blood derived normal, (",len(bdn_patients),"patients)", "\n",
          len(solid_tissue_normals), "samples correspond to solid tissue normal, (",len(stn_patients),"patients)")
    print("Normal samples are associated with", len(normal_associated_tumors), "tumoral samples.")



if __name__ == "__main__":

    main()
