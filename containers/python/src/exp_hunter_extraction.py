# import json


# def extract_info_from_json(file_path):
#     with open(file_path, "r") as file:
#         data = json.load(file)

#     result = {"Metadata": {}, "Variants": []}
#     allele_count = data["LocusResults"]["SNV_AND_STR"]["AlleleCount"]
#     coverage = data["LocusResults"]["SNV_AND_STR"]["Coverage"]
#     sex = data["SampleParameters"]["Sex"]

#     result["Metadata"].update(
#         {
#             "Allele Count": allele_count,
#             "Coverage": coverage,
#             "Sex": sex,
#         }
#     )
#     print("Allele Count:", allele_count)
#     print("Coverage:", coverage)
#     print("Sex:", sex)
#     print("Variants:")
#     for variant_id, variant_info in data["LocusResults"]["SNV_AND_STR"][
#         "Variants"
#     ].items():
#         genotype = variant_info["Genotype"]
#         variant_id = variant_info["VariantId"]
#         variant_type = variant_info["VariantType"]
#         repeat_unit = (
#             variant_info["RepeatUnit"] if "RepeatUnit" in variant_info else "/"
#         )
#         result["Variants"].append(
#             {
#                 "VariantId": variant_id,
#                 "Genotype": genotype,
#                 "VariantType": variant_type,
#                 "RepeatUnit": repeat_unit,
#             }
#         )
#         print("\tVariant ID:", variant_id)
#         print("\tGenotype:", genotype)
#         print("\tVariant Count:", variant_id)
#         print("\tVariant Type:", variant_type)
#         print("\tRepeat Unit:", repeat_unit)
#         print()
#     return result


# def extract_info_from_vcf(vcf_file):
#     variant_info = []

#     with open(vcf_file, "r") as f:
#         for line in f:
#             if line.startswith("#"):
#                 continue
#             fields = line.strip().split("\t")
#             chrom = fields[0]
#             pos = fields[1]
#             ref = fields[3]
#             alt = fields[4]
#             info = fields[7]
#             format_fields = fields[8].split(":")
#             sample_fields = fields[9].split(":")

#             info_dict = {}
#             info_items = info.split(";")
#             for item in info_items:
#                 key, value = item.split("=")
#                 info_dict[key] = value

#             genotype_index = format_fields.index("GT")
#             genotype = sample_fields[genotype_index]

#             repeat_unit = info_dict.get("RU", "N/A")
#             repeat_copy_number = info_dict.get("REPCN", "N/A")
#             locus_coverage = sample_fields[format_fields.index("LC")]

#             variant_info.append(
#                 {
#                     "Chrom": chrom,
#                     "Pos": pos,
#                     "Ref": ref,
#                     "Alt": alt,
#                     "Genotype": genotype,
#                     "RepeatUnit": repeat_unit,
#                     "RepeatCopyNumber": repeat_copy_number,
#                     "LocusCoverage": locus_coverage,
#                 }
#             )
#             print(variant_info)


# json_file_path = "exp_hunter_output/repeats.json"
# vcf_file = "exp_hunter_output/repeats.vcf"

import json
import sys


def extract_info_from_json(file_path, out_path):
    with open(file_path, "r") as file:
        data = json.load(file)

    result = {"Metadata": {}, "Locus Results": []}

    result["Metadata"].update(
        {
            "Sample ID": data["SampleParameters"]["SampleId"],
            "Sex": data["SampleParameters"]["Sex"],
        }
    )

    for locus_id, locus_info in data["LocusResults"].items():
        allele_count = locus_info["AlleleCount"]
        coverage = locus_info["Coverage"]
        fragment_length = locus_info["FragmentLength"]
        read_length = locus_info["ReadLength"]

        print(f"Locus ID: {locus_id}")
        print("Allele Count:", allele_count)
        print("Coverage:", coverage)
        print("Fragment Length:", fragment_length)
        print("Read Length:", read_length)

        # TODO: check if variants is always a single item
        print("Variants:")
        for variant_id, variant_info in locus_info["Variants"].items():
            reference_region = variant_info["ReferenceRegion"]
            variant_id = variant_info["VariantId"]
            variant_type = variant_info["VariantType"]
            repeat_unit = (
                variant_info["RepeatUnit"] if "RepeatUnit" in variant_info else "/"
            )
            genotype = variant_info.get("Genotype", "/")
            genotype_confidence_interval = variant_info.get("GenotypeConfidenceInterval", "/")
            print("\tVariant ID:", variant_id)
            print("\tReference Region:", reference_region)
            print("\tVariant Type:", variant_type)
            print("\tRepeat Unit:", repeat_unit)
            print()
            result["Locus Results"].append(
                {
                    "Locus ID": locus_id,
                    # "Variant ID": variant_id,
                    "Allele Count": allele_count,
                    "Coverage": coverage,
                    "Fragment Length": fragment_length,
                    "Read Length": read_length,
                    "Reference Region": reference_region,
                    "Variant Type": variant_type,
                    "Repeat Unit": repeat_unit,
                    "Genotype": genotype,
                    "Genotype Confidence Interval": genotype_confidence_interval,
                }
            )
    with open(out_path, "w") as file:
    # with open("result_from_json.json", "w") as file:
        json.dump(result, file, indent=4)
    return result


# Example usage:
# file_path = "/home/nevena/random/persida/prph_nf/input/DEMO-1X-1_exph.json"

def main():
    file_path = sys.argv[1]
    out_path = sys.argv[2]
    extract_info_from_json(file_path, out_path)

if __name__ == "__main__":
    main()


# file_path = sys.argv[1]
# extract_info_from_json(file_path)
