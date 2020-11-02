import sys
import csv
import gzip
from snakemake.utils import min_version

csv.field_size_limit(100000000)
csv.field_size_limit()


def main(mode, out_file, file_list  ):
    with open(out_file, 'w') as out:

        for file in file_list:

            with gzip.open(file, mode="rt") as f:
                
                if mode=="Isoform" or mode=="Gene":
                    header = ["Sample", mode, "TpM", "Read_Counts"]
                elif mode=="Isoform_salmon":
                    header = ["Sample", 'Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads']

                writer = csv.DictWriter(out, fieldnames=header, extrasaction='ignore', delimiter="\t")
                writer.writeheader()

                sample = file.split("/")[-1].split(".")[0]
                reader = csv.DictReader(f, delimiter="\t")

                for row in reader:

                    row["Sample"] = sample
                    writer.writerow(row)
                
if __name__ == '__main__':
    main(snakemake.params["feature"], snakemake.output["merged"],  snakemake.input["files"])

        
        
    
