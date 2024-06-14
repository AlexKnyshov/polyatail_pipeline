import pysam
import sys
import csv

def create_bed(bam_file):
    bed_file = bam_file.replace(".bam", ".bed")
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(bed_file, "w") as bed:
            for read in bam.fetch():
                readname = read.query_name
                read_length = read.query_length
                chrom = read.reference_name
                start = read.reference_start
                end = read.reference_end
                strand = "-" if read.is_reverse else "+"
                bed_line = f"{chrom}\t{start}\t{end}\t{readname}\t{strand}\t{read_length}\n"
                bed.write(bed_line)  

def get_soft_clipped(bam_file):
    output_file = bam_file.replace(".bam", ".csv")
    output_file = output_file.replace("stp2", "clipped_bases")
    bam_file = pysam.AlignmentFile(bam_file, "rb")
    with open(output_file, "w", newline="") as csvfile:
        fieldnames = ["readname", "clipped_bases"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for read in bam_file.fetch():
            soft_clipped_bases = 0
            for cigar in read.cigartuples:
                operation, length = cigar
                if operation == 4:  # Soft clipping operation
                    soft_clipped_bases += length

            writer.writerow({"readname": read.query_name, "clipped_bases": soft_clipped_bases})
              
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("provide the path to the BAM file.")
        sys.exit(1)

    bam_file_path = sys.argv[1]
    create_bed(bam_file_path)
    get_soft_clipped(bam_file_path)
