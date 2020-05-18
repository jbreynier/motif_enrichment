import csv
import exceptions
import os

class ReferenceGenome:
    def __init__(self, genome_fasta, genome_length_file, genome_include_file):
        for file_path in [genome_fasta, genome_length_file, genome_include_file]:
            if not os.path.isfile(file_path):
                message = ("Error: the following file path <{path}> "
                        "is missing").format(path=file_path)
                raise exceptions.IncorrectPathError(message)
        self.path_fasta = genome_fasta
        self.length_dict = generate_genome_dict(genome_length_file)
        self.path_include = genome_include_file

def generate_genome_dict(genome_length_file):
    '''Create dictionary of reference genome chromosome lengths'''
    length = {}
    with open(genome_length_file, "r") as genome_len:
        reader = csv.DictReader(genome_len, delimiter="\t")
        for row in reader:
            chrom = "chr" + row["chrom"]
            length[chrom] = int(row["length"])
    return length