from io import TextIOWrapper

# Reads FASTA file from path into dictionary
# Params:
# - fasta_file_path: FASTA file path to read from
def fasta_to_dict(fasta_file_path: str):
    d = {}
    for l in open(fasta_file_path).read().splitlines():
        if l.startswith(">"):
            k = l[1:]
            d[k] = ""
        else: 
            d[k] += l
    return d

# Write a record to a file in FASTA format
# Params:
# - FASTA_file: FASTA file to write to
# - descr: FASTA description line
# - data: FASTA data line 
def write_fasta_record(fasta_file: TextIOWrapper, descr: str, data: str): 
    [fasta_file.write(x) for x in (">" + descr + "\n", data + "\n")]
    return fasta_file