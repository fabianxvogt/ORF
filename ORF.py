import re

# Regex search pattern for Amino acid triplets
pattern = re.compile(r"(?=(ATG(?:...){" + min + r"," + max+ r"}?(?<=TAA|TAG|TGA)))")   

# Reverse complement
def revcomp(dna):
    return dna[::-1].translate(str.maketrans("ATGC","TACG"))

# Get ORFS
def orfs(dna):
    return set(pattern.findall(dna) + pattern.findall(revcomp(dna)))

# Load codesun (Amino Acid translation table)
def load_codesun():
    with open('data/codesun.txt') as f:
        return dict(x.split(" ") for x in f.read().splitlines())

# ==> findORFs
# -> infile: FASTA-file containing the genome sequence
# -> outfile_dna: FASTA-file containing DNA sequences
# -> outfile_aa: FASTA-file containing amino acids
# -> minlen: Minimal length of ORF
# -> maxlen: Maximal length of ORF
def findORFs(infile, outfile_dna, outfile_aa, minlen, maxlen):
    
    pass
