import re

# Load codesun (Amino Acid translation table)
def load_codesun():
    with open('data/codesun.txt') as f:
        return dict(x.split(" ") for x in f.read().splitlines()) 
# Reverse complement
def get_revcomp(dna):
    return dna[::-1].translate(str.maketrans("ATGC","TACG"))

# ==> findORFs
#   -> infile: FASTA-file containing the genome sequence
#   -> outfile_dna: FASTA-file containing DNA sequences
#   -> outfile_aa: FASTA-file containing amino acids
#   -> minlen: Minimal length of ORF
#   -> maxlen: Maximal length of ORF
def findORFs(infile, outfile_dna, outfile_aa, minlen, maxlen):
    # Open input file
    with open(infile) as f:
        dna = f.read().upper()
    print("==> Reading input genome: " + dna)
    # Open output files
    dna_file = open(outfile_dna, "w")
    aa_file = open(outfile_aa, "w")
    # Load codesun (Amino Acid translation table)
    codesun = load_codesun()
    # Generate reverse complement
    revcomp = get_revcomp(dna)
    # Regex search pattern for Amino acid triplets
    pattern = re.compile(r"(?=(ATG(?:...){" + str(minlen) + r"," + str(maxlen) + r"}?(?<=TAA|TAG|TGA)))")   
    # Find ORFs
    orfs = set(pattern.findall(dna) + pattern.findall(revcomp))
    print("==> Found "+str(len(orfs))+" ORFs! Begin translation using CodeSun...")
    # Loop through ORFs
    for orf in orfs:
        # translate amino acids
        transl = ""
        for i in range(0, len(orf)-2, 3):
            transl += codesun[orf[i:i+3]]
        # determine start & end positions
        pos = dna.find(orf) 
        if (pos == -1): 
            end = revcomp.find(orf) + 3
            start = end + len(orf) - 1
        else:
            start = pos + 1
            end = pos + len(orf)
        # write data
        fasta_descr = ">ORF_" + str(start) + "_" + str(end) + "\n"
        dna_file.write(fasta_descr)
        dna_file.write(orf + "\n")
        aa_file.write(fasta_descr)
        aa_file.write(transl + "\n")
    # close output files
    print("==> Translation done. Saving output files...")
    dna_file.close()
    aa_file.close()
    print("==> DONE")

def print_output():
    print("--> DNA sequences:")
    with open("output/orfs_dna.fasta") as f:
        print(f.read())
    print("--> Amino acids:")
    with open("output/orfs_aa.fasta") as f:
        print(f.read())

def test_findORFs(min, max):
    print("==> Testing ORF length " + str(min) + "-" + str(max) + "...")
    # FIND THE ORFs !!
    findORFs("input/in.fasta", "output/orfs_dna.fasta", "output/orfs_aa.fasta", min, max)
    print("\n")
    print("==> Results:")
    print_output()

# MAIN CALL 
def main():
    test_findORFs(1,10)
    print("\n")
    test_findORFs(4,10)
main()