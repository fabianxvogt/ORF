from re import finditer, compile

"""
CONSTANTS
"""
# File names
IN_FNAME = "input/in.fasta"
OUT_DNA_FNAME = "output/orfs_dna.fasta"
OUT_AA_FNAME = "output/orfs_aa.fasta"
IUPAC_FNAME = "data/iupac.txt"
# IUPAC codes (amino acid translation table)
IUPAC = dict(x.split(" ") for x in open(IUPAC_FNAME).read().splitlines())

"""
HELPER FUNCTIONS
"""
# RegEx search pattern to find all ORFs in a given string
# Params: min, max: Minimum and maximum length of amino acid
# Info: 
# - Begin of an ORF: 'ATG' (M)
# - End of an ORF: 'TAA', 'TAG' or 'TGA' (*)
# - Some number of amino acids (triplets) in between
def orf_regex(min, max):
    return compile(r"(?=(ATG(?:...){" + str(min) + r"," + str(max) + r"}?(?<=TAA|TAG|TGA)))")   

# Translate from DNA strings to UIPAC amino acid abbreviations
# Calculate start and end positions
# Params:
# - dna: Input DNA string OR the reverse complement
# - orf: ORF to translate and measure
def iupac_translate(orf):
    return ''.join([IUPAC[orf[i:i+3]] for i in range(0, len(orf)-2, 3)])

# Write a record to a file in FASTA format
# Params:
# - FASTA_file: FASTA filte to write to
# - descr: FASTA description line
# - data: FASTA data line 
def write_fasta_record(fasta_file, descr, data):
    fasta_file.write(">" + descr + "\n")
    fasta_file.write(data + "\n")

# Write ORFs to FASTA DNA file / Write amino acids to FASTA AA file
# Params:
# - dna_FASTA: FASTA DNA output file
# - aa_FASTA: FASTA AA output file
# - orf: ORF sequence
# - transl: Amino Acid abbreviations
# - orf_descr: ORF description line 
def write_orf(dna_fasta, aa_fasta, orf, start, end):
    orf_descr = "ORF_" + str(start) + "_" + str(end)
    write_fasta_record(dna_fasta, orf_descr, orf)
    write_fasta_record(aa_fasta, orf_descr, iupac_translate(orf))

# Generate the reverse complement of a DNA string
def get_revcomp(dna):
    return dna[::-1].translate(str.maketrans("ATGC","TACG"))

# Print the output files
def print_output():
    print("--> DNA sequences:")
    print(open(OUT_DNA_FNAME).read())
    print("--> Amino acids:")
    print(open(OUT_AA_FNAME).read())

"""
FIND ORFs  
"""
# Find Open Reading Frames (ORFs) in DNA sequences
# Params:
# - infile: FASTA-file containing the genome sequence
# - outfile_dna: FASTA-file containing DNA sequences
# - outfile_aa: FASTA-file containing amino acids
# - minlen: Minimal length of ORF
# - maxlen: Maximal length of ORF
def findORFs(infile, outfile_dna, outfile_aa, minlen, maxlen):
    # Open input file
    try:
        dna = open(infile).read().upper() 
    except (FileNotFoundError):
        print("ERROR: Input file could not be found!")
        return
    print("==> Reading input genome: " + dna)
    # Generate reverse complement
    revcomp = get_revcomp(dna)
    # Regex search pattern for amino acids
    regex_pattern = orf_regex(minlen, maxlen)
    # Find ORFs (RegEx matches) in input and reverse complement
    matches_input, matches_revcomp = tuple(finditer(regex_pattern, dna)), tuple(finditer(regex_pattern, revcomp))
    print("==> Found " + str(len(matches_input) + len(matches_revcomp)) + " ORFs. Saving to FASTA files...")
    # Open output files
    dna_file = open(outfile_dna, "w")
    aa_file = open(outfile_aa, "w")
    # Write input sequence ORFs
    for m in matches_input:
        write_orf(dna_file, aa_file, dna[m.regs[1][0]:m.regs[1][1]], m.regs[1][0] + 1, m.regs[1][1])
    # Write reverse complement ORFs
    for m in matches_revcomp:
        write_orf(dna_file, aa_file, revcomp[m.regs[1][0]:m.regs[1][1]], len(dna) - m.regs[1][0], len(dna) - m.regs[1][1] + 1)
    # Close output files
    dna_file.close()
    aa_file.close()
    print("==> DONE")
    return True

"""
TEST & MAIN
"""
# Test the 'findORFs' call with different ORF lengths
def test_findORFs(min, max):
    print("==> Testing ORF length " + str(min) + "-" + str(max) + "...")
    # FIND THE ORFs !!
    if (findORFs(IN_FNAME, OUT_DNA_FNAME, OUT_AA_FNAME, min, max)):
        print("\n==> Results:")
        print_output()
      
# MAIN CALL 
def main():
    test_findORFs(1,10)
    print("\n")
    test_findORFs(4,10)
main()