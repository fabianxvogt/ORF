from re import finditer, compile
from fasta_utils import fasta_to_dict, write_fasta_record

"""
CONSTANTS
"""
DOT_FASTA = ".fasta"
# File names
IN_FNAME = "input/in" + DOT_FASTA
OUT_DNA_FNAME = "output/orfs_dna"
OUT_AA_FNAME = "output/orfs_aa"
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
def orf_regex(min, max):#
    return compile(r"(?=(ATG([ACTG]{3}(?<!ATG|TAA|TAG|TGA)){"+str(min-1)+r","+str(max-1)+r"}(TAA|TAG|TGA)))")

# Translate from DNA strings to UIPAC amino acid abbreviations
# Calculate start and end positions
# Params:
# - dna: Input DNA string OR the reverse complement
# - orf: ORF to translate and measure
def iupac_translate(orf):
    return ''.join([IUPAC[orf[i:i+3]] for i in range(0, len(orf)-2, 3)])
        
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
# Param: dna: DNA string to get the Revcomp from
def get_revcomp(dna):
    return dna[::-1].translate(str.maketrans("ATGC","TACG"))

# Print the output files
def print_output(output_files): 
    for gen in output_files.values():
        for f in gen:
            print(f)
            print(open(f).read())

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
        input_genoms = fasta_to_dict(infile)
    except (FileNotFoundError):
        print("ERROR: Input file could not be found!")
        return
    print("==> " + str(len(input_genoms)) + " Genoms found. Begin searching for ORFs...")
    # Variable to store output file names
    out_files = {} 
    for i, (genome_name, dna) in enumerate(input_genoms.items()):
        print("==> Reading input genome "+str(i+1)+" ("+genome_name+"): " + dna)
        # Generate reverse complement
        revcomp = get_revcomp(dna)
        print("==> Reverse Complement "+str(i+1)+" ("+genome_name+"): " + revcomp)
        # Regex search pattern for amino acids
        regex_pattern = orf_regex(minlen, maxlen)
        # Find ORFs (RegEx matches) in input and reverse complement
        matches_input, matches_revcomp = tuple(finditer(regex_pattern, dna)), tuple(finditer(regex_pattern, revcomp))
        print("==> Found " + str(len(matches_input) + len(matches_revcomp)) + " ORFs. Saving to FASTA files...")
        # Open output files
        out_files[genome_name] = [outfile_dna + "_" + genome_name + DOT_FASTA, outfile_aa + "_" + genome_name + DOT_FASTA] # Store output file names     
        dna_file = open(out_files[genome_name][0], "w")
        aa_file = open(out_files[genome_name][1], "w")
        # Write input sequence ORFs
        for m in matches_input:
            print(m.span())
            write_orf(dna_file, aa_file, dna[m.regs[1][0]:m.regs[1][1]], m.regs[1][0] + 1, m.regs[1][1])
        # Write reverse complement ORFs
        for m in matches_revcomp:
            write_orf(dna_file, aa_file, revcomp[m.regs[1][0]:m.regs[1][1]], len(dna) - m.regs[1][0], len(dna) - m.regs[1][1] + 1)
        # Close output files
        dna_file.close()
        aa_file.close()
    print("==> DONE ("+str(i+1)+" genoms processed!)")
    return out_files

"""
TEST & MAIN
"""
# Test the 'findORFs' call with different ORF lengths
def test_findORFs(min, max):
    print("==> Testing ORF length " + str(min) + "-" + str(max) + "...")
    # FIND THE ORFs !!
    output_files = findORFs(IN_FNAME, OUT_DNA_FNAME, OUT_AA_FNAME, min, max)
    if len(output_files) > 0:
        print("\n==> Results:")
        print_output(output_files)
      
# MAIN CALL 
def main():
    test_findORFs(1,10)
    print("\n")
    test_findORFs(4,10)
main()