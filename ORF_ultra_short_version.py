from re import finditer, compile
from fasta_utils import fasta_to_dict, write_fasta_record
IN_DNA, OUT_DNA, OUT_AA, IUPAC_FNAME = "input/in.fasta", "output/orfs_dna", "output/orfs_aa", "data/iupac.txt" # FILENAMES
IUPAC = dict(x.split(" ") for x in open(IUPAC_FNAME).read().splitlines()) # Prepares IUPAC dictionary
def iupac_transl(orf): return ''.join([IUPAC[orf[i:i+3]] for i in range(0, len(orf)-2, 3)]) # Translates an ORF using UIPAC dict
def orf_regex(min, max): return compile(r"(?=(ATG([ACTG]{3}(?<!ATG|TAA|TAG|TGA)){"+str(min-1)+r","+str(max-1)+r"}(TAA|TAG|TGA)))") # ORF RegEx pattern with min and max length
def revcomp(dna): return dna[::-1].translate(str.maketrans("ATGC","TACG")) # Generates the reverse complement of a DNA string
def write_orf(dna_f, aa_f, orf, start, end): [write_fasta_record(f, "ORF_"+str(start)+"_"+str(end), d) for f,d in ((dna_f,orf),(aa_f, iupac_transl(orf)))] # Saves an ORF
def print_output(output_files): [[print(l) for l in (f, open(f).read())] for f in [f for gen in output_files.values() for f in gen]] # Print the output files
def findORFs(infile, outfile_dna, outfile_aa, minlen, maxlen): # FindORFs Function
    inp_genoms = fasta_to_dict(infile) # Open input file
    out_files = {} # Variable to store output file names
    regex_pattern = orf_regex(minlen, maxlen) # Build the RegEx pattern
    for gen_id, dna in inp_genoms.items(): # Loop over all genoms 
        rev = revcomp(dna) # Get the reverse complement
        matches_dna, matches_rev = [tuple(finditer(regex_pattern, g)) for g in (dna, rev)] # Find all matches 
        out_files[gen_id] = [f + "_" + gen_id + ".fasta" for f in (outfile_dna,outfile_aa)]
        dna_f, aa_f = [open(out_files[gen_id][i], "w") for i in (0,1)] # Open output files
        [[[write_orf(dna_f, aa_f,g,s,e) for g,s,e in [(gen[x[0]:x[1]], x[0]+1, x[1]) if r else (gen[x[0]:x[1]], len(gen)-x[0], len(gen)-x[1]+1) for x in [m.regs[1] for m in ms]]]] for ms, gen, r in [(matches_dna, dna, False),(matches_rev, rev, True)]] # Write ORFs
        [f.close() for f in (dna_f, aa_f)] # Close files
    return out_files # Success
[(print("Testing ORF length "+str(l)+"-"+str(h)), print_output(findORFs(IN_DNA, OUT_DNA, OUT_AA, l, h))) for l,h in ((1,10), (4,10))] # Test with lengths 1-10 and 4-10 