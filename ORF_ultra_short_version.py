from re import finditer, compile
IN_DNA, OUT_DNA, OUT_AA, IUPAC_FNAME = "input/in.fasta", "output/orfs_dna.fasta", "output/orfs_aa.fasta", "data/iupac.txt" # FILENAMES
IUPAC = dict(x.split(" ") for x in open(IUPAC_FNAME).read().splitlines()) # Prepares IUPAC dictionary
def iupac_transl(orf): return ''.join([IUPAC[orf[i:i+3]] for i in range(0, len(orf)-2, 3)]) # Translates an ORF using UIPAC dict
def orf_regex(min, max): return compile(r"(?=(ATG(?:...){" + str(min) + r"," + str(max) + r"}?(?<=TAA|TAG|TGA)))") # ORF RegEx pattern
def revcomp(dna): return dna[::-1].translate(str.maketrans("ATGC","TACG")) # Generates the reverse complement of a DNA string
def write_fasta_record(fasta_file, descr, data): [fasta_file.write(t) for t in (">" + descr + "\n", data + "\n")] # Writes to a file in FASTA format
def write_orf(dna_f, aa_f, orf, start, end): [write_fasta_record(f, "ORF_"+str(start)+"_"+str(end), d) for f,d in ((dna_f,orf),(aa_f, iupac_transl(orf)))] # Saves an ORF
def print_output(): [print(open(f).read()) for f in (OUT_DNA, OUT_AA)] # Print the output files
def findORFs(infile, outfile_dna, outfile_aa, minlen, maxlen): # FindORFs Function
    dna = open(infile).read().upper() # Open input DNA and translate to uppercase
    rev = revcomp(dna) # Get the reverse complement
    regex_pattern = orf_regex(minlen, maxlen) # Build the RegEx pattern
    matches_dna, matches_rev = tuple(finditer(regex_pattern, dna)), tuple(finditer(regex_pattern, rev)) # Find all matches 
    dna_f, aa_f = open(outfile_dna, "w"), open(outfile_aa, "w") # Open output files
    for m in matches_dna: write_orf(dna_f, aa_f, dna[m.regs[1][0]:m.regs[1][1]], m.regs[1][0] + 1,        m.regs[1][1]               ) # Write all input DNA matches
    for m in matches_rev: write_orf(dna_f, aa_f, rev[m.regs[1][0]:m.regs[1][1]], len(dna) - m.regs[1][0], len(dna) - m.regs[1][1] + 1) # Write all revcomp matches
    dna_f.close(), aa_f.close() # Close output files 
    return True # Success
[(print("Testing ORF length "+str(l)+"-"+str(h)), print_output() if findORFs(IN_DNA, OUT_DNA, OUT_AA, l, h) else None) for l,h in ((1,10), (4,10))] # Test with lengths 1-10 and 4-10 