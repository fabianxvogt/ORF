from fnmatch import translate
import re
import pathlib

min = "1"
max = "10"
genome = "CATGTTATGCTAGCCACATAA"

# ORF DNA 
pattern = re.compile(r"(?=(ATG(?:...){" + min + r"," + max+ r"}?(?<=TAA|TAG|TGA)))")
def revcomp(dna):
    return dna[::-1].translate(str.maketrans("ATGC","TACG"))
def orfs(dna):
    return set(pattern.findall(dna) + pattern.findall(revcomp(dna)))

orfs = orfs(genome)
print(orfs)

# Load codesun (Amino Acid translation table)
codesun = {}
with open('data/codesun.txt') as f:
    codesun = dict(x.split(" ") for x in f.read().splitlines())

for orf in orfs:
    # translate amino acids
    translation = ""
    for i in range(0, len(orf)-2, 3):
        translation += codesun[orf[i:i+3]]
    pos = genome.find(orf) 
    if (pos == -1): 
        end = revcomp(genome).find(orf) + 3
        start = end + len(orf) - 1
    else:
        start = pos + 1
        end = pos + len(orf)
    
    print(orf)
    print(translation)
    print(str(start) + "_" + str(end))