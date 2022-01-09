#Link to github-repository: https://github.com/Oliwn/phylogenetics.git
from Bio import Entrez, SeqIO, Phylo, Seq
from Bio.SeqUtils import GC
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#If not stated otherwise, I found most of the commands in this "documentation": https://biopython.org/docs/1.76/api/Bio.html
#I worked on this program together with Theresa Plumer and Gralf Antiga in a call and we helped each other out when needed,
#so our code will show some similarities.

Entrez.email = "be19b066@technikum-wien.at"

#I got the following 5 lines from: https://stackoverflow.com/questions/13557981/getting-a-gene-sequence-from-entrez-using-biopython
#It searches for the gene name and saves the found IDs for further search
GENE_NAME = "KCNT1"
handle = Entrez.esearch(db="nucleotide", term=GENE_NAME, retmax=100)
records = Entrez.read(handle)
handle.close()
ids = records["IdList"]

#Searching the DB for the found IDs. The sequence gets saved if the organism of the sequence has not already been found
print("Searching for sequences...")
orthologs = []
organisms = []
for i in range(len(ids)):
    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", retmode="text")
    records = SeqIO.read(handle, "genbank")
    handle.close()
    organism = records.annotations["organism"]
    if not organism in organisms:
        orthologs.append(records)
        organisms.append(organism)
        print("Found a sequence in organism: " + organism)
    if len(orthologs) == 5:
        break

GCs = []
for rec in orthologs:
    GCs.append(GC(rec.seq))
print("Calculated GCs")

#Since the sequences are not the same length, the shorter ones get padded with dots
#Source: https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
maxlen = max(len(rec.seq) for rec in orthologs)
for rec in orthologs:
    rec.id = rec.annotations["organism"]
    if len(rec.seq) != maxlen:
        sequence = str(rec.seq).ljust(maxlen, ".")
        rec.seq = Seq.Seq(sequence)
align = MultipleSeqAlignment(orthologs)
print("Sequences aligned")

#I used the following video and the tutorial which is linked in the description for creating the tree: https://www.youtube.com/watch?v=wBdz3vFQ4Ks
calculator = DistanceCalculator("identity")
distanceMatrix = calculator.get_distance(align)
constructor = DistanceTreeConstructor(calculator)
tree = constructor.build_tree(align)
fig = Phylo.draw(tree)

#The nucleotide sequences are translated and saved. Other parameters are calculated and partially saved.
proteins = []
for rec in orthologs:
    proteins.append(rec.seq.translate(to_stop=True))
print("Translating the nucletoide sequences in protein ones...")
arom = []
instab = []
for pro in proteins:
    print("Seq: " + str(pro))
    temp = ProteinAnalysis(str(pro))
    arom.append(temp.aromaticity())
    instab.append(temp.instability_index())
    print("Aromaticity: " + str(temp.aromaticity()))
    print("Instability: " + str(temp.instability_index()))
    print("Isoelectric point: " + str(temp.isoelectric_point()))
    print("Secondary structure fraction: " + str(temp.secondary_structure_fraction()) + "\n")

#The required parameters are written into a csv-file
with open("orthologs.csv", "w") as f:
    f.write("Accession number,Organism,Length,GC%,Instability,Aromaticity\n")
    for i in range(len(orthologs)):
        f.write(orthologs[i].name + "," + orthologs[i].annotations["organism"] + "," + str(len(orthologs[i].seq)) + "," + str(GCs[i]) + "," + str(instab[i]) + "," + str(arom[i]) + "\n")

print("\nSaved information in 'orthologs.csv'")