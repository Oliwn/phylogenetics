import Bio
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import Seq
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

Entrez.email = "be19b066@technikum-wien.at"
#handle = Entrez.einfo(db="gene")

#I got the following 5 lines from: https://stackoverflow.com/questions/13557981/getting-a-gene-sequence-from-entrez-using-biopython
GENE_NAME = "KCNT1"
handle = Entrez.esearch(db="nucleotide", term=GENE_NAME, retmax=100)
records = Entrez.read(handle)
handle.close()
ids = records["IdList"]
#testprint
print(ids)

orthologs = []
organisms = []
for i in range(len(ids)):
    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", retmode="text")
    #print(handle.read())
    records = SeqIO.read(handle, "genbank")
    handle.close()
    organism = records.annotations["organism"]
    if not organism in organisms:
        orthologs.append(records)
        organisms.append(organism)
        print(organism)
    if len(orthologs) == 5:
        break


print(orthologs)



GCs = []
for rec in orthologs:
    GCs.append(GC(rec.seq))
print(GCs)

#Since the sequences are not the same length, the shorter ones get padded with dots
#Source: https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
maxlen = max(len(rec.seq) for rec in orthologs)
for rec in orthologs:
    rec.id = rec.annotations["organism"]
    if len(rec.seq) != maxlen:
        sequence = str(rec.seq).ljust(maxlen, ".")
        rec.seq = Seq.Seq(sequence)
print(orthologs)
align = MultipleSeqAlignment(orthologs)

#AlignIO.write(align, "alignments.phy", "clustal")

#aligner = PairwiseAligner()
#alignments = []
# for i in range(len(orthologs)):
#     j = i + 1
#     while j <= 4:
#         print(str(i) + " " + str(j))
#         alignment = aligner.align(orthologs[i].seq, orthologs[j].seq)
#         AlignIO.write(alignment, "alignments.phy", "phylip")
#         alignments.append(alignment[0])
#         j += 1


#with open("alignments.phy", "r") as aln:
#    alignment = AlignIO.read(aln, "clustal")
#print(type(alignment))

#I used the following video and the tutorial which is linked in the description for creating the tree: https://www.youtube.com/watch?v=wBdz3vFQ4Ks

calculator = DistanceCalculator("identity")
distanceMatrix = calculator.get_distance(align)
#print(distanceMatrix)

constructor = DistanceTreeConstructor(calculator)
tree = constructor.build_tree(align)
#print(tree)

fig = Phylo.draw(tree)

proteins = []
for rec in orthologs:
    proteins.append(rec.seq.transcribe().translate(to_stop=True))
print("--------------------")
for pro in proteins:
    print(str(pro))
    temp = ProteinAnalysis(str(pro))
    print(temp.aromaticity())
    print(temp.instability_index())
    print(temp.isoelectric_point())
    print(temp.secondary_structure_fraction())


#calculator = DistanceCalculator("identity")
#distanceMatrix = calculator.get_distance(alignments)
#print(distanceMatrix)
#alignments = pairwise2.align.globalxx(list(orthologs.values())[0].seq, list(orthologs.values())[1].seq)
#print(format_alignment(*alignments[0]))