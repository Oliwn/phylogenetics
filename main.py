import Bio
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "be19b066@technikum-wien.at"
#handle = Entrez.einfo(db="gene")

#I got the following 5 lines from: https://stackoverflow.com/questions/13557981/getting-a-gene-sequence-from-entrez-using-biopython
GENE_NAME = "KCNT1"
handle = Entrez.esearch(db="nucleotide", term=GENE_NAME, retmax=100)
record = Entrez.read(handle)
handle.close()
ids = record["IdList"]
#testprint
print(ids)

orthologs = {}
for i in range(len(ids)):
    handle = Entrez.efetch(db="nucleotide", id=ids[i], rettype="gb", retmode="text")
    #print(handle.read())
    record = SeqIO.read(handle, "genbank")
    handle.close()
    key = record.annotations["organism"]
    if not key in orthologs:
        orthologs[key] = record
        print(key)
    if len(orthologs) == 5:
        break

print(orthologs)
print(len(orthologs))
#    record.
#    print(record.seq)

#handle = Entrez.efetch(db="nucleotide", id="AL158822", rettype="gb", retmode="text")

#records = SeqIO.read(handle, "gb")
#seq = records.seq
#for record in records:
 #   print(record)
