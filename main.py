import Bio
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "be19b066@technikum-wien.at"
#handle = Entrez.einfo(db="gene")
handle = Entrez.efetch(db="nucleotide", id="AL158822", rettype="gb", retmode="text")

records = SeqIO.read(handle, "gb")
seq = records.seq
#for record in records:
 #   print(record)
