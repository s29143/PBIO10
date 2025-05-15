from Bio import Entrez,SeqIO
import matplotlib.pyplot as plt
import csv
Entrez.email=input("mail")
Entrez.api_key=input("APIkey")
t=input("taxid")
m=int(input("Minseq"))
l=int(input("Maxseq"))
h=Entrez.esearch(db="nucleotide",term=f"txid{t}[Organism]",usehistory="y",retmax=10000)
s=Entrez.read(h)
h=Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=0,retmax=min(5, 500),webenv=s["WebEnv"],query_key=s["QueryKey"])
w=[r for r in list(SeqIO.parse(h,"genbank"))if m<=len(r.seq)<=l]
with open(f"tid_{t}.csv","w",newline="") as c:
    writer=csv.writer(c)
    writer.writerow(["Accession","Length","Description"])
    for r in w:
        writer.writerow([r.id,len(r.seq),r.description])
w=sorted(w,key=lambda x:len(x.seq),reverse=True)
plt.plot([r.id for r in w],[len(r.seq)for r in w],marker='o')
plt.savefig(f"tid_{t}.png")