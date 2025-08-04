## This script takes each read in the SAM file and tallies them to each TA site in the genome.
## This is counting (Tallying) the reads at each TA site.
## input:  python tally.py argv[1] argv[2] argv[3] > output.txt
## argv[1] is TA sites list for genome
## arvg[2] is mapped reads SAM format (single end mapping)
## argv[3] is gene list (no header)
## argv[4] is the genome fasta file (for contig lengths)
## The output file has 6 columns:  1) name of chromosome; 2) start position of TA site in genome; 3) end position of TA site in chromosome;
## 4) the locus in which each TA site sits; 5) The reads from insertions at that TA site on the Forward (+) strand; 6) reads from insertions at that TA site on the Reverse (-) strand

#python3 compatible, configured for multiple contigs7/22/25
from Bio import SeqIO  # make sure you have Biopython
import sys

chrITA={}

for line in open(sys.argv[1]):
    split= line.split()         #splits string into columns--treates consecutive spaces as 1 single tab
    if split[0]: #determine if col 1 is not empty
        chrom = split[0]
        TAsite=int(split[1]) #if yes, goes to col 2 and finds the TA position as defined by TAfinder.py
        chrITA[(chrom, TAsite)]=[0,0,'None'] #deposits all TA sites into a table, but with 0 reads

##counts reads at each TA site from SAM output
for line in open(sys.argv[2]):
    split=line.split('\t')
    if split[0][0] !="@": #ignores headers (analyzes rows don't have an @ symbol)
        contig = split[2]
        if split[1]=='16':  #if mapping is to the plus strand
            site=int(split[3])
            key = (contig, site)
            if key in chrITA:
                chrITA[key][0] += 1
        if split[1]=="0": #if mapping is on the minus strand
            site=int(split[3])-2 + len(split[9])
            key = (contig, site)
            if key in chrITA:
                chrITA[key][1] += 1  #if read has been found before, tally 1 more


# assign gene names
gene_dict = {}
prev_end = 0
IG = 'IG_1'
prev_chrom = 2

for line in open(sys.argv[3]):
    split = line.split('\t')
    #chrom = int(split[0]) #original code
    chrom = split[0] #changed to string
    gene = split[1]
    start = int(split[5])
    end = int(split[6])
    gene_dict[gene]=[start, end, chrom]
    IG = 'IG_' + gene
    if chrom != prev_chrom:
        #IG = 'IG_1' #I'm not sure why reset the name, since the intergenic region is before the gene it names, but leaving it in case it's needed for beginning of chromosome
        prev_end = 0
    if start > prev_end+1:
        gene_dict[IG]=[prev_end+1,start-1,chrom]
    prev_end = end
    prev_chrom = chrom
    
#This was hardcoded, need to define for each contig
# gene_dict['IG_chrm1end'] = [6537457,6537648,1]

# Get contig lengths from fasta using biopython
fasta_file = sys.argv[4]  # make sure this is the path to the .fasta file
contig_lengths = {record.id: len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# Add terminal IG region to end of each contig
# First, group all genes by contig
from collections import defaultdict
contig_ends = defaultdict(int)

for gene, (start, end, chrom) in gene_dict.items():
    if end > contig_ends[chrom]:
        contig_ends[chrom] = end

for chrom in contig_lengths:
    contig_len = contig_lengths[chrom]
    last_gene_end = contig_ends.get(chrom, 0)
    if last_gene_end < contig_len:
        IG_name = f'IG_{chrom}_end'
        gene_dict[IG_name] = [last_gene_end + 1, contig_len, chrom]

# store gene names into TA site dictionary
for k,v in gene_dict.items():
    start = v[0]
    end = v[1]
    chrom = v[2]
    #if chrom == 1: #why was this originally only for chromosome 1?
    for i in range(start,end+1):
        if (chrom, i) in chrITA: 
            #chrITA[i].append(k)
            chrITA[(chrom, i)][2] = k #assigns gene name to TA site
            
#print
CI_sites = list(chrITA.keys())
CI_sites.sort()

for i in CI_sites:
    print("%s\t%d\t%d\t%s\t%d\t%d" % (i[0],i[1],i[1]+1,chrITA[i][2],chrITA[i][0],chrITA[i][1])) #chromosome was hardcoded before to Chr I
    
