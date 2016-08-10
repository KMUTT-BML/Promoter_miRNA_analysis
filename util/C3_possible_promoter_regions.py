#!/usr/bin/python
import os
import sys
import re
from optparse import OptionParser
from seq_manage import Fasta_manager
from seq_manage import Gff_manager
from seq_manage import Genome_manager

if "--version" in sys.argv[1:]:
	# TODO - Capture version of Select_representative_miRNA
	print("Representative_miRNA v1.0")
	sys.exit(0)

# Parse Command Line
usage = """

Description:
This script designed for possible region of intergenic miRNA promoter
The information of intergenic miRNA that retrieving sequences related to each 
pre-miRNA downstream (x) and downstream (y) from the first nucleotide of pre-miRNA, 

with criteria
1) x = 200
2) 200 ≤ y ≤ 2,000 when sequence is not overlapped with protein-coding and miRNA gene
3) no N-gap in the sequence
(possible_promoter_iden.py)


Use as follows:
1) If you merge GFF before, Should type command following
$ python C3_proximal_promoter_iden.py \
	-i intergenic_miRNA.gff3 \
	--genome data/genome/Mesculenta_147_v4.1.fa \
	--gene_gff data/genome/Mesculenta_147_v4.1.gene.gff3 \
	-o out/step_3/possible_intergenic_miRNA_promoter.fa
"""

parser = OptionParser(usage=usage)
parser.add_option("-n", "--default", dest="default",
	default=None,
	help="Default input form thesis, please type 'thesis'")
parser.add_option("-i", "--intergenic_miRNA_gff", dest="miRNA_gff",
	default=None, metavar="FILE",
	help="GFF3 of miRNA, Gene feature file format (required)")
parser.add_option("-g", "--genome", dest="genome",
	default=None, metavar="FILE",
	help="Genome sequence in fasta format (required)")
parser.add_option("-m", "--gene_gff", dest="gene_gff",
	default=None, metavar="FILE",
	help="GFF3 of protein-coding genes, Gene feature file format (required)")
parser.add_option("-o", "--out_fa", dest="out_file",
	default=None, metavar="FILE",
	help="Sequence of intergenic promoter output file as fasta format (required)")
parser.add_option("-u", "--upstream", dest="upstream_len",
	default=2000, type="int",
	help="Length of upstream")
parser.add_option("-d", "--downstream", dest="downstream_len",
	default=200, type="int",
	help="Length of upstream")
parser.add_option("-l", "--min_len", dest="min_len",
	default=200, type="int",
	help="Minimum length of promoter")
parser.add_option("-f", "--Nfilter", dest="filterN",
	default='Y', type="str",
	help="Filtering cutting N out [Y/N](default: Y)")
options,args = parser.parse_args()

if options.default == 'thesis':
	miRNA_gff = "out/step_2/intergenic_miRNA.gff3"
	gene_gff = "data/genome/Mesculenta_147_v4.1.gene.gff3"
	genome = "data/genome/Mesculenta_147_v4.1.fa"
	promoter_miRNA_fa = "out/step_3/intergenic_miRNA_promoter.fa"
else:
	if not options.miRNA_gff:
		sys.exit("Missing miRNA GFF3 input file, -i ")
	if not os.path.isfile(options.miRNA_gff):
		sys.exit("No found miRNA Gff3: %r" % options.miRNA_gff)
	miRNA_gff = options.miRNA_gff
	if not options.genome:
		sys.exit("Missing genome sequence input file, -g ")
	if not os.path.isfile(options.genome):
		sys.exit("Not found genome sequence : %r" % options.genome)
	genome = options.genome
	if not options.gene_gff:
		sys.exit("Missing gene GFF  input file, -m ")
	if not os.path.isfile(options.gene_gff):
		sys.exit("Not found gene_gff sequence : %r" % options.gene_gff)
	gene_gff = options.gene_gff
	if not options.out_file:
		sys.exit("Missing output fasta file, -o ")
	out_file = open(options.out_file, 'w')
	out_file_pos = open(options.out_file + ".pos", 'w')


intergenic_miRNAs = Gff_manager(miRNA_gff, "miRNA_primary_transcript")
gene = Gff_manager(gene_gff, "gene")
genome = Fasta_manager(genome)


def filterUpstream(seq, min_len=400):
	foundN = seq.rfind('N')
	if len(seq)-foundN >= min_len:
		seq = seq[foundN + 1:]
	else:
		seq = ''
	return seq

def filterUpstreamHaveN(seq, min_len=100):
#Finding N from start sequence
	seq = seq.upper()
	if seq.find('N')==0:
		nu_start = re.search('[ATGC]+',seq).start()
		seq = seq[nu_start:]
	if len(seq)>min_len:
		return seq
	else:
		return ''


count_minus=0
count_plus=0

# Classification miRNA into 3 types and count their
for miRNA in intergenic_miRNAs.getTable():
	if(miRNA[2] == 'miRNA_primary_transcript'):
		miRName = miRNA[8]['Name']
		miRID = miRNA[8]['ID']
		stand = miRNA[6]
		chromosome = miRNA[0]
		start = miRNA[3]
		stop = miRNA[4]
		ref = miRNA[1]
		list = gene.isIntergenic(chromosome, start, stop, stand)
		list_miRNAs = intergenic_miRNAs.isIntergenic(chromosome, start, stop, stand)
		
		if(stand == '+'):
			forward_gene_stop = 1
		else:
			forward_gene_stop = genome.getChromosomeLength(chromosome)
		if(list[0] == ''):
			type = "Intergenic"
			if not list[1]:
				forward_gene = ""
				if list_miRNAs[1]:
					forward_gene = list_miRNAs[1][1][3:]
					forward_gene_stop = list_miRNAs[1][2]
			elif list[1] and not list_miRNAs[1]:
				forward_gene = list[1][1][3:-5]
				forward_gene_stop = list[1][2]
			elif list[1][2] >= list_miRNAs[1][2]:
				forward_gene = list[1][1][3:-5]
				forward_gene_stop = list[1][2]
			elif list[1][2] < list_miRNAs[1][2]:
				forward_gene = list_miRNAs[1][1][3:]
				forward_gene_stop = list_miRNAs[1][2]
				
			if not list[2]:
				backward_gene = ""
			else:
				backward_gene = list[2][1][:-5]
		else:
			type = list[0][0][3]
			forward_gene = list[0][0][0][:-5]
			backward_gene = ""

		# print(miRID, ref, type)
		upstrem_len = options.upstream_len
		downstream_len = options.downstream_len
		sequence = ''
		filterSeq = ''
		out_file_pos.write("\t".join([miRID, miRName, type, chromosome, stand, str(start), str(stop), forward_gene, str(forward_gene_stop)]) + "\t")
		if(stand == '+'):
			if(start + downstream_len - 1 <= genome.getChromosomeLength(chromosome)):
				if(start - forward_gene_stop > upstrem_len):
					sequence = genome.getSequence(chromosome, start - upstrem_len, start + downstream_len - 1, '+')
					out_file_pos.write("\t".join([str(start - upstrem_len), str(start + downstream_len - 1), str(len(sequence))])+ "\t")
				else:
					sequence = genome.getSequence(chromosome, forward_gene_stop + 1, start + downstream_len - 1, '+')
					out_file_pos.write("\t".join([str(forward_gene_stop + 1), str(start + downstream_len - 1), str(len(sequence))])+ "\t")
			# else:
				# print("Can't retrieve sequence because sequence length in ", chromosome, "downstream have length", genome.getChromosomeLength(chromosome), "less than", upstrem_len)	 
		else:
			if(stop - downstream_len + 1 >= 1):
				if(forward_gene_stop - stop > upstrem_len):
					sequence = genome.getSequence(chromosome, stop - downstream_len + 1, stop + upstrem_len, '-')
					out_file_pos.write("\t".join([str(stop - downstream_len + 1), str(stop + upstrem_len), str(len(sequence))]) + "\t")

				else:
					sequence = genome.getSequence(chromosome, stop - downstream_len + 1, forward_gene_stop - 1, '-')
					out_file_pos.write("\t".join([str(stop - downstream_len + 1), str(forward_gene_stop - 1), str(len(sequence))]) + "\t")
			# else:
				# print("Can't retrieve sequence because sequence length in downstream less than chromosome", chromosome)
		
		if options.filterN == 'Y':
			filterSeq = filterUpstream(sequence, options.min_len)
		else:
			filterSeq = filterUpstreamHaveN(sequence, options.min_len)

		if(len(filterSeq) > options.min_len):
			out_file.write(">" + "|".join([miRID, stand, chromosome, str(start), str(stop)]) + "|")
			if(stand == '+'):
				out_file.write("|".join([str(start + downstream_len - len(filterSeq)), str(start + downstream_len - 1), str(len(filterSeq))])+ "\n")
				out_file_pos.write("\t".join([str(start + downstream_len - len(filterSeq)), str(start + downstream_len - 1), str(len(filterSeq))])+ "\n")
				count_plus +=1
			elif(stand == '-'):
				out_file.write("|".join([str(stop - downstream_len + 1), str(stop - downstream_len + len(filterSeq)), str(len(filterSeq))])+ "\n")
				out_file_pos.write("\t".join([str(stop - downstream_len + 1), str(stop - downstream_len + len(filterSeq)), str(len(filterSeq))])+ "\n")
				count_minus +=1
			out_file.write(filterSeq + "\n")
		else:
			out_file_pos.write('\n')
print("Number of proximal intergenic miRNA promoter regions		    : ", count_plus + count_minus)
print("Number of proximal intergenic miRNA promoter regions in plus : ", count_plus)
print("Number of proximal intergenic miRNA promoter regions in minus: ", count_minus)

# print(countDownstreamErr, "can't retrived upstream becuase downstream less than", upstrem_len)
# print("Pass upstream filtered ", countFilteredSeq_plus + countFilteredSeq_minus,
# 	"[", countFilteredSeq_plus, "plus stand|", countFilteredSeq_minus, "minus stand]")
