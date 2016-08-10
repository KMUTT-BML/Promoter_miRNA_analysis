import os
import sys
from seq_manage import Fasta_manager
from seq_manage import Gff_manager
from seq_manage import Genome_manager
import re

#!/usr/bin/python
from optparse import OptionParser

if "--version" in sys.argv[1:]:
	# TODO - Capture version of Select_representative_miRNA
	print("Representative_miRNA v1.0")
	sys.exit(0)

# Parse Command Line
usage = """

Description:
This script designed for classcification of pre-miRNA into 3 types based on genome location
1. Intergenic miRNA: miRNA located between genes
2. Intronic miRNA: miRNA located within gene not overlaping to coding sequence (CDS)
3. Exonic miRNA: miRNA located overlaping coding sequence (CDS)

Use as follows:
1) If you merge GFF before, Should type command following
$ python classify_miRNA.py -i rep_miRNA_out.gff3 -o classified
"""

parser = OptionParser(usage=usage)
parser.add_option("-d", "--default", dest="default",
	default=None,
	help="Default input form thesis, please type 'thesis'")
parser.add_option("-i", "--miRNA_gff", dest="miRNA_gff",
	default=None, metavar="FILE",
	help="GFF3 of miRNA, Gene feature file format (required)")
parser.add_option("-g", "--genome", dest="genome",
	default=None, metavar="FILE",
	help="Genome sequence in fasta format (required)")
parser.add_option("-m", "--gene_gff", dest="gene_gff",
	default=None, metavar="FILE",
	help="GFF3 of protein-coding genes, Gene feature file format (required)")
options,args = parser.parse_args()

if options.default == 'thesis':
	miRNA_gff = "cassava_genome/Mes_pri_miRNA_loci.gff3"
	gene_gff = "cassava_genome/Mesculenta_147_v4.1.gene.gff3"
	genome = "cassava_genome/Mesculenta_147_v4.1.fa"
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

miRNAs = Gff_manager(miRNA_gff, "miRNA_primary_transcript")
gene = Gff_manager(gene_gff, "gene")
genome = Fasta_manager(genome)

os.makedirs(os.path.dirname('out_classify'), exist_ok=True)
out_intergenic_mi_file_name = "out_classify/intergenic_miRNA.gff3"
out_exonic_mi_file_name = "out_classify/exonic_miRNA.gff3"
out_intronic_mi_file_name = "out_classify/intronic_miRNA_.gff3"

out_file_inter = open(out_intergenic_mi_file_name, 'w')
out_file_exo = open(out_exonic_mi_file_name, 'w')
out_file_intro = open(out_intronic_mi_file_name, 'w')

def filterUpstreamHaveN(seq, limited=100):
    #Finding N from start sequence
    seq = seq.upper()
    if seq.find('N')==0:
        nu_start = re.search('[ATGC]+',seq).start()
        seq = seq[nu_start:]
    if len(seq)>limited:
        return seq
    else:
        return ''

# Classification miRNA into 3 types and count their
countIntergenic = countIntronic = countExonic = 0
countDownstreamErr = 0
countIntergenic_plus = countIntergenic_minus = 0
countIntronic_plus = countIntronic_minus = 0
countExonic_plus = countExonic_minus = 0
countFilteredSeq_plus = countFilteredSeq_minus = 0

for miRNA in miRNAs.getTable():
    if(miRNA[2] == 'miRNA_primary_transcript'):
        miRName = miRNA[8]['Name']
        miRID = miRNA[8]['ID']
        stand = miRNA[6]
        chromosome = miRNA[0]
        start = miRNA[3]
        stop = miRNA[4]
        ref = miRNA[1]
        list = gene.isIntergenic(chromosome, start, stop, stand)

        list_miRNAs = miRNAs.isIntergenic(chromosome, start, stop, stand)
        
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
                backward_gene = list[2][1][3:-5]
        else:
            type = list[0][0][3]
            forward_gene = list[0][0][0][3:-5]
            backward_gene = ""
       
        # print(chromosome, ref, type + "_miRNA", start, stop, ".", stand, ".", "miRNA_id=" + miRID + ";Name=" + miRName + ";type=" + type 
        #         + ";font_gene="+ forward_gene + ";back_gene="+backward_gene, genome.getSequence(chromosome, start, stop, stand) , sep="\t")
   
    
        if(type == 'Intergenic'):
            countIntergenic += 1
            # print(miRID, ref, type)
            upstrem_len = 1500
            downstream_len = 0
            sequence = ''
            filterSeq = ''
            print(miRID, miRName, type, chromosome, stand, start, stop, forward_gene, forward_gene_stop, sep=',', end=',')
            if(stand == '+'):
                countIntergenic_plus += 1
                if(start + downstream_len - 1 <= genome.getChromosomeLength(chromosome)):
                    if(start - forward_gene_stop > upstrem_len):
                        sequence = genome.getSequence(chromosome, start - upstrem_len, start + downstream_len - 1, '+')
                        print(start - upstrem_len, start + downstream_len - 1, len(sequence), sep=',', end=',')
                    else:
                        sequence = genome.getSequence(chromosome, forward_gene_stop + 1, start + downstream_len - 1, '+')
                        print(forward_gene_stop + 1, start + downstream_len - 1, len(sequence), sep=',' , end=',')
                else:
                    # print("Can't retrieve sequence because sequence length in ", chromosome, "downstream have length", genome.getChromosomeLength(chromosome), "less than", upstrem_len)
                    countDownstreamErr += 1
                      
            else:
                countIntergenic_minus += 1
                if(stop - downstream_len + 1 >= 1):
                    if(forward_gene_stop - stop > upstrem_len):
                        sequence = genome.getSequence(chromosome, stop - downstream_len + 1, stop + upstrem_len, '-')
                        print((stop - downstream_len + 1), str(stop + upstrem_len), len(sequence), sep=',', end=',')
                    else:
                        sequence = genome.getSequence(chromosome, stop - downstream_len + 1, forward_gene_stop - 1, '-')
                        print((stop - downstream_len + 1), (forward_gene_stop - 1), len(sequence), sep=',', end=',')
                else:
                    # print("Can't retrieve sequence because sequence length in downstream less than chromosome", chromosome)
                    countDownstreamErr += 1
        
            filterSeq = filterUpstreamHaveN(sequence, 50)
            # filterSeq = filterUpstream(sequence, downstream_len + 2000)
            if(stand == '+' and filterSeq != ''):
                countFilteredSeq_plus += 1
                print(start + downstream_len - len(filterSeq), start + downstream_len - 1, len(filterSeq), filterSeq, sep=',')
            elif(stand == '-' and filterSeq != ''):
                countFilteredSeq_minus += 1
                print(stop - downstream_len + 1, stop - downstream_len + len(filterSeq), len(filterSeq), filterSeq, sep=',')
            else:
                print()
        elif(type == 'Intronic'):
            countIntronic += 1
            if stand == '+':
                countIntronic_plus+=1
            else:
                countIntronic_minus+=1
            # print(miRID, ref, type)
            # print(miRID, miRName, type, chromosome, stand, start, stop, forward_gene, '-', sep=',')
            # print(chromosome, ref, type + "_miRNA", start, stop, ".", stand, ".", "miRNA_id=" + miRID + ";Name=" + miRName + ";type=" + type 
            #         + ";font_gene="+ forward_gene + ";back_gene=", genome.getSequence(chromosome, start, stop, stand) , sep="\t")
         
         
         
        elif(type == 'Exonic'):
            countExonic += 1
            if stand == '+':
                countExonic_plus+=1
            else:
                countExonic_minus+=1
            # print(miRID, ref, type)
            # print(chromosome, ref, type + "_miRNA", start, stop, ".", stand, ".", "miRNA_id=" + miRID + ";Name=" + miRName + ";type=" + type 
            #         + ";font_gene="+ forward_gene + ";back_gene=", genome.getSequence(chromosome, start, stop, stand) , sep="\t")
   
 

print(countDownstreamErr, "can't retrived upstream becuase downstream less than", upstrem_len)
print("Pass upstream filtered ", countFilteredSeq_plus + countFilteredSeq_minus,
      "[", countFilteredSeq_plus, "plus stand|", countFilteredSeq_minus, "minus stand]")

print("count intergenic:", countIntergenic,
      "[", countIntergenic_plus, "plus stand|", countIntergenic_minus, "minus stand]")
print("count introgenic:", countIntronic,
      "[", countIntronic_plus, "plus stand|", countIntronic_minus, "minus stand]")
print("count exonic:", countExonic,
      "[", countExonic_plus, "plus stand|", countExonic_minus, "minus stand]")
