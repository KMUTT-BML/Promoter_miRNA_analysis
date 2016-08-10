# #!/usr/bin/python
# import os
# import sys
# import re
# from optparse import OptionParser

# from seq_manage import Fasta_manager
# from seq_manage import Gff_manager
# from seq_manage import Genome_manager

# if "--version" in sys.argv[1:]:
# 	# TODO - Capture version of predict proximal promoter region and retriving promoter sequence.
# 	print("Representative_miRNA v1.0")
# 	sys.exit(0)

# # Parse Command Line
# usage = """

# Description:
# This script designed for retrived possible region of intergenic miRNA and retrived promoter sequence  

# with criteria
# 1) upstream 1.5 kb not overlap with other gene
# 2) 200 ≤ y ≤ 2,000 when sequence is not overlapped with protein-coding and miRNA gene
# 3) no N-gap in the sequence
# (possible_promoter_iden.py)


# Use as follows:
# 1) If you merge GFF before, Should type command following
# $ python C3_proximal_promoter_iden.py \
# 	-i intergenic_miRNA.gff3 \
# 	--genome data/genome/Mesculenta_147_v4.1.fa \
# 	--gene_gff data/genome/Mesculenta_147_v4.1.gene.gff3 \
# 	-o out/step_3/possible_intergenic_miRNA_promoter.fa
# """

# parser = OptionParser(usage=usage)
# parser.add_option("-d", "--default", dest="default",
# 	default=None,
# 	help="Default input form thesis, please type 'thesis'")
# parser.add_option("-i", "--intergenic_miRNA_gff", dest="miRNA_gff",
# 	default=None, metavar="FILE",
# 	help="GFF3 of miRNA, Gene feature file format (required)")
# parser.add_option("-g", "--genome", dest="genome",
# 	default=None, metavar="FILE",
# 	help="Genome sequence in fasta format (required)")
# parser.add_option("-m", "--gene_gff", dest="gene_gff",
# 	default=None, metavar="FILE",
# 	help="GFF3 of protein-coding genes, Gene feature file format (required)")
# parser.add_option("-o", "--out_fa", dest="out_file",
# 	default=None, metavar="FILE",
# 	help="Sequence of intergenic promoter output file as fasta format (required)")
# options,args = parser.parse_args()

# if options.default == 'thesis':
# 	miRNA_gff = "out/step_2/intergenic_miRNA.gff3"
# 	gene_gff = "data/genome/Mesculenta_147_v4.1.gene.gff3"
# 	genome = "data/genome/Mesculenta_147_v4.1.fa"
# 	promoter_miRNA_fa = "out/step_3/intergenic_miRNA_promoter.fa"
# else:
# 	if not options.miRNA_gff:
# 		sys.exit("Missing miRNA GFF3 input file, -i ")
# 	if not os.path.isfile(options.miRNA_gff):
# 		sys.exit("No found miRNA Gff3: %r" % options.miRNA_gff)
# 	miRNA_gff = options.miRNA_gff
# 	if not options.genome:
# 		sys.exit("Missing genome sequence input file, -g ")
# 	if not os.path.isfile(options.genome):
# 		sys.exit("Not found genome sequence : %r" % options.genome)
# 	genome = options.genome
# 	if not options.gene_gff:
# 		sys.exit("Missing gene GFF  input file, -m ")
# 	if not os.path.isfile(options.gene_gff):
# 		sys.exit("Not found gene_gff sequence : %r" % options.gene_gff)
# 	gene_gff = options.gene_gff
# 	if not options.out_file:
# 		sys.exit("Missing output fasta file, -o ")
# 	out_file = open(options.out_file, 'w')
# 	out_file_pos = open(options.out_file + ".pos", 'w')


# intergenic_miRNAs = Gff_manager(miRNA_gff, "miRNA_primary_transcript")
# gene = Gff_manager(gene_gff, "gene")
# genome = Fasta_manager(genome)

file_PPde_name = 'out/step_3/possible_intergenic_miRNA_promoter_PPde.txt'

file = open(file_PPde_name, 'r')
writeFile = open('out/step_3/primary_miRNA_loci.gff3','w')
data = file.read().split('\nID\t')

data = data[1:]
corePromoter = []
list_found_promoter = []
downstream_len = 200
numberOfCorePromoterInPlus = 0
numberOfCorePromoterInMinus = 0
numberOfCorePromoterPredicted = 0
print("pre_name","stand", "numberOfPromoter","prom_len","distant_from_pre","Not_overlap", sep="\t")
for i in data:
	i = i.split('\n>')
	#name = [id, stand, chromosome, pre_start, pre_end, up_start, up_end, up_length]
	name = i[0][:i[0].find('\t')].split('|')
	# print(name)
	# exit()
	#change string to intiger
	for j in range(3,len(name)):
		name[j] = int(name[j])
	if(len(i)>1):
		last = len(i)-1
		numberOfCorePromoterPredicted +=1
		print(name[0], name[1], len(i)-1, sep='\t', end="\t")
	else:
		last = 1
		print(name[0], name[1], len(i)-1, 0, sep='\t')
	for seq in i[last:]:
		prom = seq.split('\t')
		prom[0] = prom[0].split('..')
		prom[0][0] = int(prom[0][0])
		prom[0][1] = int(prom[0][1])
		if(name[1]=='+'):
			GenomePos = [name[5] + prom[0][0] - 1, name[5] + int(prom[0][1]) - 1]
			prom.insert(0,GenomePos)
			# print("\n",name)
			# print(prom)
			# distant from precursor microRNA
			print(prom[2], name[7]-prom[1][1]-downstream_len, sep = "\t", end = "\t")
			print(name[3]-prom[0][1]>=1)
			writeFile.write(name[2]+"\tmiRNA\tmiRNA_primary_transcript\t"+str(prom[0][1]-1)+"\t"+str(name[4])+'\t.\t' +name[1]+ '\t.\t' +
				'ID='+name[0]+'_pri'+';Alias=;Name='+name[0]+'_pri\n')
			# writeFile.write(name[2]+"\tmiRNA\tmiRNA_precursor\t"+ str(name[3]) +"\t"+str(name[4]) + '\t.\t' +name[1] +'\t.\t' + 
			# 	'ID='+name[0]+'_pre'+';Parent='+ name[0] +'_pri;Name='+name[0]+'_pre\n')
			list_found_promoter.append(name[0])
			numberOfCorePromoterInPlus +=1
		else:
			GenomePos = [name[6] - int(prom[0][1])+1, name[6] - int(prom[0][0]) + 1]
			prom.insert(0,GenomePos)
			print(prom[2],name[7]-prom[1][1]-downstream_len, sep = "\t", end="\t")
			print(prom[2],name[7]-prom[1][1]-downstream_len, sep = "\t", end="\t")
			print(prom[0][0]-name[5]>=1)
			# print(name)
			# print(prom)
			writeFile.write(name[2] +"\tmiRNA\tmiRNA_primary_transcript\t" + str(name[3]) +"\t"+str(prom[0][0]-1) +'\t.\t'+ name[1] + '\t.\t' + 
				'ID='+name[0]+'_pri'+';Alias=;Name='+name[0]+'_pri\n')
			# writeFile.write(name[2]+"\tmiRNA\tmiRNA_precursor\t"+ str(name[3]) +"\t"+str(name[4]) + '\t.\t' +name[1] +'\t.\t' + 
			# 	'ID='+name[0]+'_pre'+';Parent='+ name[0] +'_pri;Name='+name[0]+'_pre\n')
			list_found_promoter.append(name[0])
			numberOfCorePromoterInMinus +=1
			corePromoter.append(name[0])

file_miRNA_gff_all = 'out/step_2/intergenic_miRNA.gff3'
file_miRNA_all = open(file_miRNA_gff_all,'r')
all_miRNA = file_miRNA_all.read().split('\n')
count = 0
for i in all_miRNA[:-1]:
	id = i[i.find('ID=')+3:i.find(';',i.find('ID=')+3)]
	if( id not in list_found_promoter):
		count+=1
		writeFile.write(i.replace(';Alias=;','_pri;Alias=;').replace('total', 'miRNA')+'_pri\n')
		# writeFile.write(i.replace(';Alias=;','_pre;Parent='+id+'_pri;').replace('total', 'miRNA').replace('miRNA_primary_transcript','miRNA_precursor')+'_pre\n')

print("Number of core promoter in plus stand: ", numberOfCorePromoterInPlus)
print("Number of core promoter in minus stand: ", numberOfCorePromoterInMinus)
print("Number of core promoter:", numberOfCorePromoterPredicted)
print("Additional miRNA gene not have promoter:", count)