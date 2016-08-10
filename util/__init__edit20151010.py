import re
import gzip
import codecs
from operator import itemgetter, attrgetter

utf8 = codecs.getreader('UTF-8')

class Fasta_manager(object):
	def __init__(self, fastaFile):
		self.chromosomeLength = {}
		self.chromosomeSeq = {}
		self.chromosomeStatistics = {}	# Length, GC, AT, N
		sumGC = sumAT = sumN = sumLength = 0
		
		if(fastaFile.find('.gz') > 0):
			filegz = gzip.open(fastaFile, 'rb')
			self.file = utf8(filegz)
		else:
			self.file = open(fastaFile, 'r')
		fasta = self.file.read().split('>')
		fasta = fasta[1:]
		for chromosome in fasta:
			if (chromosome[:50].find(' ') < 0):
				header = chromosome[:chromosome[:50].find('\n')]
			else:
				header = chromosome[:chromosome[:50].find(' ')]
			sequence = chromosome[chromosome.find('\n'):-1].replace('\n', '')
			
			length = len(sequence)
			self.chromosomeSeq[header] = sequence
			self.chromosomeLength[header] = length

	
	
	def getGCcontent(self, sequence):
		GC = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
		AT = sequence.count('A') + sequence.count('T') + sequence.count('a') + sequence.count('t')
		return float(GC) * 100 / (AT + GC)
	def getGC(self, sequence):
		return sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
	def getStatisticSequence(self, sequence):
		GC = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
		AT = sequence.count('A') + sequence.count('T') + sequence.count('a') + sequence.count('t')
		N = sequence.count('N') + sequence.count('n')
		return [len(sequence), GC, AT, N, float(GC) * 100 / (AT + GC)]
	def getStatisticSeqFromGenome(self, chromosome, start, end, strand):
		seqLength = self.getChromosomeLength(chromosome)
		if (start > 0 and start < seqLength + 1 and end < seqLength + 1):
			if(strand == '+'):
				return self.getStatisticSequence(self.chromosomeSeq[chromosome][start - 1:end])
			else:
				reverse = self.chromosomeSeq[chromosome][start - 1:end]
				reverse = self.complementary(reverse[::-1])
				return self.getStatisticSequence(reverse)
		else:
			print("Out of length in seq please check again")
			print("chromosome", chromosome, "length:", seqLength)
			print("gene position:", start, "to", end, "on", strand, "strand")
			exit()
	def getChromosomeLength(self, chromosome_name):
		return self.chromosomeLength[chromosome_name]
	def getSequence(self, chromosome, start, end, strand):
		seqLength = self.getChromosomeLength(chromosome)
		if (start > 0 and start < seqLength + 1 and end < seqLength + 1):
			if(strand == '+'):
				return self.chromosomeSeq[chromosome][start - 1:end]
			else:
				reverse = self.chromosomeSeq[chromosome][start - 1:end]
				reverse = self.complementary(reverse[::-1])
				return reverse
		else:
			return False
			print("\nOut of length in seq please check again")
			print("chromosome length:", seqLength)
			print("command: getSequence(", chromosome, start, end, strand, ") error!!!!", sep=', ')
			exit()
	def complementary(self, seq):
		new = ""
		for base in seq:
			if(base == 'A'):
				new = new + 'T'
			elif(base == 'T'):
				new = new + 'A'
			elif(base == 'G'):
				new = new + 'C'
			elif(base == 'C'):
				new = new + 'G'
			elif(base == 'a'):
				new = new + 't'
			elif(base == 't'):
				new = new + 'a'
			elif(base == 'g'):
				new = new + 'c'
			elif(base == 'c'):
				new = new + 'g'
			else:
				new = new + base
		return new
	def searchSeqInChromosome(self, chromosome_name, pattern):
		pattern = pattern.upper()
		len_pattern = len(pattern)
		index_found = []
		# Search pattern in plus strand
		index = self.chromosomeSeq[chromosome_name].find(pattern)
		while(index > -1):
			index_found.append([index + 1, index + len_pattern, '+'])
			index = self.chromosomeSeq[chromosome_name].find(pattern, index + 1)
		# Search pattern in minus strand
		pattern = self.complementary(pattern)[::-1]
		index = self.chromosomeSeq[chromosome_name].find(pattern)
		while(index > -1):
			index_found.append([index + 1, index + len_pattern, '-'])
			index = self.chromosomeSeq[chromosome_name].find(pattern, index + 1)
		# Return [fistMatch,endMatch,strand]
		return index_found
	def searchSeqInGenome(self, pattern):
		pattern = pattern.upper()
		len_pattern = len(pattern)
		index_found = []
		for chromosome_name, seq in sorted(self.chromosomeSeq.items()):
			# Search pattern in plus strand
			index = seq.find(pattern)
			while(index > -1):
				index_found.append([chromosome_name , index + 1, index + len_pattern, '+'])
				index = seq.find(pattern, index + 1)
			# Search pattern in minus strand
			pattern = self.complementary(pattern)[::-1]
			index = seq.find(pattern)
			while(index > -1):
				index_found.append([chromosome_name, index + 1, index + len_pattern, '-'])
				index = seq.find(pattern, index + 1)	
		return index_found

class Gff_manager(object):
	def __init__(self, file_name):
		self.data = []
		if(file_name.find('.gz') > 0):
			filegz = gzip.open(file_name, 'rb')
			gff_file = utf8(filegz)
		else:
			gff_file = open(file_name, 'r')
		for line in gff_file:
			if(line[0] != '#'):
				line = line.split()
				line[3] = int(line[3])
				line[4] = int(line[4])
				line[8] = line[8].split(';')
				if(line != ''):
					self.data.append(line)
	def getLastLineInGFF(self):
		return len(self.data)
	def printdata(self):
		countLine = 0
		for line in self.data:
			if(line[2] == 'five_prime_UTR' or line[2] == 'miRNA_primary_transcript'):
				print(line[0] + "\t" + line[2] + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + line[6] + "\t" + line[8][0])
				countLine += 1
	def getTable(self):
		return self.data
	def getTableSpecificType(self, geneType):
		table = []
		for line in self.data:
			if(line[2] == geneType):
				table.append(line)
		return table
	def getTableDataOfGene(self, geneName):
		table = []
		found = False
		breaker = False
		for line in self.data:
			if (line[2] == 'gene'):
				if(line[8][0].find(geneName) > 0):
					found = True
					breaker = True
				elif(breaker == True): 
					found = False
					break
			if(found == True):
				table.append(line)
		return table
	def getTableDataOfGeneAndType(self, geneName, type):
		table = []
		found = False
		breaker = False
		for line in self.data:
			if (line[2] == 'mRNA'):
				if(line[8][4].find(geneName,7) > 0 and line[8][3][8].find('1')):
					found = True
					breaker = True
				elif(breaker == True): 
					found = False
					break
			if(found == True):
				if(line[2] == type):
					table.append(line)
		return sorted(table,key=itemgetter(4,5))
	def getTranscripthave5UTR(self):
		print("gene", "transcript", "label5UTR", "lengthOf5UTR", "strand", "start", "stop", sep='\t')
		for line in self.data:
			if(line[2] == 'gene'):
				geneName = line[8][0][3:]
			elif(line[2] == 'five_prime_UTR'):
				transcriptName = line[8][0][3:26]
				label5UTR = line[8][0][-1:]
				start5UTR = int(line[3])
				stop5UTR = int(line[4])
				len5UTR = stop5UTR - start5UTR + 1
				strand = line[6]
#				 print(line)
				print(geneName, transcriptName, label5UTR, len5UTR, strand, start5UTR, stop5UTR, sep='\t')
	def getTranscripthave5UTRKU50(self):
		print("gene", "transcript", "label5UTR", "lengthOf5UTR", "strand", "start", "stop", sep='\t')
		for line in self.data:
			if(line[2] == 'gene'):
				geneName = line[8][0][3:]
			elif(line[2] == '5-UTR'):
				transcriptName = line[8][0][3:11]
				label5UTR = line[8][0][-1:]
				start5UTR = int(line[3])
				stop5UTR = int(line[4])
				len5UTR = stop5UTR - start5UTR + 1
				strand = line[6]
#				 print(line)
				print(geneName, transcriptName, label5UTR, len5UTR, strand, start5UTR, stop5UTR, sep='\t')
	def getGene(self):
		for line in self.data:
			if(line[2] == 'gene'):
				print(line[8][0][3:])
	def getDataSpecificType(self,gene_component):
		table = []
		for line in self.data:
			if(line[2] == gene_component):
				table.append(line)
		return table
	def getTranscript(self):
		for line in self.data:
			if(line[2] == 'mRNA'):
				print(line[8][0][3:])
	def isIntergenic(self, chromosome, start, stop, strand):
		isIntergenic = True
		listIntergenic = []
		if(strand == '+'):
			nearlyStart = [0, ""]
			nearlyEnd = [1000000000, ""]
			for gene in self.data:
				if(gene[6] == '+' and chromosome == gene[0] and gene[2] == 'gene'):
					# Check near front gene of miRNA
					if(gene[3] < start and gene[4] < stop):
						if(nearlyStart[0] < gene[4]):
							nearlyStart = [gene[4], gene[8][0], gene[4]]
					# Check near end gene of miRNA
					elif(gene[3] > start and gene[4] > stop):
						if(nearlyEnd[0] > gene[3]):
							nearlyEnd = [gene[3], gene[8][0], gene[3]]
					# Is Overlap with host gene show gene, start, stop
					else:
					# print(gene[8][0])
						listIntergenic.append([gene[8][0], gene[3], gene[4]])
		else:
			nearlyStart = [1000000000, ""]
			nearlyEnd = [0, ""]
			for gene in self.data:
				if(gene[6] == '-' and chromosome == gene[0] and gene[2] == 'gene'):
					# Check near front gene of miRNA
					if(gene[3] < start and gene[4] < stop):
						if(nearlyEnd[0] < gene[4]):
							nearlyEnd = [gene[4], gene[8][0], gene[4]]
					# Check near end gene of miRNA
					elif(gene[3] > start and gene[4] > stop):
						if(nearlyStart[0] > gene[3]):
							nearlyStart = [gene[3], gene[8][0], gene[3]]
					# Is Overlap with host gene show gene, start, stop
					else:
					# print(gene[8][0])
						listIntergenic.append([gene[8][0], gene[3], gene[4]])

		if(len(listIntergenic) == 0):
			if(nearlyStart[0] == 0 or nearlyStart[0] == 1000000000):
				nearlyStart = ''
			if(nearlyEnd[0] == 1000000000 or nearlyEnd[0] == 0):
				nearlyEnd = ''
			# Return intergenic list of naibor miRNA
			return ["", nearlyStart, nearlyEnd]
		else:
			return [listIntergenic, "", ""]
	def selectingPreMicroRNA(self):
		# diviced into 2 part
		# 1) selecting representative of pre-miRNA based on mature overlaping
		#	When overlaing between mature, then grouping into representative mature
		# 2) selecting representative of pre-miRNA based on precursor overlaping
		#	When overlaping between precursor, then grouping into representative precursor miRNA

		# First part,
		# 1.1) Create 'mi_loc' data table of of miRNA precursor and mature location on genomes
		# |	 0	 |	1	|	2	 |	3	|	 4	 |	5	 |	 6	|	7	|	 8		|
		# |P_Ma_ID|chr	|source|pre_start|pre_end|ma_start|ma_end|strand|miRNA_name|
		mi_loc = []
		for g in self.data:
			if(g[2] == 'miRNA_primary_transcript'):
				data = [g[8][0][3:],
						g[0],
						g[1],
						g[3],
						g[4],
						0,
						0,
						g[6],
						g[8][2][5:]]
			elif(g[2] == 'miRNA_mature'):
				data[0] += '_' + g[8][0][3:]
				data[5] = g[3]
				data[6] = g[4]
				mi_loc.append(data)
		mi_loc = sorted(mi_loc, key=itemgetter(7, 1, 5))	#order by strand,chr,ma_start
		# # Print list input mature and precursor data
		# for i in range(len(mi_loc)):
		# 	 print(i+1, mi_loc[i][0], mi_loc[i][2], mi_loc[i][1], mi_loc[i][7],sep='\t')
		
		# 1.2) Create 'ma_group' table of each group of mature overlaped
		# Grouping mature-miRNA 17 nt overlapped with other mature annotated
		# |	 0	 |	1	|	2	 |	3	|	 4	 |	5	 |	 6	|	7	|	 8		|
		# |P_Ma_ID|chr	|source|pre_start|pre_end|ma_start|ma_end|strand|miRNA_name|
		ma_group = []
		count = 0
		i = 0
		while(i < len(mi_loc)):
			if(mi_loc[i][5] < mi_loc[i][3] or mi_loc[i][6] > mi_loc[i][4]):
				# print("Warning! locate of mature miRNA incorrect:", mi_loc[i])
				count += 1
				ma_group.append([mi_loc[i]])
				i += 1
			else:
				j = i + 1
				in_group = [mi_loc[i]]
				while(j < len(mi_loc) and mi_loc[i][7] == mi_loc[j][7] and mi_loc[i][1] == mi_loc[j][1] and mi_loc[i][6] - mi_loc[j][5] + 1 >= 15):
					in_group.append(mi_loc[j])
					cover = mi_loc[i][6] - mi_loc[j][5] + 1
					if(cover > (mi_loc[j][6] - mi_loc[j][5] + 1)):
						cover = mi_loc[j][6] - mi_loc[j][5] + 1
					j += 1
				ma_group.append(in_group)
				count += 1
				i = j
		for i in ma_group:
			print(i)
		# 1.3) Selecting representative of miRNA precursor based on mature overlaping
		# By create 'list_pre_miRNA' table of best precursor can be represent representative mature region
		# |	 0	| 1 |	2	|	3	|	 4	 |		 5		 |			6		|		7		 |	 8		 |
		# |pre_id|chr|strand|pre_start|pre_end|[[ma_start,ma_end]]|[list_of_mature_ID]|source_from_pre_ID|length_of_pre|
		list_pre_miRNA = []
		for i in range(len(ma_group)): 
			if(len(ma_group[i]) == 1):
				# if unique of mature and precursor miRNA
				list_pre_miRNA.append([i,
										 ma_group[i][0][1],
										 ma_group[i][0][7],
										 ma_group[i][0][3],
										 ma_group[i][0][4],
										 [ma_group[i][0][5], ma_group[i][0][6]],
										 [ma_group[i][0][0]],
										 ma_group[i][0][0],
										 ma_group[i][0][4] - ma_group[i][0][3] + 1,
										 [ma_group[i][0][8]],
										 [ma_group[i][0][0]]
										 ]
										)
			else:
				# if multiple of mature and precursor miRNA
				strnad_of_precursor = ma_group[i][0][7]
				ma_min = ma_group[i][0][5]
				ma_max = ma_group[i][0][6]
				list_order_nearest_with_pre_miRNA = []
				list_of_mature_ID = []
				for j in range(len(ma_group[i])):
					if ma_min > ma_group[i][j][5]: ma_min = ma_group[i][j][5]
					if ma_max < ma_group[i][j][6]: ma_max = ma_group[i][j][6]
					list_of_mature_ID.append(ma_group[i][j][0])
				# selected pre-miRNA nearlest pre-miRNA to representative pre-miRNA
				list_ma_name_in_group = []
				for pre in ma_group[i]:
					if pre[8] not in list_ma_name_in_group:
						list_ma_name_in_group.append(pre[8])
					if(strnad_of_precursor=='+'):
						length = ma_min - pre[3]
					else:
						length = pre[4] - ma_max
					if(pre[3] <= ma_min and pre[4] >= ma_max):
						list_order_nearest_with_pre_miRNA.append([length, pre[0], pre[3], pre[4]])
				# sorted precursor which nearest with pre-miRNA
				# print("representative mature location min:", ma_min, ", max:", ma_max)
				if(len(list_order_nearest_with_pre_miRNA)==0):
					print("cann't be selecting representative pre-miRNA")
				else:
					list_order_nearest_with_pre_miRNA = sorted(list_order_nearest_with_pre_miRNA, key=itemgetter(0, 1))
				
				list_pre_miRNA.append([i,
										 ma_group[i][0][1],
										 ma_group[i][0][7],
										 list_order_nearest_with_pre_miRNA[0][2],
										 list_order_nearest_with_pre_miRNA[0][3],
										 [ma_min, ma_max],
										 list_of_mature_ID,
										 list_order_nearest_with_pre_miRNA[0][1],
										 list_order_nearest_with_pre_miRNA[0][3] - list_order_nearest_with_pre_miRNA[0][2] + 1,
										 list_ma_name_in_group,
										 list_of_mature_ID
										 ])
		# # show data selected from step overlaping mature on precursor
		# for i in list_pre_miRNA:
		#	 if(len(i[10])>2):
		#		 print(i[0]+1, i[7], sep='\t')
		#	 else:
		#		 print(i[0]+1, i[7], sep='\t')

		# select the best precursor and merge precursor
		# |	 0	 |	1	|	2	 |	3	|	 4	 |	5	 |	 6	|	7	|	 8		|
		# |PRE_ID|chr	|source|pre_start|pre_end|ma_start|ma_end|strand|miRNA_name|
		list_hairpin = []
		count = 1
		i = 0
		while(i < len(list_pre_miRNA)):
			if(list_pre_miRNA[i][1] == list_pre_miRNA[i - 1][1] and 
				list_pre_miRNA[i][2] == list_pre_miRNA[i - 1][2] and
				list_pre_miRNA[i - 1][5][0] >= list_pre_miRNA[i - 1][3] and
				list_pre_miRNA[i][5][0] >= list_pre_miRNA[i][3] and
				(list_pre_miRNA[i - 1][3] - list_pre_miRNA[i][4] <= 10 and list_pre_miRNA[i - 1][4] - list_pre_miRNA[i][3] >= 20)
				):
				list_name_of_pre = []
				for name in list_pre_miRNA[i-1][9]:
					if name not in list_name_of_pre:
						list_name_of_pre.append(name)
				for name in list_pre_miRNA[i][9]:
					if name not in list_name_of_pre:
						list_name_of_pre.append(name)
				list_id_of_pre = []
				for id in list_pre_miRNA[i-1][10]:
					list_id_of_pre.append(id)
				for id in list_pre_miRNA[i][10]:
					list_id_of_pre.append(id)
					
				list_select = sorted([list_pre_miRNA[i - 1], list_pre_miRNA[i]], key=itemgetter(8, 2))
				
				ma_min = list_pre_miRNA[i - 1][5][0]
				ma_max = list_pre_miRNA[i][5][1]
				# print(count)
				# count += 1
				# print(list_pre_miRNA[i - 1])
				# print(list_pre_miRNA[i])
				if(list_select[0][3] <= ma_min and list_select[0][4] >= ma_max):
					merge = list_select[0]
					merge[9] = list_name_of_pre
					merge[10] = list_id_of_pre
					merge[5].extend(list_select[1][5])
					merge[6].extend(list_select[1][6])
					if(list_hairpin[len(list_hairpin)-1][0]==merge[0] or list_hairpin[len(list_hairpin)-1][0]==list_select[1][0]):
						list_hairpin.pop()
						# print("pop",merge[0])
					list_hairpin.append(merge)
					# print(merge)
				elif(list_select[1][3] <= ma_min and list_select[1][4] >= ma_max):
					merge = list_select[1]
					merge[9] = list_name_of_pre
					merge[10] = list_id_of_pre
					merge[6].extend(list_select[0][6])
					merge[5].extend(list_select[0][5])
					if(list_hairpin[len(list_hairpin)-1][0]==merge[0] or list_hairpin[len(list_hairpin)-1][0]==list_select[0][0]):
						list_hairpin.pop()
#						 print("pop",merge[0])
					list_hairpin.append(merge)
				else:
					print("break")
					break
				i += 1
			else:
				list_hairpin.append(list_pre_miRNA[i])
				i += 1

		# # Print list of pre-miRNA selected
		# for i in list_hairpin:
		#	 for j in i[10]:
		#		 print(i[7], j, sep='\t')

		
		# Print to GFF format
		# chromosome00007	Ballen-Taborda_2013	miRNA_primary_transcript	3876	3971	.	+	.	ID=b415;Alias=;Name=sp-0426
		i = 0
		sum = 0
		for list in list_hairpin:
			i=i+1
			sum += len(list[6])
			print(list[1],"total", "miRNA_primary_transcript", list[3], list[4], '.',list[2], '.', "ID="+list[7]+";Alias=;Name="+list[7], sep="\t", end="\t")
			for name in list[9][:-1]: print(name, end=",")
			print(list[9][len(list[9])-1], end="\t")
			
			# for name in list[10][:-1]: print(name, end=",")
			# print(list[10][len(list[10])-1])
			for source in list[10]:
				print(source, end=",")

			print()
		# print(sum)
class Genome_manager(Fasta_manager, Gff_manager):
	def __init__(self, fastaFile, GffFile):
		self.fastaFile = fastaFile
		Fasta_manager.__init__(self, fastaFile)
		Gff_manager.__init__(self, GffFile)
	def getGCcontentInTranscript(self, type):
		sumGC = 0
		sumAT = 0
		for line in self.data:
			if(line[2] == type):
#				 print(line[8][0][3:], line[0], line[3], line[4] , line[6], sep='\t',end = '\t')
				statistic = Fasta_manager.getStatisticSeqFromGenome(self, line[0], line[3], line[4] , line[6])
#				 print(statistic[0], statistic[1], statistic[2], sep='\t')
				sumGC += statistic[1]
				sumAT += statistic[2]
		print("Summary GC content in", type, ":", float(sumGC) * 100 / (sumGC + sumAT))
	def selectedTSSProtein(self, upstream, downstream):
		file_write = open("%s_upstream_-%dto+%d.fa" % (self.fastaFile[:-6], upstream, downstream), 'w')
		statistic_of_5_prime_length = []
		geneListSelected = []
		geneCount = 0
		transcriptName = geneName = ''
		five_prime_UTR = []
		three_prime_UTR = []
		CDS = []
		count_five_prime_UTR_selected = 0
		count_five_prime_UTR_total = 0
		count_upstream_out_of_criteria = 0
		count_seq = 0
		for line in self.data:
			if(line[2] == 'gene'):
				geneName = line[8][0][3:]
				geneCount += 1
			elif(line[2] == 'mRNA'):
				count_five_prime = len(five_prime_UTR)
				if(count_five_prime > 0):
					count_five_prime_UTR_selected += 1
					count_five_prime_UTR_total += count_five_prime
					if geneName not in geneListSelected:
						geneListSelected.append(geneName)
					if(five_prime_UTR[0][6] == '+'):
						five_prime_UTR.sort(key=itemgetter (3, 4))
						selected_five_prime = five_prime_UTR[count_five_prime - 1]
					else:
						five_prime_UTR.sort(key=itemgetter (4, 3))
						selected_five_prime = five_prime_UTR[0]
					sequence = Fasta_manager.getSequence(self, selected_five_prime[0], selected_five_prime[3], selected_five_prime[4], selected_five_prime[6])
					statistic_of_5_prime_length.append(len(sequence))
#					 print(">", transcriptName, sep="")
#					 print(sequence)
					text = self.upstreamRetrive(upstream, downstream, selected_five_prime)
					if(text == False):
						count_upstream_out_of_criteria += 1
					else:
						file_write.writelines(text)
						count_seq += 1

				transcriptName = line[8][0][3:]
				five_prime_UTR = []
				three_prime_UTR = []
				CDS = []
			elif(line[2] == 'five_prime_UTR' or line[2] == '5-UTR'):
				five_prime_UTR.append(line)
			elif(line[2] == 'tree_prime_UTR' or line[2] == '3-UTR'):
				three_prime_UTR.append(line)
			elif(line[2] == 'CDS'):
				CDS.append(line)
		# lastLine imporve data
		count_five_prime = len(five_prime_UTR)
		if(count_five_prime > 0):
			count_five_prime_UTR_selected += 1
			count_five_prime_UTR_total += count_five_prime
			if geneName not in geneListSelected:
				geneListSelected.append(geneName)
			if(five_prime_UTR[0][6] == '+'):
				five_prime_UTR.sort(key=itemgetter (3, 4))
				selected_five_prime = five_prime_UTR[count_five_prime - 1]
			else:
				five_prime_UTR.sort(key=itemgetter (4, 3))
				selected_five_prime = five_prime_UTR[0]
			sequence = Fasta_manager.getSequence(self, selected_five_prime[0], selected_five_prime[3], selected_five_prime[4], selected_five_prime[6])
			statistic_of_5_prime_length.append(len(sequence))
#			 print(">", transcriptName, sep="")
#			 print(sequence)
			text = self.upstreamRetrive(upstream, downstream, selected_five_prime)
			if(text == False):
				count_upstream_out_of_criteria += 1
			else:
				file_write.writelines(text)
				count_seq += 1

		# Get statistic
		print("Statistic of genome", "%s_upstream_-%dto+%d.fa" % (self.fastaFile[:-6], upstream, downstream))
		print("Number of annotated gene:", geneCount)
		print("Number of 5'UTR of known gene:", len(geneListSelected))
		print("Number of alternative 5'UTR transcript:", count_five_prime_UTR_total)
		print("Number of selected 5'UTR transcript (unique):", count_five_prime_UTR_selected)
		print("Upstream correct:", count_seq)
		print("Upstream out of criteria:", count_upstream_out_of_criteria)
		# Number of 5'UTR of selected transcript
	def upstreamRetrive(self, upstream, downstream, five_prime_UTR):
		if(five_prime_UTR[6] == '+'):
			seq = Fasta_manager.getSequence(self, five_prime_UTR[0], five_prime_UTR[3] - upstream, five_prime_UTR[3] + downstream, five_prime_UTR[6])
		else:
			seq = Fasta_manager.getSequence(self, five_prime_UTR[0], five_prime_UTR[4] - downstream, five_prime_UTR[4] + upstream, five_prime_UTR[6])
			
		if(seq == False):
			return False
		else:
			if(seq.count('N') == 0):
				if(five_prime_UTR[6] == '+'):
					text = ">" + five_prime_UTR[8][0][3:] + "|" + str(five_prime_UTR[3] - upstream) + "|" + str(five_prime_UTR[3] + downstream) + "|+\n"
				else:
					text = ">" + five_prime_UTR[8][0][3:] + "|" + str(five_prime_UTR[4] - downstream) + "|" + str(five_prime_UTR[4] + upstream) + "|-\n"
				
				if(len(seq) != upstream + downstream + 1):
					print("\nLength of sequence not correct please check code it again.")
					exit()
				text = text + str(seq) + "\n"
#				 print(text)
				return text
			else:
				return False		
	def upstreamSeq(self, upstream, downstream):
		not_selected = 0
		not_selected_polyN = 0
		count_seq = 0
		for line in self.data:
			if(line[2] == 'five_prime_UTR'):
				if(line[6] == '+'):
					seq = Fasta_manager.getSequence(self, line[0], line[3] - upstream, line[3] + downstream, line[6])
				else:
					seq = Fasta_manager.getSequence(self, line[0], line[4] - downstream, line[4] + upstream, line[6])
				
				if(seq == False):
					not_selected += 1
				else:
					if(seq.count('N') == 0):
						if(line[6] == '+'):
							pass
#							 print(">", line[8][0][3:],"|",line[3]-upstream,"|",line[3]+downstream,"|+" ,sep='')
						else:
							pass
#							 print(">", line[8][0][3:],"|",line[4]-downstream,"|",line[4]+upstream, "|-" ,sep='')
						if(len(seq) != upstream + downstream + 1):
							print("\nLength of sequence not correct please check code it again.")
							exit()
#						 print(seq)
						count_seq += 1
					else:
						not_selected_polyN += 1
		print("not selected sequence:", not_selected)
		print("not selected sequence because N:", not_selected_polyN)
		print("It including ", count_seq, "sequences for next step")

# miRNAall = Gff_manager('Mesculenta_147_v4.1.miRNAall-testing.gff3')
miRNAall = Gff_manager('input/Mesculenta_147_v4.1.miRNAall.gff3')
miRNAall.selectingPreMicroRNA()
	




# class microRNA():
#	 def __init__(self, Genome_annotation, microRNA_annotation):
#		 self.genes = Gff_manager(Genome_annotation)
#		 self.miRNAs = Gff_manager(microRNA_annotation)
#		 self.intragenic_miRNA_identify('cassava4.1_008691m.g', 95500, 96000)
#	 def intragenic_miRNA_identify(self,geneName, pre_start, pre_end):
#		 table = self.genes.getTableDataOfGeneAndType(geneName, 'CDS')
#		 for i in range(len(table)): 
#			 print(table[i])
#			 if(table[i][4]<=pre_end):
#				 if table[i][6]=='+': return ['five_prime',{'start':}]
#			 else: before.append(before)
#
#			 print(x)
			
		
		# return 4 type
		# [type, {type_detail:[start,end]}]
		# type1 exonic_miRNAs: {in_exon:[start,stop]}
		# type2 intronic_miRNAs: {between_exonA:[], between_exonB:[]
		# type3 five_prime: {five_prime:[start,end]}
		# type4 tree_prime: {tree_prim:[start,end]}

# AM560MicroRNA = Gff_manager('Mesculenta_147_v4.1.miRNAall.gff3')
# AM560MicroRNA.selectingPreMicroRNA()

# AM560MicroRNA = microRNA('Mesculenta_147_v4.1.gene.gff3', 'Mesculenta_147_v4.1.miRNA.nature.gff3')
# AM560MicroRNA.classification_miRNA_loci()

# AM560v4 = Fasta_manager('Mesculenta_147_v4.1.fa')
# print(AM560v4.getSequence('chromosome04165', 19988, 20417, '-'))

# AM560Gffv4 = Gff_manager('Mesculenta_147_v4.1.gene.gff3')
# AM560Gffv4.getTable()
# AM560genome = Genome_manager('Mesculenta_147_v4.1.fa', 'Mesculenta_147_v4.1.gene.gff3.gz')
# AM560genome.getSequence(chromosome, start, end, strand)
# AM560genome.getGCcontentInTranscript('five_prime_UTR')
# AM560genome.getGCcontentInTranscript('three_prime_UTR')
# AM560genome.getGCcontentInTranscript('CDS')
# gene = Gff_manager('MW_v1.gene.gff')
# AM560genome.selectedTSSProtein(500, 500)
# AM560genome.selectedTSSProtein(1000, 1000)
# AM560genome.selectedTSSProtein(2000, 2000)
# AM560genome.selectedTSSProtein(5000, 5000)

# AM560v5 = Fasta_manager('Mesculenta_147_v5.1.fa')
# AM560Gffv5 = Gff_manage('Mesculenta_147_v5.1.gene.gff3')
# AM560genomev5 = Genome_manager('Mesculenta_147_v5.1.fa', 'Mesculenta_147_v5.1.gene.gff3')
# AM560genomev5.getGCcontentInTranscript('five_prime_UTR')
# AM560genomev5.getGCcontentInTranscript('three_prime_UTR')
# AM560genomev5.getGCcontentInTranscript('CDS')
# AM560genomev5.upstreamSeq(2000,2000)
# AM560genomev5.selectedTSSProtein()
# AM560genomev5.selectedTSSProtein(500, 500)
# AM560genomev5.selectedTSSProtein(1000, 1000)
# AM560genomev5.selectedTSSProtein(2000, 2000)
# AM560genomev5.selectedTSSProtein(5000, 5000)

# AM560v6 = Fasta_manager('cassavaV6_0.hardmasked.fasta')

# gene.getTable()
# gene.getTranscripthave5UTR()
# gene.getGene()
