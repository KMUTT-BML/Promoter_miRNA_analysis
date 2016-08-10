#!/usr/bin/python

import os
import sys
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
This script designed for selecting representative miRNA from multiple source

Use as follows:
1) If you merge GFF before, Should type command following
$ python Representative_miRNA.py -o out.gff3 -i merge.gff3
2) Or for multiple GFF files
$ python Representative_miRNA.py -o out.gff3 -i source1.gff3 source2.gff3 source3.gff3 source4.gff3
"""

parser = OptionParser(usage=usage)
parser.add_option("-d", "--default", dest="default",
	default=None,
	help="Default input form thesis, please type 'thesis'")
parser.add_option("-o", "--output", dest="outputGff3",
	default=None, metavar="FILE",
	help="Output filename (required)")
parser.add_option("-i", "--input", dest="file_gff",
	default=None, metavar="FILE",
	help="GFF3 of miRNA, Gene feature file format (required), if you have multiple GFF please type followed")

options,args = parser.parse_args()

if options.default == 'thesis':
	outputGff3 = "input_for_classify_miRNA/rep_miRNA_out.gff3"
	file_gff = "input_find_representative_miRNA/Mesculenta_147_v4.1.miRNA_miRBased.gff3"
	args = ["input_find_representative_miRNA/Mesculenta_147_v4.1.miRNA_drougth_heat.gff3", "input_find_representative_miRNA/Mesculenta_147_v4.1.miRNA_chilling.gff3", "input_find_representative_miRNA/Mesculenta_147_v4.1.miRNA_cbb.gff3"]
else:
	if not options.outputGff3:
		sys.exit("Missing database type, -o ")
	outputGff3 = options.outputGff3
	if not options.file_gff:
		sys.exit("Missing database type, -i Mesculenta_147_v4.1.miRNAall.gff3")
	if not os.path.isfile(options.file_gff):
		sys.exit("Missing input file for Gff file: %r" % options.file_gff)
	file_gff = options.file_gff

miRNA_Genome = Gff_manager(file_gff, "miRNA_primary_transcript")
if len(args) > 1:
	for file_gff in args:
		miRNA_Genome.mergeGFFdata(file_gff, "miRNA_primary_transcript")

text = miRNA_Genome.selectingPreMicroRNA()
out_file = open(outputGff3, 'w')
out_file.write(text)




