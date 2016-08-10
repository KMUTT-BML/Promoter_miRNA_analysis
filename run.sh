echo -e "$(tput setab 7)$(tput setaf 2) =================================================================================== $(tput sgr 0)"
echo -e "$(tput setab 7)$(tput setaf 2)|                       Promoter analysis pipeline                                  |$(tput sgr 0)"
echo -e "$(tput setab 7)$(tput setaf 2) +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ $(tput sgr 0)\n"


out_dir=out
out_dir_step_1=${out_dir}/step_1
out_dir_step_2=${out_dir}/step_2
out_dir_step_3=${out_dir}/step_3

dir_soft_PromPredict=auxiliary/PromPredict_mulseq
# dir_soft_PromPredict=auxiliary/PromPredict_genome_V2
# # STEP 1 finding representative miRNA
# echo -e "$(tput setab 7)$(tput setaf 1) STEP 1 Finding representative miRNA $(tput sgr 0)" 
# mkdir -p ${out_dir}
# mkdir -p ${out_dir_step_1}
# python3 util/C1_representative_miRNA.py -o ${out_dir_step_1}/out.gff3 -i data/miRNA_annotation/Mesculenta_147_v4.1.miRNAall.gff3
# echo " + output: ${out_dir_step_1}/out.gff3"

# # STEP 2 classifying representative miRNA based on miRNA gene location
# echo -e "$(tput setab 7)$(tput setaf 1) STEP 2 classifying miRNA based on gene location $(tput sgr 0)" 
# mkdir -p ${out_dir_step_2}
# python3 util/C2_classify_miRNA.py \
# 	-i out/step_1/out.gff3 \
# 	--genome data/genome/Mesculenta_147_v4.1.fa \
# 	--gene_gff data/genome/Mesculenta_147_v4.1.gene.gff3 \
# 	-o out/step_2
# echo " + output: ${out_dir_step_2}/intergenic_miRNA.gff3"
# echo " + output: ${out_dir_step_2}/intronic_miRNA_.gff3"
# echo " + output: ${out_dir_step_2}/exonic_miRNA.gff3"


# STEP 3 Proximal promoter identification classifying
echo -e "$(tput setab 7)$(tput setaf 1) STEP 3 Proximal promoter identification$(tput sgr 0)" 
mkdir -p ${out_dir_step_3}
python3 util/C3_possible_promoter_regions.py \
	-i out/step_2/intergenic_miRNA.gff3 \
	--genome data/genome/Mesculenta_147_v4.1.fa \
	--gene_gff data/genome/Mesculenta_147_v4.1.gene.gff3 \
	-o out/step_3/possible_intergenic_miRNA_promoter.fa \
	--upstream 2000 \
	--downstream 200 \
	--min_len 400 \
	--Nfilter Y
echo " + output: out/step_3/possible_intergenic_miRNA_promoter.fa"
echo " + output: out/step_3/possible_intergenic_miRNA_promoter.fa.pos"

echo "$(tput setaf 2)Promoter region prediction by PromPredict$(tput sgr 0)"
p_E1_window_size=100
p_Genome_GC_content=35.29
${dir_soft_PromPredict} 1>&2 out/step_3/out_prompredict.log << EOF
out/step_3/possible_intergenic_miRNA_promoter.fa
${p_E1_window_size} 
${p_Genome_GC_content} 
EOF

python3 util/C4_parse_primary_miRNA.py 

python3 util/C3_possible_promoter_regions.py \
	-i out/step_3/primary_miRNA_loci.gff3 \
	--genome data/genome/Mesculenta_147_v4.1.fa \
	--gene_gff data/genome/Mesculenta_147_v4.1.gene.gff3 \
	-o out/step_3/proximal_intergenic_miRNA_promoter.fa \
	--upstream 1500 \
	--downstream 0 \
	--min_len 100 \
	--Nfilter N
echo " + output: out/step_3/proximal_intergenic_miRNA_promoter.fa"
echo " + output: out/step_3/proximal_intergenic_miRNA_promoter.fa.pos"