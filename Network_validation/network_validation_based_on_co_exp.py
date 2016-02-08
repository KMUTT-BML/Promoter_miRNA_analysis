# co-expression network between TF gene and miRNA
# Create by Nattawet Sriwichai, KMUTT, Thailand
# Release v.1
# Date 2015-12-27
# 1) finding cutoff DEG of TF and miRNA
# 2) calculate PCC and filltered |PCC| >= 0.8 to high correlation between TF and miRNA

import math
import statistics
import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set(color_codes=True)

node_tf = []
node_miRNA = []

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearson_cal(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for i in range(n):
        xdiff = x[i] - avg_x
        ydiff = y[i] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    if(xdiff2 * ydiff2 == 0):
    	return 0
    else:
    	return diffprod / math.sqrt(xdiff2 * ydiff2)

def cutoff_sd(list_sd, cutoff_sd_percentile):
	# # Plot distribution of SD
	# np.random.seed(sum(map(ord, "distributions of SD")))
	# x = list_sd
	# sns.distplot(x);
	# plt.show()
	# print(max(list_sd), min(list_sd))
	return np.percentile(list_sd, cutoff_sd_percentile)

def convert_csv_to_exp_table(file_name, num_condition, skip_header=-1):
	table = open(file_name).read().split("\n")[:skip_header]
	table[0] = table[0].split(",")
	table[0].append('stdev')
	dict_table = {}
	for i in range(1,len(table)):
		table[i] = table[i].split(",")
		for j in range(1,len(table[i])):
			table[i][j] = float(table[i][j])
		std = statistics.stdev(table[i][1:num_condition+1])
		table[i].append(std)
		dict_table[table[i][0]] = table[i][1:]
	table = []
	return dict_table
def convert_csv_to_interaction(file_name, skip_header=-1):
	# structure of csv file "miRNA_name|TF_protein_name"
	table = open(file_name).read().split("\n")[:skip_header]
	for i in range(len(table)):
		table[i] = table[i].split(',')[:2]
		if(table[i][0] not in node_miRNA):
			node_miRNA.append(table[i][0])
		if(table[i][1] not in node_tf):
			node_tf.append(table[i][1])
	return table

# 1) Read iteraction file and process
file_template_network = "interaction_info.csv"
network_interaction = convert_csv_to_interaction(file_template_network, -1)
print("number of interaction:", len(network_interaction))
print("Number of miRNA node:", len(node_miRNA))
print("Number of TF node   :", len(node_tf))
# print(node_tf)

# 2) Read expression profile and filltering good expression

# 2.1) TF expression profile filltering
# cutoff DEG based on 50 percentides of SD in tf_expression
n = 3 # number of time series condition
file_name_tf_expression = "mRNA_chilling_expression_all.csv"
tf_expression_table = convert_csv_to_exp_table(file_name_tf_expression, n, -1)
DEG_TF_sd = []
for tf in tf_expression_table:
	DEG_TF_sd.append(tf_expression_table.get(tf)[-1])
cutoff_DEG_TF_sd = cutoff_sd(DEG_TF_sd, 80)
print("Cutoff SD of mRNA DEG:", cutoff_DEG_TF_sd, "of ", len(DEG_TF_sd), "mRNA transcript")

define_CPM_value_expressed = 10
filltered_TF = {}
for tf in tf_expression_table:
	selected_TF = False
	for i in tf_expression_table[tf]:
		if i >= 10 :
			selected_TF = True
	if selected_TF and tf in node_tf and tf_expression_table[tf][-1]>=cutoff_DEG_TF_sd:
		filltered_TF[tf] = tf_expression_table[tf]
print("Number of filltered TF:", len(filltered_TF))

# 2.2) miTNA expression profile filltering
file_name_miRNA_expression = "miRNA_chilling_expression_qPCR_53.csv"
miRNA_expression_table = convert_csv_to_exp_table(file_name_miRNA_expression, n, -1)
# cutoff DEG based on 50 percentide of SD in miRNA_expression
defind_SD_percentide = 50
DEG_miRNA_sd = []
for key in miRNA_expression_table:
	DEG_miRNA_sd.append(miRNA_expression_table.get(key)[-1])
cutoff_DEG_miRNA_sd = cutoff_sd(DEG_miRNA_sd, defind_SD_percentide)
print("Cutoff SD of miRNA DEG:", cutoff_DEG_miRNA_sd)
# filltering miRNA expression profile
filltered_miRNA = {}
for miRNA in miRNA_expression_table:
	if miRNA in node_miRNA and miRNA_expression_table[miRNA][-1]>= cutoff_DEG_miRNA_sd:
		filltered_miRNA[miRNA] = miRNA_expression_table[miRNA]
print("Number of filltered miRNA:", len(filltered_miRNA))
# for i in filltered_miRNA:
# 	print(i, filltered_miRNA[i])


# 3) calculate co-expression correlation
pcc_cutoff = 0.8 # absolute of PCC cutoff
# calculate PCC of each TF and miRNA
count = 0
list_of_pcc = []
for interaction in network_interaction[1:]:
	if interaction[0] in filltered_miRNA and interaction[1] in filltered_TF:
		pcc = pearson_cal(filltered_miRNA.get(interaction[0])[:-2], filltered_TF.get(interaction[1])[:-2])
		list_of_pcc.append(pcc)
		# tf_DEG = "YES" if tf[n]>=cutoff_DEG_TF_sd else "NO"
		# miRNA_DEG = "YES" if miRNA[n]>=cutoff_DEG_miRNA_sd else "NO"
		correlation_type = "positive" if pcc>0 else "negative"
		# if ((tf_DEG == "YES" and miRNA_DEG == "YES") and abs(pcc)>=pcc_cutoff):
		# 	print(tf[0], miRNA[0], correlation_type, tf_DEG, miRNA_DEG, sep="\t")
		if (abs(pcc)>=0.8):
			print(interaction[0], interaction[1], pcc, correlation_type, sep="\t")
			count += 1
		else:
			print(interaction[0], interaction[1], pcc, "no correlation", sep="\t")
print("Number of confirmed interaction TF and miRNA: ", count, "of ", len(list_of_pcc), "interaction")

# Plot distribution of PCC
# np.random.seed(sum(map(ord, "PCC value TF-miRNA distribution ")))
# sns.distplot(list_of_pcc);
# plt.show()


