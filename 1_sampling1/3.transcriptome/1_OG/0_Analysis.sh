echo "==========================================================="
echo "step 0 : KO-abundance to Pathway-enrich-score by ssGSEA"
echo "==========================================================="

gsva.R --input transcript_tpm_filt.tsv --geneset /public2/zhushen/7_diff_point_PA/3_transcriptome/02_gene_kegg_annotation/merge/KEGG_pame_gene_2_pathway_wide.txt --output Total_transcriptome_KEGG_ssGSEA_pathway.tsv --method ssgsea


echo "==========================================================="
echo "step 1 : 'each group' VS 'the other two group'"
echo "==========================================================="

for i in Human Wild Lab
do 
	limma.R --input Total_transcriptome_KEGG_ssGSEA_pathway.tsv --group-file group_${i}.txt --group1 Other --group2 ${i} --output ${i}_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv
done


echo "==========================================================="
echo "step 2 : 'each group' VS 'each group'"
echo "==========================================================="


for i in Human Wild Lab
do
	for j in Human Wild Lab
	do if [ "$i" != "$j" ]; then 
		limma.R --input Total_transcriptome_KEGG_ssGSEA_pathway.tsv --group-file group_merge.txt --group1 ${i} --group2 ${j} --output ${j}_vs_${i}_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv
	fi
	done
done


echo "==========================================================="
echo "step 3 : intersection"
echo "==========================================================="

for i in up down;
do 
	paste	<(grep -w "${i}" Wild_vs_Lab_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iWild_vs_Lab_${i}") \
			<(grep -w "${i}" Wild_vs_Human_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iWild_vs_Human_${i}") \
			> Wild_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done


for i in up down
do
	intersection.py --input Wild_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Wild_pathway_${i}_intersection.tsv
done

for i in up down
do 
	csvtk join -t -f '1;1;1' <(cat Wild_pathway_${i}_intersection.tsv | sed '1ipathway') Wild_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv  Total_transcriptome_KEGG_ssGSEA_pathway.tsv  > Wild_pathway_${i}_intersection_pathway.tsv
done
	
echo "====================================================="
echo "Wild WGCNA module regulation analysis finished "
echo "====================================================="


for i in up down
do
	paste	<(grep -w "${i}" Lab_vs_Wild_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iLab_vs_Wild_${i}") \
			<(grep -w "${i}" Lab_vs_Human_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iLab_vs_Human_${i}") \
			> Lab_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done

for i in up down
do
	intersection.py --input Lab_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Lab_pathway_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Lab_pathway_${i}_intersection.tsv | sed '1ipathway') Lab_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv  Total_transcriptome_KEGG_ssGSEA_pathway.tsv > Lab_pathway_${i}_intersection_pathway.tsv
done

echo "===================================================="
echo "Lab WGCNA module regulation analysis finished "
echo "===================================================="


for i in up down
do
	paste	<(grep -w "${i}" Human_vs_Wild_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iHuman_vs_Wild_${i}") \
			<(grep -w "${i}" Human_vs_Lab_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iHuman_vs_Lab_${i}") \
			> Human_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done

for i in up down
do
	intersection.py --input Human_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Human_pathway_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Human_pathway_${i}_intersection.tsv | sed '1ipathway') Human_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv  Total_transcriptome_KEGG_ssGSEA_pathway.tsv  > Human_pathway_${i}_intersection_pathway.tsv
done

echo "======================================================"
echo "Human WGCNA module regulation analysis finished "
echo "======================================================"




echo "==============================================================================="
echo "step 4 : extract enrich & deplete modules(information) in three group cockroach"
echo "==============================================================================="


cat	<(csvtk -t sort -k 8:n Human_pathway_up_intersection_pathway.tsv | sed 's/up/Human_up/g' | sed 's/none/Human_up/g') \
	<(csvtk -t sort -k 8:n Lab_pathway_up_intersection_pathway.tsv | sed 's/up/Lab_up/g' | sed 's/none/Lab_up/g' |sed '1d') \
	<(csvtk -t sort -k 8:n Wild_pathway_up_intersection_pathway.tsv | sed 's/up/Wild_up/g' | sed 's/none/Wild_up/g' |sed '1d') \
	> 1_all_up_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 1_all_up_pathway.tsv) /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 1_all_up_pathway_information.tsv

cat	<(csvtk -t sort -k 8:nr Human_pathway_down_intersection_pathway.tsv | sed 's/down/Human_down/g' | sed 's/none/Human_down/g') \
	<(csvtk -t sort -k 8:nr Lab_pathway_down_intersection_pathway.tsv | sed 's/down/Lab_down/g' | sed 's/none/Lab_down/g' | sed '1d') \
	<(csvtk -t sort -k 8:nr Wild_pathway_down_intersection_pathway.tsv | sed 's/down/Wild_down/g' | sed 's/none/Wild_down/g' | sed '1d') \
	> 2_all_down_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 2_all_down_pathway.tsv) /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 2_all_down_pathway_information.tsv



cat	<(head -n 1 Wild_pathway_up_intersection_pathway.tsv) \
	<(csvtk -t sort -k 8:n Human_pathway_up_intersection_pathway.tsv | sed 's/up/Human_up/g' | sed 's/none/Human_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Lab_pathway_up_intersection_pathway.tsv | sed 's/up/Lab_up/g' | sed 's/none/Lab_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Wild_pathway_up_intersection_pathway.tsv | sed 's/up/Wild_up/g' | sed 's/none/Wild_up/g' | sed '1d' | tail -n 10) \
	> 3_each_top10_up_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 3_each_top10_up_pathway.tsv) /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 3_each_top10_up_pathway_information.tsv


cat	<(head -n 1 Wild_pathway_down_intersection_pathway.tsv) \
	<(csvtk -t sort -k 8:nr Human_pathway_down_intersection_pathway.tsv | sed 's/down/Human_down/g' | sed 's/none/Human_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Lab_pathway_down_intersection_pathway.tsv | sed 's/down/Lab_down/g' | sed 's/none/Lab_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Wild_pathway_down_intersection_pathway.tsv | sed 's/down/Wild_down/g' | sed 's/none/Wild_down/g' | sed '1d' | tail -n 10) \
	> 4_each_top10_down_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 4_each_top10_down_pathway.tsv) /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 4_each_top10_down_pathway_information.tsv

csvtk join -t -f '1;1' group_Wild.txt  <(tsv-utils transpose Total_transcriptome_KEGG_ssGSEA_pathway.tsv) > 5_lasso_input.tsv

limma.R --input Total_transcriptome_KEGG_ssGSEA_pathway.tsv --group-file group_merge.txt --group1 Wild --group2 Human --output 6_Human_vs_Wild_transcriptome_KEGG_ssGSEA_DEP_limma.tsv

limma.R --input Total_transcriptome_KEGG_ssGSEA_pathway.tsv --group-file group_merge.txt --group1 Wild --group2 Lab --output 6_Lab_vs_Wild_transcriptome_KEGG_ssGSEA_DEP_limma.tsv

csvtk join -t -f 1 \
	<(cut -f 1,8 6_Human_vs_Wild_transcriptome_KEGG_ssGSEA_DEP_limma.tsv | sed '1s/Directionality x -log10(p.adj)/Human_vs_Wild/') \
	<(cut -f 1,8 6_Lab_vs_Wild_transcriptome_KEGG_ssGSEA_DEP_limma.tsv | sed '1s/Directionality x -log10(p.adj)/Lab_vs_Wild/') \
	> 6_double_volcano_log2trans.tsv

