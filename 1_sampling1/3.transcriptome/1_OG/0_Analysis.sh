echo "==========================================================="
echo "step 0 : KO-abundance to Pathway-enrich-score by ssGSEA"
echo "==========================================================="

convert_table.py --input <(csvtk -t join -f '2;1' /public2/zhushen/4_four_cockroach/3_transcriptome/5_OG/5_ssGSEA/OG_2_KO.txt /public2/zhushen/database/KEGG/1_KEGG_pathway/Insect_pathway/KEGG_pame_KO_2_pathway.txt | cut -f 1,3) --output /public2/zhushen/4_four_cockroach/3_transcriptome/5_OG/5_ssGSEA/KEGG_pathway_2_OG_wide.txt --mode long_to_wide  --gene-col protein --set-col pathway

gsva.R --input OG_tpm_salmon_QN.tsv --geneset /public2/zhushen/4_four_cockroach/3_transcriptome/5_OG/5_ssGSEA/KEGG_pathway_2_OG_wide.txt --output transcriptome_KEGG_pathway_ssGSEA_matrix.tsv --method ssgsea --normalize log


echo "==========================================================="
echo "step 1 : 'each group' VS 'the other three group'"
echo "==========================================================="

for i in Bger Esin Pame Pful
do 
	limma.R --input transcriptome_KEGG_pathway_ssGSEA_matrix.tsv --group-file group_${i}.tsv --group1 other --group2 ${i} --output ${i}_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv
done


echo "==========================================================="
echo "step 2 : 'each group' VS 'each group'"
echo "==========================================================="


for i in Bger Esin Pame Pful
do
	for j in Bger Esin Pame Pful
	do if [ "$i" != "$j" ]; then 
		limma.R --input transcriptome_KEGG_pathway_ssGSEA_matrix.tsv --group-file group.tsv --group1 ${i} --group2 ${j} --output ${j}_vs_${i}_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv
	fi
	done
done


echo "==========================================================="
echo "step 3 : intersection"
echo "==========================================================="

for i in up down;
do 
	paste	<(grep -w "${i}" Bger_vs_Esin_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iBger_vs_Esin_${i}") \
			<(grep -w "${i}" Bger_vs_Pame_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iBger_vs_Pame_${i}") \
			<(grep -w "${i}" Bger_vs_Pful_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iBger_vs_Pful_${i}") \
			> Bger_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done


for i in up down
do
	intersection.py --input Bger_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Bger_pathway_${i}_intersection.tsv
done

for i in up down
do 
	csvtk join -t -f '1;1;1' <(cat Bger_pathway_${i}_intersection.tsv | sed '1ipathway') Bger_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv transcriptome_KEGG_pathway_ssGSEA_matrix.tsv > Bger_pathway_${i}_intersection_pathway.tsv
done
	
echo "==========================================================="
echo "Bger KEGG ssGSEA pathway regulation analysis finished "
echo "==========================================================="


for i in up down
do
	paste	<(grep -w "${i}" Esin_vs_Bger_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iEsin_vs_Bger_${i}") \
			<(grep -w "${i}" Esin_vs_Pame_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iEsin_vs_Pame_${i}") \
			<(grep -w "${i}" Esin_vs_Pful_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iEsin_vs_Pful_${i}") \
			> Esin_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done

for i in up down
do
	intersection.py --input Esin_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Esin_pathway_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Esin_pathway_${i}_intersection.tsv | sed '1ipathway') Esin_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv transcriptome_KEGG_pathway_ssGSEA_matrix.tsv > Esin_pathway_${i}_intersection_pathway.tsv
done

echo "==========================================================="
echo "Esin KEGG ssGSEA pathway regulation analysis finished "
echo "==========================================================="


for i in up down
do
	paste	<(grep -w "${i}" Pame_vs_Bger_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iPame_vs_Bger_${i}") \
			<(grep -w "${i}" Pame_vs_Esin_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iPame_vs_Esin_${i}") \
			<(grep -w "${i}" Pame_vs_Pful_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iPame_vs_Pful_${i}") \
			> Pame_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done

for i in up down
do
	intersection.py --input Pame_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Pame_pathway_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Pame_pathway_${i}_intersection.tsv | sed '1ipathway') Pame_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv transcriptome_KEGG_pathway_ssGSEA_matrix.tsv > Pame_pathway_${i}_intersection_pathway.tsv
done

echo "==========================================================="
echo "Pame KEGG ssGSEA pathway regulation analysis finished "
echo "==========================================================="


for i in up down
do
	paste	<(grep -w "${i}" Pful_vs_Bger_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iPful_vs_Bger_${i}") \
			<(grep -w "${i}" Pful_vs_Esin_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iPful_vs_Esin_${i}") \
			<(grep -w "${i}" Pful_vs_Pame_transcriptome_KEGG_ssGSEA_pathway_limma_diff.tsv | cut -f 1 | sed "1iPful_vs_Pame_${i}") \
			> Pful_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv
done

for i in up down
do 
	intersection.py --input Pful_transcriptome_KEGG_ssGSEA_pathway_each_${i}.tsv --output Pful_pathway_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Pful_pathway_${i}_intersection.tsv | sed '1ipathway') Pful_enrich_transcriptome_KEGG_ssGSEA_pathway_limma.tsv transcriptome_KEGG_pathway_ssGSEA_matrix.tsv > Pful_pathway_${i}_intersection_pathway.tsv
done

echo "==========================================================="
echo "Pful KEGG ssGSEA pathway regulation analysis finished "
echo "==========================================================="



echo "================================================================================="
echo "step 4 : extract enrich & deplete pathways(information) in four species cockroach"
echo "================================================================================="


cat	<(csvtk -t sort -k 8:n Pful_pathway_up_intersection_pathway.tsv | sed 's/up/Pful_up/g' | sed 's/none/Pful_up/g') \
	<(csvtk -t sort -k 8:n Pame_pathway_up_intersection_pathway.tsv | sed 's/up/Pame_up/g' | sed 's/none/Pame_up/g' |sed '1d') \
	<(csvtk -t sort -k 8:n Bger_pathway_up_intersection_pathway.tsv | sed 's/up/Bger_up/g' | sed 's/none/Bger_up/g' |sed '1d') \
	<(csvtk -t sort -k 8:n Esin_pathway_up_intersection_pathway.tsv | sed 's/up/Esin_up/g' | sed 's/none/Esin_up/g' |sed '1d') \
	> 1_all_up_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 1_all_up_pathway.tsv)  /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 1_all_up_pathway_information.tsv

cat	<(csvtk -t sort -k 8:nr Pful_pathway_down_intersection_pathway.tsv | sed 's/down/Pful_down/g' | sed 's/none/Pful_down/g') \
	<(csvtk -t sort -k 8:nr Pame_pathway_down_intersection_pathway.tsv | sed 's/down/Pame_down/g' | sed 's/none/Pame_down/g' | sed '1d') \
	<(csvtk -t sort -k 8:nr Bger_pathway_down_intersection_pathway.tsv | sed 's/down/Bger_down/g' | sed 's/none/Bger_down/g' | sed '1d') \
	<(csvtk -t sort -k 8:nr Esin_pathway_down_intersection_pathway.tsv | sed 's/down/Esin_down/g' | sed 's/none/Esin_down/g' | sed '1d') \
	> 2_all_down_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 2_all_down_pathway.tsv)  /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 2_all_down_pathway_information.tsv



cat	<(head -n 1 Pful_pathway_up_intersection_pathway.tsv) \
	<(csvtk -t sort -k 8:n Pful_pathway_up_intersection_pathway.tsv | sed 's/up/Pful_up/g' | sed 's/none/Pful_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Pame_pathway_up_intersection_pathway.tsv | sed 's/up/Pame_up/g' | sed 's/none/Pame_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Bger_pathway_up_intersection_pathway.tsv | sed 's/up/Bger_up/g' | sed 's/none/Bger_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Esin_pathway_up_intersection_pathway.tsv | sed 's/up/Esin_up/g' | sed 's/none/Esin_up/g' | sed '1d' | tail -n 10) \
	> 3_each_top10_up_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 3_each_top10_up_pathway.tsv)  /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 3_each_top10_up_pathway_information.tsv


cat	<(head -n 2 Pful_pathway_down_intersection_pathway.tsv) \
	<(csvtk -t sort -k 8:nr Pful_pathway_down_intersection_pathway.tsv | sed 's/down/Pful_down/g' | sed 's/none/Pful_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Pame_pathway_down_intersection_pathway.tsv | sed 's/down/Pame_down/g' | sed 's/none/Pame_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Bger_pathway_down_intersection_pathway.tsv | sed 's/down/Bger_down/g' | sed 's/none/Bger_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Esin_pathway_down_intersection_pathway.tsv | sed 's/down/Esin_down/g' | sed 's/none/Esin_down/g' | sed '1d' | tail -n 10) \
	> 4_each_top10_down_pathway.tsv

csvtk join -t -f '1' <(cut -f 1 4_each_top10_down_pathway.tsv)  /public2/zhushen/database/KEGG/1_KEGG_pathway/3_KEGG_pathway_information.tsv > 4_each_top10_down_pathway_information.tsv

csvtk join -t -f '1;1' group_Esin.tsv  <(tsv-utils transpose transcriptome_KEGG_pathway_ssGSEA_matrix.tsv) > 5_lasso_input.tsv


limma.R --input transcriptome_KEGG_pathway_ssGSEA_matrix.tsv --group-file group_ternary.tsv --group1 Esin --group2 Peri --output 6_Peri_vs_Esin_metagenome_KEGG_ssGSEA_DEP_limma.tsv

limma.R --input transcriptome_KEGG_pathway_ssGSEA_matrix.tsv --group-file group_ternary.tsv --group1 Esin --group2 Bger --output 6_Bger_vs_Esin_metagenome_KEGG_ssGSEA_DEP_limma.tsv

csvtk join -t -f 1 \
	<(cut -f 1,8 6_Peri_vs_Esin_metagenome_KEGG_ssGSEA_DEP_limma.tsv | sed '1s/Directionality x -log10(p.adj)/Peri_vs_Esin/') \
	<(cut -f 1,8 6_Bger_vs_Esin_metagenome_KEGG_ssGSEA_DEP_limma.tsv | sed '1s/Directionality x -log10(p.adj)/Bger_vs_Esin/') \
	> 6_double_volcano_log2trans.tsv

