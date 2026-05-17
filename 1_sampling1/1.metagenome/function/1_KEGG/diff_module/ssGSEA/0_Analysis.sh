echo "==========================================================="
echo "step 0 : KO-abundance to Module-enrich-score by ssGSEA"
echo "==========================================================="

gsva.R --input Total_KO_TPM_filt.tsv --geneset /public2/zhushen/database/KEGG/3_KEGG_module/Metagenome_KEGG_module_2_KO_wide.tsv --output Total_metagenome_KEGG_ssGSEA_module.tsv --method ssgsea --normalize log


echo "==========================================================="
echo "step 1 : 'each group' VS 'the other three group'"
echo "==========================================================="

for i in Bger Esin Pame Pful
do 
	limma.R --input Total_metagenome_KEGG_ssGSEA_module.tsv --group-file group_${i}.tsv --group1 other --group2 ${i} --output ${i}_enrich_metagenome_KEGG_ssGSEA_module_limma.tsv
done


echo "==========================================================="
echo "step 2 : 'each group' VS 'each group'"
echo "==========================================================="


for i in Bger Esin Pame Pful
do
	for j in Bger Esin Pame Pful
	do if [ "$i" != "$j" ]; then 
		limma.R --input Total_metagenome_KEGG_ssGSEA_module.tsv --group-file group.tsv --group1 ${i} --group2 ${j} --output ${j}_vs_${i}_metagenome_KEGG_ssGSEA_module_limma_diff.tsv
	fi
	done
done


echo "==========================================================="
echo "step 3 : intersection"
echo "==========================================================="

for i in up down;
do 
	paste	<(grep -w "${i}" Bger_vs_Esin_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iBger_vs_Esin_${i}") \
			<(grep -w "${i}" Bger_vs_Pame_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iBger_vs_Pame_${i}") \
			<(grep -w "${i}" Bger_vs_Pful_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iBger_vs_Pful_${i}") \
			> Bger_metagenome_KEGG_ssGSEA_module_each_${i}.tsv
done


for i in up down
do
	intersection.py --input Bger_metagenome_KEGG_ssGSEA_module_each_${i}.tsv --output Bger_module_${i}_intersection.tsv
done

for i in up down
do 
	csvtk join -t -f '1;1;1' <(cat Bger_module_${i}_intersection.tsv | sed '1imodule') Bger_enrich_metagenome_KEGG_ssGSEA_module_limma.tsv Total_metagenome_KEGG_ssGSEA_module.tsv > Bger_module_${i}_intersection_module.tsv
done
	
echo "==========================================================="
echo "Bger KEGG ssGSEA module regulation analysis finished "
echo "==========================================================="


for i in up down
do
	paste	<(grep -w "${i}" Esin_vs_Bger_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iEsin_vs_Bger_${i}") \
			<(grep -w "${i}" Esin_vs_Pame_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iEsin_vs_Pame_${i}") \
			<(grep -w "${i}" Esin_vs_Pful_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iEsin_vs_Pful_${i}") \
			> Esin_metagenome_KEGG_ssGSEA_module_each_${i}.tsv
done

for i in up down
do
	intersection.py --input Esin_metagenome_KEGG_ssGSEA_module_each_${i}.tsv --output Esin_module_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Esin_module_${i}_intersection.tsv | sed '1imodule') Esin_enrich_metagenome_KEGG_ssGSEA_module_limma.tsv Total_metagenome_KEGG_ssGSEA_module.tsv > Esin_module_${i}_intersection_module.tsv
done

echo "==========================================================="
echo "Esin KEGG ssGSEA module regulation analysis finished "
echo "==========================================================="


for i in up down
do
	paste	<(grep -w "${i}" Pame_vs_Bger_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iPame_vs_Bger_${i}") \
			<(grep -w "${i}" Pame_vs_Esin_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iPame_vs_Esin_${i}") \
			<(grep -w "${i}" Pame_vs_Pful_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iPame_vs_Pful_${i}") \
			> Pame_metagenome_KEGG_ssGSEA_module_each_${i}.tsv
done

for i in up down
do
	intersection.py --input Pame_metagenome_KEGG_ssGSEA_module_each_${i}.tsv --output Pame_module_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Pame_module_${i}_intersection.tsv | sed '1imodule') Pame_enrich_metagenome_KEGG_ssGSEA_module_limma.tsv Total_metagenome_KEGG_ssGSEA_module.tsv > Pame_module_${i}_intersection_module.tsv
done

echo "==========================================================="
echo "Pame KEGG ssGSEA module regulation analysis finished "
echo "==========================================================="


for i in up down
do
	paste	<(grep -w "${i}" Pful_vs_Bger_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iPful_vs_Bger_${i}") \
			<(grep -w "${i}" Pful_vs_Esin_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iPful_vs_Esin_${i}") \
			<(grep -w "${i}" Pful_vs_Pame_metagenome_KEGG_ssGSEA_module_limma_diff.tsv | cut -f 1 | sed "1iPful_vs_Pame_${i}") \
			> Pful_metagenome_KEGG_ssGSEA_module_each_${i}.tsv
done

for i in up down
do 
	intersection.py --input Pful_metagenome_KEGG_ssGSEA_module_each_${i}.tsv --output Pful_module_${i}_intersection.tsv
done

for i in up down
do
	csvtk join -t -f '1;1;1' <(cat Pful_module_${i}_intersection.tsv | sed '1imodule') Pful_enrich_metagenome_KEGG_ssGSEA_module_limma.tsv Total_metagenome_KEGG_ssGSEA_module.tsv > Pful_module_${i}_intersection_module.tsv
done

echo "==========================================================="
echo "Pful KEGG ssGSEA module regulation analysis finished "
echo "==========================================================="



echo "================================================================================="
echo "step 4 : extract enrich & deplete modules(information) in four species cockroach"
echo "================================================================================="


cat	<(csvtk -t sort -k 8:n Pful_module_up_intersection_module.tsv | sed 's/up/Pful_up/g' | sed 's/none/Pful_up/g') \
	<(csvtk -t sort -k 8:n Pame_module_up_intersection_module.tsv | sed 's/up/Pame_up/g' | sed 's/none/Pame_up/g' |sed '1d') \
	<(csvtk -t sort -k 8:n Bger_module_up_intersection_module.tsv | sed 's/up/Bger_up/g' | sed 's/none/Bger_up/g' |sed '1d') \
	<(csvtk -t sort -k 8:n Esin_module_up_intersection_module.tsv | sed 's/up/Esin_up/g' | sed 's/none/Esin_up/g' |sed '1d') \
	> 1_all_up_module.tsv

csvtk join -t -f '1' <(cut -f 1 1_all_up_module.tsv)  /public2/zhushen/database/KEGG/3_KEGG_module/Metagenome_KEGG_module_information_Total.tsv > 1_all_up_module_information.tsv

cat	<(csvtk -t sort -k 8:nr Pful_module_down_intersection_module.tsv | sed 's/down/Pful_down/g' | sed 's/none/Pful_down/g') \
	<(csvtk -t sort -k 8:nr Pame_module_down_intersection_module.tsv | sed 's/down/Pame_down/g' | sed 's/none/Pame_down/g' | sed '1d') \
	<(csvtk -t sort -k 8:nr Bger_module_down_intersection_module.tsv | sed 's/down/Bger_down/g' | sed 's/none/Bger_down/g' | sed '1d') \
	<(csvtk -t sort -k 8:nr Esin_module_down_intersection_module.tsv | sed 's/down/Esin_down/g' | sed 's/none/Esin_down/g' | sed '1d') \
	> 2_all_down_module.tsv

csvtk join -t -f '1' <(cut -f 1 2_all_down_module.tsv)  /public2/zhushen/database/KEGG/3_KEGG_module/Metagenome_KEGG_module_information_Total.tsv > 2_all_down_module_information.tsv



cat	<(head -n 1 Pful_module_up_intersection_module.tsv) \
	<(csvtk -t sort -k 8:n Pful_module_up_intersection_module.tsv | sed 's/up/Pful_up/g' | sed 's/none/Pful_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Pame_module_up_intersection_module.tsv | sed 's/up/Pame_up/g' | sed 's/none/Pame_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Bger_module_up_intersection_module.tsv | sed 's/up/Bger_up/g' | sed 's/none/Bger_up/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:n Esin_module_up_intersection_module.tsv | sed 's/up/Esin_up/g' | sed 's/none/Esin_up/g' | sed '1d' | tail -n 10) \
	> 3_each_top10_up_module.tsv

csvtk join -t -f '1' <(cut -f 1 3_each_top10_up_module.tsv)  /public2/zhushen/database/KEGG/3_KEGG_module/Metagenome_KEGG_module_information_Total.tsv > 3_each_top10_up_module_information.tsv


cat	<(head -n 2 Pful_module_down_intersection_module.tsv) \
	<(csvtk -t sort -k 8:nr Pful_module_down_intersection_module.tsv | sed 's/down/Pful_down/g' | sed 's/none/Pful_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Pame_module_down_intersection_module.tsv | sed 's/down/Pame_down/g' | sed 's/none/Pame_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Bger_module_down_intersection_module.tsv | sed 's/down/Bger_down/g' | sed 's/none/Bger_down/g' | sed '1d' | tail -n 10) \
	<(csvtk -t sort -k 8:nr Esin_module_down_intersection_module.tsv | sed 's/down/Esin_down/g' | sed 's/none/Esin_down/g' | sed '1d' | tail -n 10) \
	> 4_each_top10_down_module.tsv

csvtk join -t -f '1' <(cut -f 1 4_each_top10_down_module.tsv)  /public2/zhushen/database/KEGG/3_KEGG_module/Metagenome_KEGG_module_information_Total.tsv > 4_each_top10_down_module_information.tsv
