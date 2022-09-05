# CNNantigentic

Prerequisites
----------
TensorFlow 2.7.0 <br>
sklearn 0.19.2 <br>
pandas 1.3.5 <br>
numpy 1.21.6 <br>
keras 2.7.0 <br>


Getting started
----------
<br>
1. Run input_matrix_generation by CNNantigentic:<br><br>

`samtools mpileup -B -d 100 -f /path/to/ref.fasta [-l] [-r] -q 10 -O -s -a /path/to/tumor.bam /path/to/normal.bam | bgzip > /path/to/mixed_pileup_file`<br><br>Note: For the case of applying TransSSVs on a part of the whole genome, increase the BED entry by n (the number of flanking genomic sites to the left or right of the candidate somatic site) base pairs in each direction, and specify the genomic region via the option -l or -r. <br><br>
2. Run identi_candi_sites.py to identify candidate somatic small variants from the mixed pileup file: <br><br>` python3 identi_candi_sites.py
 --Tumor_Normal_mpileup /path/to/mixed_pileup_file
 --Candidate_somatic_sites /path/to/candidate_sites`<br><br>
3. Run mapping_infor_candi_sites.py to create a file with mapping information for candidate somatic small variant sites as input for trained TransSSVs, or to create a file with mapping information for validated somatic sites for training TransSSVs:<br><br> `python3 mapping_infor_candi_sites.py --Tumor_Normal_mpileup /path/to/mixed_pileup_file --Candidate_validated_somatic_sites /path/to/candidate_sites --number_of_columns N --length L --path /path/to/save --filename_1 filename_1 `<br><br>
4. Run model_train.py to train TransSSVs:<br><br>`CUDA_VISIBLE_DEVICES='' python3 model_train.py --input_dir /path/to/input --filename filename --vaild_dir /path/to/vaild_dir`<br><br>
5. Run model_infer.py to predict somatic small variants:<br><br>`CUDA_VISIBLE_DEVICES='' python3 model_infer.py --weights /path/to/weights --input_dir /path/to/input --filename filename --save_dir /results/TransSSVs`<br><br>`python3 write_vcf.py --vcf_file /path/to/vcf_file --pred_class /results/TransSSVs/y_pred_all.txt --Candidate_somatic_sites /path/to/candidate_sites`<br><br>


Example of the Validated_labels file for validated sites with labels (1: somatic site, 0: non-somatic site):<br><br>
        `chr1    790265  0`<br>
        `chr1    1595272 1`<br>
        `chr1    2312314 1`<br>
        `chr1    5006153 0`<br>


Please help us improve TransSSVs by reporting bugs or ideas on how to make things better. You can submit an issue or send me an email.<br>

Jing Meng, Jiangyuan Wang<br>

jing.mengrabbit@outlook.com<br>
wjy_bazi@hotmail.com<br>
