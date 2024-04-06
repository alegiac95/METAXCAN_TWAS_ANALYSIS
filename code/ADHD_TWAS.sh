#!/bin/bash
# :author: Alessio Giacomel
# :date: 2024-18-03
# :version: 1.0
# :disease: ADHD
# This scritp is used to run the TWAS prediction using S-PrediXcan and S-MetaXcan to predict the gene expression using the GTEx dataset.
# The GTEx data used are prediciton models from 12 brain tissues using both mashr models models.

# Directories setup
DATA_DIR="/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/data"
OUTPUT_DIR="/Users/alessiogiacomel/Dropbox/PhD/Analysis/transcriptomics_gio/MetaXcan_TWAS_Analysis/results/ADHD"
METAXCAN="/Users/alessiogiacomel/MetaXcan/software"
MODELS_DIR="${DATA_DIR}/models/elastic_net_models"
GWAS_FILE="${DATA_DIR}/GWAS/ADHD/ADHD2022_iPSYCH_deCODE_PGC.assoc.gz"

# Code setup 
# Check if the output folder already exists
if [ ! -d $OUTPUT_DIR ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p $OUTPUT_DIR
fi
# S-Predixcan of the models using both mashr models
# Loop over all the models in each model directory and run the S-PrediXcan
echo "<<<<<<<===========================================================================================================>>>>>>>"
echo "Running S-PrediXcan..."
for model in ${MODELS_DIR}/eqtl/*.db; do
    # save the file name without the extension
    model_name=$(basename "$model")
    model_name="${model_name%.*}"

    # Check if the output file already exists and if this exists skip the Spredixcan
    if [ ! -f $OUTPUT_DIR/en_PGC_ADHD_${model_name}.csv ]; then
        echo "INFO - Output file does not exist. Running S-PrediXcan for model: $model_name"
        echo ""
        python $METAXCAN/SPrediXcan.py \
            --model_db_path $MODELS_DIR/"eqtl"/$model_name".db" \
            --covariance $MODELS_DIR/"eqtl"/$model_name".txt.gz" \
            --gwas_file $GWAS_FILE \
            --snp_column SNP \
            --chromosome_column CHR \
            --position_column BP \
            --effect_allele_column A1 \
            --non_effect_allele_column A2 \
            --or_column OR \
            --se_column SE \
            --pvalue_column P \
            --model_db_snp_key rsid \
            --keep_non_rsid \
            --additional_output \
            --output_file $OUTPUT_DIR/en_PGC_ADHD_${model_name}.csv
        echo "INFO - Completed S-PrediXcan for model: $model_name"
        echo ""
        echo "--------------------------------------------------------------------------------------------------------------------"
        echo ""
    else
        echo "INFO - Output file already exists. Skipping S-PrediXcan for model: $model_name"
        echo ""
        echo "--------------------------------------------------------------------------------------------------------------------"
    fi

done
echo "<<<<<<<========================================================================================================>>>>>>>"
echo ""
# # S-MetaXcan of the models using all mashr models
echo "INFO - Running S-MultiXcan..."
echo ""
# Note: the model name pattern and metaxcan file name must match EXACTLY for S-MultiXcan to work
python $METAXCAN/SMultiXcan.py \
    --models_folder $MODELS_DIR/"eqtl/" \
    --models_name_pattern "en_(.*).db" \
    --models_name_filter "(.*)_Brain_(.*).db" \
    --snp_covariance ${MODELS_DIR}/"gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz" \
    --metaxcan_folder ${OUTPUT_DIR} \
    --metaxcan_file_name_parse_pattern "(.*)_PGC_ADHD_en_(.*).csv" \
    --gwas_file $GWAS_FILE \
    --snp_column SNP \
    --chromosome_column CHR \
    --position_column BP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --or_column OR \
    --se_column SE \
    --pvalue_column P \
    --keep_non_rsid \
    --cutoff_condition_number 30 \
    --model_db_snp_key rsid \
    --throw \
    --verbosity 3 \
    --output $OUTPUT_DIR/PGC_ADHD_SMetaXcan_en.tsv