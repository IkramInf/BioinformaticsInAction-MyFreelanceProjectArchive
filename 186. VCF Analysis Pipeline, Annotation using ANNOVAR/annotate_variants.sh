#!/bin/bash

# Set paths
VCF_FILE="denovo_mutations.vcf"
ANNOVAR_DIR="$HOME/Installs/annovar"
DATABASE_DIR="$ANNOVAR_DIR/humandb"
THREADS=8
BUILD="hg19"

# Specify Databases
DATABASES="refGene,dbnsfp42a,clinvar_20210501"  # add/remove databases here
OPERATIONS="g,f,f"  # Adjust operations according to databases (g:gene-based, f:filter-based)

# Convert VCF into Annovar Supported Format
perl $ANNOVAR_DIR/convert2annovar.pl -format vcf4 -allsample -withfreq $VCF_FILE > denovo_mutations.avinput

# Download databases if they don't exist
IFS=',' read -ra DB_ARRAY <<< "$DATABASES"
for DB in "${DB_ARRAY[@]}"; do
    if [ ! -f "$DATABASE_DIR/${BUILD}_${DB}.txt" ]; then
        perl $ANNOVAR_DIR/annotate_variation.pl -buildver $BUILD -downdb -webfrom annovar $DB $DATABASE_DIR
    else
        echo "Database ${BUILD}_${DB}.txt exists in $DATABASE_DIR. Skipping download..."
    fi
done

# Annotate the VCF File
perl $ANNOVAR_DIR/table_annovar.pl denovo_mutations.avinput $DATABASE_DIR \
    -buildver $BUILD \
    -out denovo_mutations_annotated \
    -remove \
    -protocol $DATABASES \
    -operation $OPERATIONS \
    -nastring . \
    -vcfinput \
    -thread $THREADS

echo "Annotation completed. Check the output files: denovo_mutations_annotated.*"

