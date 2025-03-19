#!/bin/bash

# GET COMMAND LINE ARGUMENTS
FASTQ_DIR="."
genome="hg19"
leading=0
trailing=0
crop=50

if [ "$#" -ge 3 ]; then
    FASTQ_DIR="$1"
    # Remove trailing slash, if present
    FASTQ_DIR=$(echo "$FASTQ_DIR" | sed 's/\/$//')
    genome="$2"
    crop="$3"
elif [ "$#" -ge 2 ]; then
    FASTQ_DIR="$1"
    # Remove trailing slash, if present
    FASTQ_DIR=$(echo "$FASTQ_DIR" | sed 's/\/$//')
    genome="$2"
else
    echo "Usage: $0 <FASTQ_DIR> <GENOME> <CROP>"
fi

echo "Fastq Directory is: $FASTQ_DIR"

# CREATE NECESSARY DIRECTORY TO SAVE GENOME, ADAPTER SEQUENCE, GENOME INDEX, PICARD TOOLS FILE
GENOMEDIR="$HOME/GENOMEDIR"
# Directories to create
DIRECTORIES=(
    "$GENOMEDIR"
    "$GENOMEDIR/hg19_index"
    "$GENOMEDIR/hg38_index"
    "$GENOMEDIR/picard"
    "$GENOMEDIR/adapter"
    "$GENOMEDIR/RefFlat"
)

# Check and create directories
for DIR in "${DIRECTORIES[@]}"; do
    if [ ! -d "$DIR" ]; then
        mkdir -p "$DIR"
        echo "Directory $DIR created."
    fi
done

# DOWNLOAD GENOME AND ANNOTATION FILES
hg19_FILES=(
    "hg19.fa.gz"
    "hg19.refGene.gtf.gz"
)
hg19_URLS=(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz"
)
hg38_FILES=(
    "hg38.fa.gz"
    "hg38.refGene.gtf.gz"
)
hg38_URLS=(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
)

if [ $genome = "hg19" ]; then
    for ((i=0; i<${#hg19_FILES[@]}; i++)); do
        FILE="${hg19_FILES[i]}"
        URL="${hg19_URLS[i]}"

        if [ ! -e "$GENOMEDIR/$FILE" ]; then
            echo "Downloading $FILE..."
            wget "$URL" -O "$GENOMEDIR/$FILE"
            echo "Downloaded $FILE to $GENOMEDIR"
        fi
    done
else
    for ((i=0; i<${#hg38_FILES[@]}; i++)); do
        FILE="${hg38_FILES[i]}"
        URL="${hg38_URLS[i]}"

        if [ ! -e "$GENOMEDIR/$FILE" ]; then
            echo "Downloading $FILE..."
            wget "$URL" -O "$GENOMEDIR/$FILE"
            echo "Downloaded $FILE to $GENOMEDIR"
        fi
    done
fi

# DOWNLOAD PICARD.JAR FILE
# Check if picard.jar exists, and if not, download it
if [ ! -e "$GENOMEDIR/picard/picard.jar" ]; then
    echo "Downloading picard.jar..."
    wget -O $GENOMEDIR/picard/picard.jar https://github.com/broadinstitute/picard/releases/download/2.26.2/picard.jar
    chmod +x $GENOMEDIR/picard/picard.jar
    echo "Downloaded picard.jar to $GENOMEDIR/picard"
fi

# DOWNLOAD ADAPTER SEQUENCES
# Adapter filenames
ADAPTERS=("TruSeq3-SE.fa" "TruSeq3-PE.fa")
# Check if each adapter file exists, and if not, download it
for ADAPTER in "${ADAPTERS[@]}"; do
    if [ ! -e "$GENOMEDIR/adapter/$ADAPTER" ]; then
        echo "Downloading $ADAPTER..."
        wget -O "$GENOMEDIR/adapter/$ADAPTER" "https://github.com/usadellab/Trimmomatic/raw/main/adapters/$ADAPTER"
        echo "Downloaded $ADAPTER to $GENOMEDIR/adapter"
    fi
done

# PERFORM GENOME INDEXING FOR STAR MAPPING
hg19_genome="$GENOMEDIR/hg19.fa.gz"
#hg19_gtf="$GENOMEDIR/hg19.refGene.gtf.gz"
hg19_gtf="$GENOMEDIR/genes.gtf"
hg38_genome="$GENOMEDIR/hg38.fa.gz"
#hg38_gtf="$GENOMEDIR/hg38.refGene.gtf.gz"
hg38_gtf="$GENOMEDIR/genes.gtf"

# Function to decompress files with -kf option
decompress_file() {
    local file="$1"
    local decompressed_file="${file%.gz}"

    if [[ -f "$decompressed_file" ]]; then
        echo "Decompressed file $decompressed_file already exists. Skipping decompression."
    elif [[ -f "$file" && "$file" =~ \.gz$ ]]; then
        echo "Decompressing $file..."
        gunzip -kf "$file"
        echo "Decompression completed."
    fi
}

#splice junction database overhang = max(ReadLength)-1
sjdboverhang=59
#genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1)
genomesaindexnbases=14

# Get number of threads for STAR
N=$(nproc)
numthreads=$((N - N / 3))

# Function to perform genome indexing
perform_indexing() {
    local genome_dir="$1"
    local genome_fasta="$2"
    local sjdb_gtf="$3"

    # Check if the index already exists
    if [ ! -e "$genome_dir/SAindex" ]; then
        # Index does not exist, so perform indexing
        STAR --runMode genomeGenerate --genomeDir "$genome_dir" --genomeFastaFiles "$genome_fasta" --sjdbGTFfile "$sjdb_gtf" --sjdbOverhang $sjdboverhang --genomeSAindexNbases $genomesaindexnbases --runThreadN $numthreads
    else
        echo "Genome index already exists in $genome_dir. Skipping indexing."
    fi
}

# Check if gtfToGenePred is present in the directory
if [ ! -f "$GENOMEDIR/gtfToGenePred" ]; then
    # If not present, download gtfToGenePred
    echo "Downloading gtfToGenePred..."
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -O "$GENOMEDIR/gtfToGenePred"
    # Make it executable
    chmod +x "$GENOMEDIR/gtfToGenePred"
    echo "gtfToGenePred downloaded and made executable."
else
    echo "gtfToGenePred is already present in the directory. Skipping download."
fi

# Function to convert gtf into gene prediction file
get_REF_FLAT_file() {
    local genome="$1"
    local refFlatFile="$2"

    # Check if refFlat file already exists
    if [ ! -f "$refFlatFile" ]; then
        # If not present, download
        wget "https://hgdownload.cse.ucsc.edu/goldenPath/$genome/database/refFlat.txt.gz" -O "$refFlatFile.gz"
        echo "Download completed, $refFlatFile..."
        decompress_file "$refFlatFile.gz"
    else
        echo "$refFlatFile already exists. Skipping download."
    fi
}

if [ "$genome" = "hg19" ]; then
    genomeFastaFile="$hg19_genome"
    sjdbGTFfile="$hg19_gtf"
    decompress_file $genomeFastaFile
    decompress_file $sjdbGTFfile
    genomeFastaFile="${genomeFastaFile%.gz}"
    sjdbGTFfile="${sjdbGTFfile%.gz}"
    perform_indexing "$GENOMEDIR/hg19_index" "$genomeFastaFile" "$sjdbGTFfile"
    refFlat="$GENOMEDIR/RefFlat/hg19_ref_flat.txt"
    get_REF_FLAT_file "hg19" "$refFlat"
    genomeDir="$GENOMEDIR/hg19_index"
elif [ "$genome" = "hg38" ]; then
    genomeFastaFile="$hg38_genome"
    sjdbGTFfile="$hg38_gtf"
    decompress_file $genomeFastaFile
    decompress_file $sjdbGTFfile
    genomeFastaFile="${genomeFastaFile%.gz}"
    sjdbGTFfile="${sjdbGTFfile%.gz}"
    perform_indexing "$GENOMEDIR/hg38_index" "$genomeFastaFile" "$sjdbGTFfile"
    refFlat="$GENOMEDIR/RefFlat/hg38_ref_flat.txt"
    get_REF_FLAT_file "hg38" "$refFlat"
    genomeDir="$GENOMEDIR/hg38_index"
else
    echo "Unsupported genome: $genome! Choose either hg19 or hg38..."
    exit 1
fi

PICARD="$GENOMEDIR/picard/picard.jar"
TRUSEQ="$GENOMEDIR/adapter/TruSeq3-SE.fa"
MAPPED="$GENOMEDIR/$genome"_mapped.tsv

# MAIN FUNCTION
# Iterate over all Read 1 FASTQ files in the directory
for R1_FASTQ_FILE in "$FASTQ_DIR"/*_R1*.fastq*; do
    echo $R1_FASTQ_FILE
    # Extract basename from Read 1 file
    filename=$(basename "$R1_FASTQ_FILE")
    fname="${filename%%_R1*.fastq*}"
    out_dir="$FASTQ_DIR/$fname"
    echo "$out_dir $fname"
    # Form the filename for Read 2
    R2_FASTQ_FILE="${R1_FASTQ_FILE/_R1/_R2}"
    echo $R2_FASTQ_FILE
    
    # Create a directory for each sample
    mkdir -p "$out_dir"

    # Check if Read 2 file exists
    if [ -e "$R2_FASTQ_FILE" ]; then
        # Paired-end (PE) processing
        # Perform FastQC for paired end
        fastqc "$R1_FASTQ_FILE" -o "$out_dir"
        fastqc "$R2_FASTQ_FILE" -o "$out_dir"
        # Perform Trimmomatic for paired-end
        trimmomatic PE -threads $numthreads -phred33 "$R1_FASTQ_FILE" "$R2_FASTQ_FILE" \
            "$out_dir/$fname"_trimmed_R1.fastq "$out_dir/$fname"_untrimmed_R1.fastq \
            "$out_dir/$fname"_trimmed_R2.fastq "$out_dir/$fname"_untrimmed_R2.fastq \
            ILLUMINACLIP:"$TRUSEQ":2:30:10 LEADING:$leading TRAILING:$trailing CROP:$crop SLIDINGWINDOW:4:15 MINLEN:15

        # Perform STAR alignment for paired-end
        prefix="$out_dir/$fname"_PE_
        echo $prefix
        STAR --runThreadN "$numthreads" \
             --genomeDir "$genomeDir" \
             --readFilesIn "$out_dir/$fname"_trimmed_R1.fastq "$out_dir/$fname"_trimmed_R2.fastq \
             --outFileNamePrefix $prefix \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --sjdbGTFfile "$sjdbGTFfile"

        input_bam="$out_dir/$fname_PE_Aligned.sortedByCoord.out.bam"
        gene_count_file="$out_dir/$fname_PE_ReadsPerGene.out.tab"

    else
        # Perform FastQC for single end
        fastqc "$R1_FASTQ_FILE" -o "$out_dir"
        # Perform Trimmomatic for single end
        trimmomatic SE -threads $numthreads -phred33 "$R1_FASTQ_FILE" "$out_dir/$fname".fastq \
            ILLUMINACLIP:"$TRUSEQ":2:30:10 LEADING:$leading TRAILING:$trailing \
            CROP:$crop SLIDINGWINDOW:4:15 MINLEN:10
        # Perform STAR alignment for single end
        prefix="$out_dir/$fname"_SE_
        echo $prefix
        STAR --runThreadN "$numthreads" \
             --genomeDir "$genomeDir" \
             --readFilesIn "$out_dir/$fname.fastq" \
             --outFileNamePrefix $prefix \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --sjdbGTFfile "$sjdbGTFfile"

        input_bam="$out_dir/$fname_SE_Aligned.sortedByCoord.out.bam"
        gene_count_file="$out_dir/$fname_SE_ReadsPerGene.out.tab"
    fi

    # replace ensemble id with gene name
    #awk 'BEGIN{FS=OFS="\t"} NR==FNR{gene[$1]=$2; next} {if($1 in gene) $1=gene[$1]} 1' $MAPPED "$gene_count_file" > "$gene_count_file"

    # multiqc
    multiqc -f "$out_dir" -o "$out_dir"
    # picard
    java -jar $PICARD CollectRnaSeqMetrics \
      I=$input_bam \
      O="$out_dir/output.RNA_Metrics" \
      REF_FLAT=$refFlat \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      MINIMUM_LENGTH=10

    java -jar $PICARD CollectAlignmentSummaryMetrics \
      R=$genomeFastaFile \
      I=$input_bam \
      O="$out_dir/output_mapping_metrics.txt"

    java -jar $PICARD MarkDuplicates \
      I=$input_bam \
      O="$out_dir/output_marked_duplicates.bam" \
      M="$out_dir/output_duplicate_metrics.txt"
done
