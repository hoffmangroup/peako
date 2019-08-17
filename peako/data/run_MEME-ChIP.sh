#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit
function finish {
    for bname in ${TFKO%.*} ${TF_file%.*} ${KOIN_file%.*}; do
	rm "$TF_DIR/$bname-expanded-edited.bed" "$TF_DIR/$bname-expanded-sequencesRemoved.bed" "$TF_DIR/$bname-expandedTo500bpRegions-slop.bed" "$TF_DIR/$bname-pre.bed" "$TF_DIR/fastaFromBed-$bname.err" || true
    done
} 
trap finish EXIT

ERR_EXIT=109

#########################################################
##  MEME-ChIP pipeline for WT, KO, and KOIN data sets  ##
#########################################################

# Note: files should end in .narrowPeak or .narrowPeak.gz

##### PARAMETERS

WT_FILE=''
KO_FILE=''
KOIN_FILE=''

usage() {
   echo -e "Usage: $0 \n -i <narrowPeak directory full path (input)> \n" \
       "-o <output directory> \n" \
       "-M <meme-chip command path (run_meme-chip-neg)> \n" \
       "-s <chromosome sizes file (mm9, mm10, or hg38)> \n" \
       "-g <TRF-masked genome FASTA file (mm9, mm10, or hg38)> \n" \
       "-t <wt narrowPeak filename (contained in input directory)> \n" \
       "-c <ko narrowPeak filename (contained in input directory)> \n" \
       "-k <KOIN narrowPeak filename (contained in input directory)> \n" \
       "-p <number of processors per run> \n" \
       "[-h for help]" \
       1>&2
    exit 1
}

while getopts ":i:o:M:s:g:t:c:k:h" opt; do
    case $opt in
        i)
            INPUT_DIR="${OPTARG:-.}/"
            if [[ ! -d "$INPUT_DIR" ]]; then
                >&2 echo "-i requires a directory."
                exit $ERR_EXIT
            fi
            ;;
        o)
            OUTPUT_DIR="${OPTARG:-.}/"
            ;;
        M)
            MEMECHIP_PATH="${OPTARG:-}"  # run_meme-chip-neg
            ;;
        s)
            CHR_SIZES="${OPTARG:-}"  # mm9, mm10, or hg38
            ;;
	g)
	    TRF_MASKED_GENOME="${OPTARG:-}"  # mm9, mm10, or hg38
	    ;;
        t)
            WT_FILE="${OPTARG:-}"
            ;;
        c)
            KO_FILE="${OPTARG:-}"
            ;;
        k)
            KOIN_FILE="${OPTARG:-}"
            ;;
        h)
            usage
            ;;
        \?)
            >&2 echo "Invalid option: -${OPTARG:-}"
            exit $ERR_EXIT
            ;;
        :)
            >&2 echo "Option -${OPTARG:-} requires an argument."
            exit $ERR_EXIT
            ;;
    esac
done

INPUT_DIR=$(realpath $INPUT_DIR)
mkdir $OUTPUT_DIR || true
OUTPUT_DIR=$(realpath $OUTPUT_DIR)
WT_FILE="$(basename -- $WT_FILE)"
KO_FILE="$(basename -- $KO_FILE)"
KOIN_FILE="$(basename -- $KOIN_FILE)"


# $f == $1 (*.narrowPeak)
function make_masked_fasta {
    cp "$1" "$TF_DIR"  # remove this copy at the end
    if [[ "${1##*.}" = "gz" ]]; then
	gzip -d "$TF_DIR/$1"
	bname=$(echo "$1" | sed -r 's/\.narrowPeak.gz//')
    else
	bname=$(echo "$1" | sed -r 's/\.narrowPeak//')	
    fi
    cd "$TF_DIR"
    echo $bname
    bed_pre="$bname-pre.bed"; awk 'BEGIN{FS=OFS="\t"} {midPos=($2+$10); print $1,midPos,midPos+1}' "$bname.narrowPeak" > "$TF_DIR/$bed_pre"
    bed_exp="$bname-expandedTo500bpRegions-slop.bed"; bedtools slop -i "$TF_DIR/$bed_pre" -g <(sort -V "$CWD_DIR/$CHR_SIZES" | sed "1ichrom\tsize") -l 250 -r 249 > "$TF_DIR/$bed_exp"
    awk '(($3-$2) == 500)' "$TF_DIR/$bed_exp" | awk '($1 != "chrM")' > "$TF_DIR/$bname-expanded-edited.bed" 
    awk '(($3-$2) != 500); ($1 == "chrM")' "$TF_DIR/$bed_exp" | awk '!x[$0]++' > "$TF_DIR/$bname-expanded-sequencesRemoved.bed"
    fastaFromBed -fi "$CWD_DIR/$TRF_MASKED_GENOME" -bed "$TF_DIR/$bname-expanded-edited.bed" -fo "$TF_DIR/$bname.fa" 2> "$TF_DIR/fastaFromBed-$bname.err"
}

function name_fasta {
    if [[ "${1##*.}" = "gz" ]]; then
        echo "${1%.narrowPeak.gz}.fa"
    else
        echo "${1%.*}.fa"
    fi
}

##### RUN

echo $(pwd)
CWD_DIR=$(pwd)
# create output subdirectory for TF
cd $OUTPUT_DIR  # condense into one line?
TF_DIR="$OUTPUT_DIR"

# extract genomic sequences from centered and extended peaks
for file in $WT_FILE $KO_FILE $KOIN_FILE; do  # allows for some to be unset
    cd $INPUT_DIR
    make_masked_fasta $file
done

# pull out .fa names
TFKO=$(if [[ -n "$KO_FILE" ]]; then name_fasta $KO_FILE; fi)
TF_file=$(if [[ -n "$WT_FILE" ]]; then name_fasta $WT_FILE; fi)
KOIN_file=$(if [[ -n "$KOIN_FILE" ]]; then name_fasta $KOIN_FILE; fi)
