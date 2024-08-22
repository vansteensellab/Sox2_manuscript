# Call CTCF motifs from existing ChIP-Seq peaks
# Vinícius H. Franceschini-Santos, 2024-01-22
# ========================================================================================

# Here I'll use peaks that were already called by Tom and will call CTCF motifs on them
# using the MotifScan tool.
#
# First step is to convert the peaks to BED format, which is the format required by motifscan.
# But before I'll just define some variables and activate the environment.
# 0) Define variables and activate environment
# ========================================================================================
HERE=/DATA/usr/v.franceschini/Workspaces/2024_01_MATHIAS_CTCF_PEAKS/VF240122_calling_CTCF_motifs_from_previous_peaks/
PEAKS_DIR=${HERE}/01_INPUT/peaks
conda activate ctcf_peak_calling

# 1) Convert peaks to BED format
# ========================================================================================
# The output BED will contain the rows: chr, start, end, peak_name, score(q-value), strand(.)
# Note that I'll also remove the first row of the peaks which contain the sequence from ME20220215_F1CBS16_SB(lb02Fw)_Construct_KnockIN
for file in 6764_1_ME_E1_ATCACG_S1_R1_001_peaks 6764_2_ME_E2_CGATGT_S2_R1_001_peaks 6764_3_ME_E3_TTAGGC_S3_R1_001_peaks
    do
        grep -v "^#" ${PEAKS_DIR}/${file}.xls \
            | sed '1,3d' \
            | awk -v OFS="\t" -v FS="\t" '{ print $1, $2, $3, $10, $9, "." }' \
            > ${HERE}/01_INPUT/${file}.bed
    done

# 2) Call CTCF motifs using MotifScan
# ========================================================================================

# Before starting, I need to configure the motifscan tool.
# As I did it already for Lise, I'll just use the configuration file I created for her.

MOTIFSCAN_DIR="/DATA/usr/v.franceschini/Workspaces/2023_11_LISE_MOTIFS_IN_LBR/01_MOTIFSCAN_TOOL/motifscan_misc"
motifscan config --set-default-genome ${MOTIFSCAN_DIR}
motifscan config --set-default-motif ${MOTIFSCAN_DIR}

# Now, add the CTCF motif to the database
motifscan motif --install -n CTCF_mouse \
    -i ${HERE}/01_INPUT/MA0139.1.jaspar \
    -g mm10
# note that the above command takes a while to run

# Now, I'll call the CTCF motifs on the peaks
for file in 6764_1_ME_E1_ATCACG_S1_R1_001_peaks 6764_2_ME_E2_CGATGT_S2_R1_001_peaks 6764_3_ME_E3_TTAGGC_S3_R1_001_peaks
    do
        motifscan scan -i ${HERE}/01_INPUT/${file}.bed -w 0 \
            -g mm10 -m CTCF_mouse --no-enrich --site \
            -t 50 --strand both -o ${HERE}/02_OUTPUTS/intermediate_files/${file}_motifs
    done

# 3) Renaming the files
# ========================================================================================
# I'll rename the files to make it easier to work with them later
for file in 6764_1_ME_E1_ATCACG_S1_R1_001_peaks 6764_2_ME_E2_CGATGT_S2_R1_001_peaks 6764_3_ME_E3_TTAGGC_S3_R1_001_peaks
    do
        cp ${HERE}/02_OUTPUTS/${file}_motifs/motif_sites/MA0139_1_CTCF_sites.bed \
        ${HERE}/02_OUTPUTS/CTCF_sites/${file/_peaks/}_CTCF_sites.bed
    done
