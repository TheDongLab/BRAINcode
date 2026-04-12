#!/bin/bash
#SBATCH --job-name=prep_sQTL
#SBATCH --output=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.out
#SBATCH --error=/home/zw529/donglab/data/target_ALS/QTL/%x_%j.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=499G

set -euo pipefail
module load R

### Arguments: ###
if [ $# -lt 2 ]; then
    echo "ERROR: Missing arguments."
    echo "Usage: sbatch prep_sQTL.sh \"Cerebellum\" \"Male\""
    exit 1
fi

TISSUE="$1"
STRAT_SEX="$2"
TISSUE_DIR=$(echo "$TISSUE" | tr ' ' '_')

echo "============================================"
echo "  sQTL Prep for tissue : $TISSUE"
echo "  Stratification       : $STRAT_SEX"
echo "  $(date)"
echo "============================================"

### Paths: ###
BASE=/home/zw529/donglab/data/target_ALS
PLINK=$BASE/QTL/plink
SPLICING_RAW=$BASE/QTL/splicing_matrix.txt 
SPLICING_LOC_SRC=$BASE/QTL/splicing_events_hg38.bed 

RAW=$PLINK/joint_autosomes_matrixEQTL.raw
META=$BASE/QTL/expression_sample_metadata.csv
WGS_MAP=$BASE/QTL/wgs_samples_for_vcf_merge.csv
COV=$BASE/QTL/covariates.tsv
BIM=$PLINK/joint_autosomes_filtered_bed.bim

OUTDIR=$BASE/$TISSUE_DIR/sQTL/$STRAT_SEX
mkdir -p $OUTDIR

# Temp Files
TMP_META=$OUTDIR/tmp_meta.txt
TMP_WGS=$OUTDIR/tmp_wgs_map.txt
TMP_MATCHED=$OUTDIR/tmp_matched.txt
TMP_HRA=$OUTDIR/tmp_HRA_ids.txt
TMP_HRA_UNDER=$OUTDIR/tmp_HRA_ids_underscore.txt
TMP_HDA=$OUTDIR/tmp_HDA_ids.txt
CLEAN_COV=$OUTDIR/covariates_stratified.tsv

##############################################
# STEP 0: Generate BED from Matrix IDs
##############################################
echo "[0] Generating BED file from matrix IDs..."
python3 << 'EOF'
import sys

# Read splicing matrix and extract coordinates
try:
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    opened = False
    
    for enc in encodings:
        try:
            with open(sys.argv[1] if len(sys.argv) > 1 else "/dev/null", 'r', encoding=enc, errors='replace') as f:
                header = f.readline()  # Skip header
                with open(sys.argv[2] if len(sys.argv) > 2 else "/dev/null", 'w') as out:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        # Parse ID format: chr:start-end
                        parts = line.split('\t')
                        if len(parts) < 1:
                            continue
                        event_id = parts[0]
                        
                        # Handle encoding issues by filtering non-ASCII
                        event_id = ''.join(c for c in event_id if ord(c) < 128)
                        
                        try:
                            # Expected format: chr:start-end
                            if ':' in event_id and '-' in event_id:
                                chr_part, pos_part = event_id.rsplit(':', 1)
                                start, end = pos_part.split('-')
                                out.write(f"{chr_part}\t{start}\t{end}\t{event_id}\n")
                        except (ValueError, IndexError):
                            pass
            opened = True
            break
        except (UnicodeDecodeError, FileNotFoundError):
            continue
    
    if not opened:
        raise ValueError("Could not read splicing matrix with any encoding")
        
except Exception as e:
    print(f"ERROR in STEP 0: {e}", file=sys.stderr)
    sys.exit(1)

print("STEP 0: BED file generated successfully", file=sys.stderr)
EOF "$SPLICING_RAW" "$SPLICING_LOC_SRC"

##############################################
# STEP 1: Build sample mapping
##############################################
echo "[1] Building sample mapping..."
python3 << 'EOF'
import pandas as pd
import sys

def try_read_csv(path, **kwargs):
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    for enc in encodings:
        try:
            return pd.read_csv(path, encoding=enc, on_bad_lines='skip', **kwargs), enc
        except UnicodeDecodeError:
            continue
    return None, None

try:
    # Read covariates
    df_cov, enc = try_read_csv(sys.argv[1], sep='\t')
    if df_cov is None:
        raise ValueError("Could not read covariates file")
    print(f"Successfully read covariates with encoding: {enc}", file=sys.stderr)
    
    # Filter out bad samples
    bad_samples = df_cov[df_cov['rin'].isna() & (df_cov['rna_skew'] > 1.0)]['externalsampleid'].tolist()
    df_cov = df_cov[~df_cov['externalsampleid'].isin(bad_samples)]
    
    # Stratify by sex
    df_cov = df_cov[df_cov['sex'] == sys.argv[2]]
    
    # Collapse subject group to binary
    def collapse(g):
        g_str = str(g)
        if 'ALS' in g_str:
            return 1
        if 'Non-Neurological' in g_str:
            return 0
        return None
    
    df_cov['is_als'] = df_cov['subject_group'].apply(collapse)
    df_cov = df_cov[df_cov['is_als'].notna()]
    
    # Save cleaned covariates
    df_cov.to_csv(sys.argv[3], sep='\t', index=False)
    
    # Get list of samples to keep
    keep_samples = df_cov['externalsampleid'].tolist()
    with open(sys.argv[4], 'w') as f:
        for s in keep_samples:
            f.write(str(s).strip() + "\n")
    
    print(f"Filtered to {len(keep_samples)} samples after stratification", file=sys.stderr)
    
    # Read and filter metadata
    df_meta, enc = try_read_csv(sys.argv[5])
    if df_meta is None:
        raise ValueError("Could not read metadata file")
    print(f"Successfully read metadata with encoding: {enc}", file=sys.stderr)
    
    # Read WGS mapping
    df_wgs, enc = try_read_csv(sys.argv[6])
    if df_wgs is None:
        raise ValueError("Could not read WGS map file")
    print(f"Successfully read WGS map with encoding: {enc}", file=sys.stderr)
    
    # Filter metadata by tissue
    sample_col = df_meta.columns[0]  # externalsampleid
    subject_col = df_meta.columns[1]  # externalsubjectid
    tissue_col = df_meta.columns[2]  # tissue
    
    df_meta_filt = df_meta[df_meta[tissue_col] == sys.argv[7]]
    df_meta_filt = df_meta_filt[df_meta_filt[sample_col].isin(keep_samples)]
    
    # Write tmp_meta: sample_id,subject_id
    with open(sys.argv[8], 'w') as f:
        for _, row in df_meta_filt.iterrows():
            f.write(f"{row[sample_col]},{row[subject_col]}\n")
    
    # Write tmp_wgs: subject_id,wgs_id
    wgs_subject_col = df_wgs.columns[0]
    wgs_id_col = df_wgs.columns[1]
    with open(sys.argv[9], 'w') as f:
        for _, row in df_wgs.iterrows():
            f.write(f"{row[wgs_subject_col]},{row[wgs_id_col]}\n")
    
    print("Metadata and WGS mapping written successfully", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 1: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF "$COV" "$STRAT_SEX" "$CLEAN_COV" "$OUTDIR/keep_samples.txt" "$META" "$WGS_MAP" "$TISSUE" "$TMP_META" "$TMP_WGS"

# Join metadata with WGS mapping
python3 << 'EOF'
import pandas as pd
import sys

try:
    meta = pd.read_csv(sys.argv[1], header=None, names=['sample_id', 'subject_id'])
    wgs = pd.read_csv(sys.argv[2], header=None, names=['subject_id', 'wgs_id'])
    
    # Merge on subject_id
    merged = meta.merge(wgs, on='subject_id', how='inner')
    
    # Write matched file
    with open(sys.argv[3], 'w') as f:
        for _, row in merged.iterrows():
            f.write(f"{row['subject_id']},{row['sample_id']},{row['wgs_id']}\n")
    
    # Write HRA IDs
    with open(sys.argv[4], 'w') as f:
        for sample_id in merged['sample_id']:
            f.write(f"{sample_id}\n")
    
    # Write HRA underscore version
    with open(sys.argv[5], 'w') as f:
        for sample_id in merged['sample_id']:
            f.write(f"{sample_id.replace('-', '_')}\n")
    
    # Write HDA IDs
    with open(sys.argv[6], 'w') as f:
        for wgs_id in merged['wgs_id']:
            f.write(f"{wgs_id}\n")
    
    print(f"Matched {len(merged)} samples across metadata and WGS", file=sys.stderr)

except Exception as e:
    print(f"ERROR joining metadata: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF "$TMP_META" "$TMP_WGS" "$TMP_MATCHED" "$TMP_HRA" "$TMP_HRA_UNDER" "$TMP_HDA"

##############################################
# STEP 2: Subset Splicing Matrix
##############################################
echo "[2] Processing Splicing Matrix..."

cat > "$OUTDIR/step2.R" << 'REOF'
library(data.table)

tryCatch({
    tmp_hra_under <- Sys.getenv("TMP_HRA_UNDER")
    splicing_raw <- Sys.getenv("SPLICING_RAW")
    outdir <- Sys.getenv("OUTDIR")
    tissue_dir <- Sys.getenv("TISSUE_DIR")
    
    keep_ids <- fread(tmp_hra_under, header=FALSE)[, V1]
    
    print(paste("Reading splicing matrix from:", splicing_raw))
    full_data <- fread(splicing_raw, header=TRUE, fill=TRUE)
    
    print(paste("Found", nrow(full_data), "rows and", ncol(full_data), "columns"))
    
    # Get column names and intersect with keep_ids
    col_names <- names(full_data)
    cols_to_keep <- c(col_names[1], intersect(col_names[-1], keep_ids))
    
    print(paste("Keeping", length(cols_to_keep), "columns"))
    
    # Subset columns
    splicing <- full_data[, ..cols_to_keep]
    rm(full_data)
    gc()
    
    # Remove low-variance rows
    data_cols <- cols_to_keep[-1]
    row_vars <- apply(splicing[, ..data_cols], 1, function(x) var(x, na.rm=TRUE))
    splicing <- splicing[which(row_vars > 1e-10), ]
    
    print(paste("After variance filtering:", nrow(splicing), "rows"))
    
    # Write output
    output_file <- paste0(outdir, "/splicing_", tissue_dir, ".txt")
    fwrite(splicing, output_file, sep="\t", quote=FALSE)
    print("Splicing matrix written successfully")
    
}, error = function(e) {
    print(paste("ERROR in STEP 2:", e$message))
    quit(status=1)
})
REOF

TMP_HRA_UNDER="$TMP_HRA_UNDER" \
SPLICING_RAW="$SPLICING_RAW" \
OUTDIR="$OUTDIR" \
TISSUE_DIR="$TISSUE_DIR" \
Rscript "$OUTDIR/step2.R"

##############################################
# STEP 3: Subset + Transpose SNPs
##############################################
echo "[3] Subsetting SNPs..."
python3 << 'EOF'
import re
import numpy as np
import sys

try:
    # Read HDA (WGS) sample IDs
    with open(sys.argv[1], 'r') as f:
        keep_ids = set(line.strip() for line in f if line.strip())
    
    print(f"Looking for {len(keep_ids)} WGS samples", file=sys.stderr)
    
    # Read BIM file to get SNP names and positions
    chrpos_names = []
    with open(sys.argv[2], 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                chrpos_names.append(f"chr{parts[0]}:{parts[3]}")
    
    print(f"Found {len(chrpos_names)} SNPs in BIM file", file=sys.stderr)
    
    # Read RAW file with encoding fallback
    sample_ids = []
    rows = []
    
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    opened = False
    
    for enc in encodings:
        try:
            with open(sys.argv[3], 'r', encoding=enc, errors='replace') as f:
                header = f.readline().split()
                n_snps = len(chrpos_names)
                
                for line_num, line in enumerate(f, 1):
                    fields = line.split()
                    if len(fields) < 6 + n_snps:
                        continue
                    
                    # Extract individual ID and remove batch suffix
                    iid = fields[1]
                    iid = re.sub(r'-b\d+$', '', iid)
                    
                    if iid in keep_ids:
                        sample_ids.append(iid)
                        
                        # Extract genotypes
                        geno_fields = fields[6:6+n_snps]
                        geno_row = []
                        for v in geno_fields:
                            try:
                                geno_row.append(float(v))
                            except (ValueError, TypeError):
                                geno_row.append(np.nan)
                        
                        rows.append(geno_row)
            
            opened = True
            print(f"Successfully read RAW file with encoding: {enc}", file=sys.stderr)
            break
        except UnicodeDecodeError:
            continue
    
    if not opened:
        raise ValueError("Could not read RAW file with any encoding")
    
    if not rows:
        print("ERROR: No matching genomic samples found.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found genotypes for {len(sample_ids)} samples", file=sys.stderr)
    
    # Convert to numpy array
    geno = np.array(rows, dtype=float)
    
    # Calculate MAF
    af = np.nanmean(geno, axis=0) / 2.0
    maf = np.minimum(af, 1 - af)
    keep_idx = np.where(np.nan_to_num(maf) >= 0.05)[0]
    
    print(f"Keeping {len(keep_idx)} SNPs with MAF >= 0.05", file=sys.stderr)
    
    # Write output
    with open(sys.argv[4], 'w') as out:
        out.write('snpid\t' + '\t'.join(sample_ids) + '\n')
        for idx in keep_idx:
            vals = []
            for v in geno[:, idx]:
                if np.isnan(v):
                    vals.append('NA')
                else:
                    vals.append(str(int(v)))
            out.write(chrpos_names[idx] + '\t' + '\t'.join(vals) + '\n')
    
    print("SNP file written successfully", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 3: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF "$TMP_HDA" "$BIM" "$RAW" "$OUTDIR/snp_${TISSUE_DIR}.txt"

##############################################
# STEP 4: Subset Covariates
##############################################
echo "[4] Subsetting Covariates..."
python3 << 'EOF'
import pandas as pd
import sys

def try_read_csv(path, **kwargs):
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    for enc in encodings:
        try:
            return pd.read_csv(path, encoding=enc, on_bad_lines='skip', **kwargs), enc
        except UnicodeDecodeError:
            continue
    return None, None

try:
    # Read HRA sample IDs
    with open(sys.argv[1], 'r') as f:
        keep_ids = set(line.strip() for line in f if line.strip())
    
    # Read cleaned covariates
    df_cov, enc = try_read_csv(sys.argv[2], sep='\t')
    if df_cov is None:
        raise ValueError("Could not read covariates file")
    
    # Filter to keep_ids
    df_cov = df_cov[df_cov['externalsampleid'].isin(keep_ids)]
    
    # Save temporary version
    df_cov.to_csv(sys.argv[3], sep='\t', index=False)
    print(f"Subset covariates to {len(df_cov)} samples", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 4: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF "$TMP_HRA" "$CLEAN_COV" "$OUTDIR/covariates_temp.txt"

##############################################
# STEP 5: ID Translation & Final Alignment
##############################################
echo "[5] Final Alignment..."
python3 << 'EOF'
import pandas as pd
import sys

def try_read_csv(path, **kwargs):
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    for enc in encodings:
        try:
            return pd.read_csv(path, encoding=enc, on_bad_lines='skip', **kwargs), enc
        except UnicodeDecodeError:
            continue
    return None, None

try:
    # Read mapping file
    mapping = pd.read_csv(sys.argv[1], header=None, names=['sub', 'hra', 'hda'])
    hda_to_hra = dict(zip(mapping['hda'], mapping['hra']))
    
    # Read SNP file to get sample order
    with open(sys.argv[2], 'r') as f:
        snp_header = f.readline().strip().split('\t')
        snp_hda_samples = snp_header[1:]
    
    print(f"SNP file has {len(snp_hda_samples)} samples", file=sys.stderr)
    
    # Convert HDA IDs to HRA IDs
    target_hra_samples = []
    for hda_id in snp_hda_samples:
        if hda_id in hda_to_hra:
            target_hra_samples.append(hda_to_hra[hda_id])
        else:
            print(f"WARNING: Could not find mapping for {hda_id}", file=sys.stderr)
            target_hra_samples.append(hda_id)
    
    # Read splicing matrix
    df_splicing, enc = try_read_csv(sys.argv[3], sep='\t', index_col=0)
    if df_splicing is None:
        raise ValueError("Could not read splicing file")
    
    # Replace underscores in HRA names
    target_hra_samples_underscore = [s.replace('-', '_') for s in target_hra_samples]
    
    # Subset and reorder splicing matrix
    cols_available = [c for c in target_hra_samples_underscore if c in df_splicing.columns]
    if not cols_available:
        print("ERROR: No matching samples found in splicing matrix", file=sys.stderr)
        sys.exit(1)
    
    df_splicing_aligned = df_splicing[cols_available]
    df_splicing_aligned.columns = snp_hda_samples[:len(cols_available)]
    df_splicing_aligned.to_csv(sys.argv[3], sep='\t')
    
    print(f"Aligned splicing matrix to {len(cols_available)} samples", file=sys.stderr)
    
    # Align covariates
    df_cov, enc = try_read_csv(sys.argv[4], sep='\t', index_col=0)
    target_hra_cov = [s.replace('_', '-') for s in cols_available]
    
    available_cov = [c for c in target_hra_cov if c in df_cov.index]
    df_cov_aligned = df_cov.loc[available_cov].T
    df_cov_aligned.columns = snp_hda_samples[:len(available_cov)]
    df_cov_aligned.to_csv(sys.argv[5], sep='\t')
    
    print(f"Aligned covariates to {len(available_cov)} samples", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 5: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF "$TMP_MATCHED" "$OUTDIR/snp_${TISSUE_DIR}.txt" "$OUTDIR/splicing_${TISSUE_DIR}.txt" "$OUTDIR/covariates_temp.txt" "$OUTDIR/covariates_${TISSUE_DIR}.txt"

##############################################
# STEP 6: Encode Covariates
##############################################
echo "[6] Encoding Covariates..."
python3 << 'EOF'
import csv
import sys

try:
    DROP = {'externalsampleid', 'externalsubjectid', 'tissue', 'subject_group', 'sex'}
    NUMERIC = {'age_at_death', 'rin', 'rna_skew', 'is_als'}
    
    out_rows = []
    header = None
    
    # Try multiple encodings
    encodings = ['utf-8', 'latin1', 'iso-8859-1', 'cp1252']
    opened = False
    
    for enc in encodings:
        try:
            with open(sys.argv[1], 'r', encoding=enc, errors='replace') as f:
                reader = csv.reader(f, delimiter='\t')
                header_row = next(reader)
                header = header_row[1:]
                
                for row in reader:
                    if len(row) < 2:
                        continue
                    
                    cov_name = row[0].lower()
                    
                    if cov_name in DROP:
                        continue
                    
                    if cov_name in NUMERIC:
                        vals = []
                        for v in row[1:]:
                            try:
                                if v and v != 'NA' and v != 'nan':
                                    vals.append(str(float(v)))
                                else:
                                    vals.append('NA')
                            except (ValueError, TypeError):
                                vals.append('NA')
                        out_rows.append((row[0], vals))
            
            opened = True
            break
        except UnicodeDecodeError:
            continue
    
    if not opened:
        raise ValueError("Could not read covariate file with any encoding")
    
    # Write encoded output
    with open(sys.argv[2], 'w') as f:
        f.write('covariate\t' + '\t'.join(header) + '\n')
        for name, vals in out_rows:
            f.write(name + '\t' + '\t'.join(vals) + '\n')
    
    print(f"Encoded {len(out_rows)} covariates", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 6: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
EOF "$OUTDIR/covariates_${TISSUE_DIR}.txt" "$OUTDIR/covariates_${TISSUE_DIR}_encoded.txt"

##############################################
# STEP 7: SNP location
##############################################
echo "[7] Generating SNP location file..."
python3 << 'EOF'
import sys

try:
    with open(sys.argv[2], 'w') as out:
        out.write("snpid\tchr\tpos\n")
        
        with open(sys.argv[1], 'r', errors='replace') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 4:
                    chr_num = parts[0]
                    pos = parts[3]
                    snpid = f"chr{chr_num}:{pos}"
                    out.write(f"{snpid}\tchr{chr_num}\t{pos}\n")
    
    print("SNP location file written", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 7: {e}", file=sys.stderr)
    sys.exit(1)
EOF "$BIM" "$OUTDIR/snp_location.txt"

##############################################
# STEP 8: Splicing location
##############################################
echo "[8] Generating splicing location file..."
python3 << 'EOF'
import sys

try:
    with open(sys.argv[2], 'w') as out:
        out.write("geneid\tchr\tleft\tright\n")
        
        with open(sys.argv[1], 'r', errors='replace') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chr_part = parts[0]
                    left = parts[1]
                    right = parts[2]
                    geneid = parts[3]
                    out.write(f"{geneid}\t{chr_part}\t{left}\t{right}\n")
    
    print("Splicing location file written", file=sys.stderr)

except Exception as e:
    print(f"ERROR in STEP 8: {e}", file=sys.stderr)
    sys.exit(1)
EOF "$SPLICING_LOC_SRC" "$OUTDIR/splicing_location.txt"

##############################################
# STEP 9: Cleanup temporary files
##############################################
echo "[9] Cleaning up temporary files..."
rm -f "$TMP_META" "$TMP_WGS" "$TMP_MATCHED" "$TMP_HRA" "$TMP_HRA_UNDER" "$TMP_HDA" "$OUTDIR/covariates_temp.txt" "$OUTDIR/step2.R"

##############################################
# STEP 10: Final Summary
##############################################
echo "============================================"
echo "  sQTL Prep complete for $TISSUE ($STRAT_SEX)"
echo "  Output directory: $OUTDIR"
echo "  $(date)"
echo "============================================"

# List output files
echo "Output files:"
ls -lh "$OUTDIR"/*.txt 2>/dev/null || echo "No output files found"
