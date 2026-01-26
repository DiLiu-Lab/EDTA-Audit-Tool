#!/bin/bash

# ==============================================================================
# Script: generate_validation_set.sh
# Purpose: Create precise synthetic datasets for EDTA-Audit paper validation
# ==============================================================================

ROOT_DIR="./validation_dataset"
rm -rf "$ROOT_DIR"  # Clean start
mkdir -p "$ROOT_DIR"

# Function to create a realistic EDTA output folder
create_mock_sample() {
    local case_id=$1     # e.g., Case1_Rice
    local te_val=$2      # e.g., 45.5
    local ltr_size=$3    # e.g., 5M
    local line_size=$4   # e.g., 50K
    local heli_size=$5   # e.g., 100K (0 for missing)
    local sine_size=$6   # e.g., 10K (0 for missing)
    
    # 1. Setup Directory Structure
    # Structure: Sample/Sample.mod.EDTA.raw/files...
    local sample_dir="$ROOT_DIR/$case_id"
    local raw_dir="$sample_dir/$case_id.mod.EDTA.raw"
    mkdir -p "$raw_dir"
    
    # 2. Simulate Completion Indicators
    touch "$sample_dir/run.log"
    echo "Total TE: $te_val %" > "$sample_dir/$case_id.TEanno.sum"
    
    # 3. Simulate Library Files (Sizes using truncate)
    # Always include TIR (assume 100K) unless specifically testing it
    truncate -s 100K "$raw_dir/$case_id.TIR.intact.raw.fa"
    
    [ "$ltr_size" != "0" ] && truncate -s $ltr_size "$raw_dir/$case_id.LTR.raw.fa"
    [ "$line_size" != "0" ] && truncate -s $line_size "$raw_dir/$case_id.LINE.raw.fa"
    [ "$heli_size" != "0" ] && truncate -s $heli_size "$raw_dir/$case_id.Helitron.intact.raw.fa"
    [ "$sine_size" != "0" ] && truncate -s $sine_size "$raw_dir/$case_id.SINE.raw.fa"
    
    echo "Created $case_id (TE: $te_val%)"
}

echo ">>> Generating Validation Dataset..."

# --- Group A: PLANTS (Target: >20% TE, >500KB LTR) ---
# Case 1: Perfect Plant (Control)
create_mock_sample "C1_Plant_Pass" "45.0" "5M" "50K" "1M" "5K"
# Case 2: Failed Plant (Low TE - Assembly Collapse)
create_mock_sample "C2_Plant_Fail_TE" "15.0" "5M" "50K" "1M" "5K"
# Case 3: Failed Plant (Small LTR - Bad Annotation)
create_mock_sample "C3_Plant_Fail_LTR" "35.0" "10K" "50K" "1M" "5K"

# --- Group B: ANIMALS (Target: >5% TE, >50KB LINE, Allow Small LTR) ---
# Case 4: Perfect Animal (Fly/Bird) - Note: LTR is small (20K) which would fail in Plant mode!
create_mock_sample "C4_Animal_Pass" "12.0" "20K" "200K" "10K" "5K"
# Case 5: Failed Animal (Missing LINE) - Mammals/Animals need LINEs
create_mock_sample "C5_Animal_Fail_LINE" "15.0" "20K" "0" "10K" "5K"
# Case 6: Low TE Animal (Fail)
create_mock_sample "C6_Animal_Fail_TE" "2.5" "20K" "100K" "10K" "5K"

# --- Group C: FUNGI (Target: >1% TE, Allow Missing Heli/SINE) ---
# Case 7: Perfect Yeast (Minimalist) - Missing Heli & SINE
create_mock_sample "C7_Fungi_Pass_Minimal" "3.5" "10K" "5K" "0" "0"
# Case 8: Failed Fungi (Empty/Too Low)
create_mock_sample "C8_Fungi_Fail_TE" "0.5" "5K" "1K" "0" "0"
# Case 9: Cross-Check (Fungi Logic) - Normal Fungi
create_mock_sample "C9_Fungi_Pass_Std" "5.0" "20K" "10K" "2K" "1K"

echo ">>> Done. Data located in $ROOT_DIR"
