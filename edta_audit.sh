#!/bin/bash

# ============================================================
# Script Name: edta_audit.sh
# Version: 1.0.0
# Description: A comprehensive QC audit tool for EDTA results.
#              Supports flat/nested directory structures and 
#              species-specific biological thresholds.
# Author: [Di Liu/dogdogdoghead]
# ============================================================

# --- Configuration ---
INPUT_DIR=""
FAIL_FILE=""
PASS_FILE=""
QUIET_MODE=0
SPECIES_TYPE="other"
ORGANIZE_MODE=0
MOVE_PASS_DIR=""
MOVE_FAIL_DIR=""

# UI Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'

# Help Function
usage() {
    echo -e "${CYAN}Usage: $0 -d <dir> -t <type> [options]${NC}" >&2
    echo -e "Required:"
    echo -e "  -d <path>  Path to EDTA results directory (Use '.' for current dir)"
    echo -e "  -t <type>  Species: plant, animal, fungi, other"
    echo -e "Options:"
    echo -e "  -O         Organize directory (Aggressively merge scattered files)"
    echo -e "  -m <path>  Move PASSED samples to <path>-edta_audit"
    echo -e "  -x <path>  Move FAILED samples to <path>-edta_audit"
    echo -e "  -o <file>  Output FAILED list (default: edta_failed_list.txt)"
    echo -e "  -p <file>  Output PASSED list (default: edta_passed_list.txt)"
    exit 1
}

# Logger
log() { if [ $QUIET_MODE -eq 0 ]; then echo -e "$@" >&2; fi; }

# Parse Arguments
while getopts "d:t:o:p:m:x:Oqh" opt; do
    case "$opt" in
        d) INPUT_DIR=$OPTARG ;;
        t) SPECIES_TYPE=$(echo "$OPTARG" | tr '[:upper:]' '[:lower:]') ;;
        o) FAIL_FILE=$OPTARG ;;
        p) PASS_FILE=$OPTARG ;;
        m) MOVE_PASS_DIR="${OPTARG}-edta_audit" ;; 
        x) MOVE_FAIL_DIR="${OPTARG}-edta_audit" ;; 
        O) ORGANIZE_MODE=1 ;;
        q) QUIET_MODE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [ -z "$INPUT_DIR" ]; then echo -e "${RED}Error: -d is required.${NC}" >&2; usage; fi
if [ -z "$FAIL_FILE" ]; then FAIL_FILE="edta_failed_list.txt"; fi
if [ -z "$PASS_FILE" ]; then PASS_FILE="edta_passed_list.txt"; fi

# --- Feature: Smart Merge (Using Move) ---
merge_move() {
    local src="$1"; local dest_parent="$2"; local dest_path="$dest_parent/$(basename "$src")"
    if [ -f "$src" ]; then 
        mv -f "$src" "$dest_parent/" 2>/dev/null
    elif [ -d "$src" ]; then
        if [ -d "$dest_path" ]; then
            find "$src" -mindepth 1 -maxdepth 1 -exec mv -f {} "$dest_path/" \;
            rmdir "$src" 2>/dev/null
        else
            mv -f "$src" "$dest_parent/"
        fi
    fi
}

organize_directory() {
    log "${YELLOW}>>> Organizing Directory (Smart Move)...${NC}"
    find "$1" -maxdepth 1 \( -name "GC*" -o -name "*.EDTA.*" -o -name "*.MAKER.*" -o -name "*.mod" \) | while read item; do
        basename_item=$(basename "$item")
        if [[ "$basename_item" == *"-edta_audit" ]]; then continue; fi
        parent_guess="${basename_item%%.fna*}" 
        if [ "$parent_guess" == "$basename_item" ]; then parent_guess="${basename_item%%.mod*}"; fi
        full_folder_path="$1/${parent_guess}"
        if [ -d "$full_folder_path" ]; then
             if [ "$item" != "$full_folder_path" ]; then merge_move "$item" "$full_folder_path"; fi
        fi
    done
    log "${GREEN}>>> Organization complete.${NC}"
    echo "" >&2
}
if [ $ORGANIZE_MODE -eq 1 ]; then organize_directory "$INPUT_DIR"; fi

# --- 🧠 Biological Thresholds Setup ---

# Initialize defaults (General Mode)
LIMIT_LTR=10240      # 10KB
LIMIT_TIR=10240      # 10KB
LIMIT_HELI=5120      # 5KB
LIMIT_LINE=5120      # 5KB
LIMIT_SINE=500       # 500B
MIN_TE_PERCENT=3.0
ALLOW_MISSING_SINE=0
ALLOW_MISSING_HELI=0

case "$SPECIES_TYPE" in
    plant)
        # === Plant Mode ===
        LIMIT_LTR=512000    # 500 KB
        LIMIT_TIR=102400    # 100 KB
        LIMIT_HELI=51200    # 50 KB
        LIMIT_LINE=30720    # 30 KB (Lower in plants)
        LIMIT_SINE=500      # 500 B
        MIN_TE_PERCENT=20.0
        TYPE_LABEL="🌱 Plant Mode (Strict)"
        ;;
    animal)
        # === Animal Mode ===
        LIMIT_LTR=10240     # 10 KB
        LIMIT_TIR=10240     # 10 KB
        LIMIT_HELI=5120     # 5 KB
        LIMIT_LINE=51200    # 50 KB (Higher in animals)
        LIMIT_SINE=500      # 500 B
        MIN_TE_PERCENT=5.0
        TYPE_LABEL="🐅 Animal Mode (Standard)"
        ;;
    fungi)
        # === Fungi Mode ===
        LIMIT_LTR=5120      # 5 KB
        LIMIT_TIR=5120      # 5 KB
        LIMIT_LINE=1024     # 1 KB
        LIMIT_HELI=100      # 100 B
        LIMIT_SINE=100      # 100 B
        MIN_TE_PERCENT=1.0
        ALLOW_MISSING_SINE=1
        ALLOW_MISSING_HELI=1
        TYPE_LABEL="🍄 Fungi Mode (Permissive)"
        ;;
    *)
        TYPE_LABEL="❓ General Mode"
        ;;
esac

# --- Init ---
TEMP_FAIL_LIST="/tmp/edta_fail_$$"; > "$TEMP_FAIL_LIST"
TEMP_PASS_LIST="/tmp/edta_pass_$$"; > "$TEMP_PASS_LIST"

if [ -n "$MOVE_PASS_DIR" ]; then mkdir -p "$MOVE_PASS_DIR"; fi
if [ -n "$MOVE_FAIL_DIR" ]; then mkdir -p "$MOVE_FAIL_DIR"; fi
CNT_TOTAL=0; CNT_FIX=0; CNT_PASS=0; CNT_SKIP=0

if [ $QUIET_MODE -eq 0 ]; then
    log "${CYAN}================================================================${NC}"
    log "${CYAN}          EDTA QC Audit Tool v1.0 (Species Aware)               ${NC}"
    log "${CYAN}================================================================${NC}"
    log "Standard    : $TYPE_LABEL"
    log "Safety      : Skipping folders ending in '-edta_audit'"
    log "${CYAN}----------------------------------------------------------------${NC}"
fi

CURRENT_FAIL_REASON=""

check_asset() {
    local raw_dir=$1; local pattern=$2; local label=$3; local thresh=$4; local allow=$5
    local target=$(find "$raw_dir" -name "$pattern" | head -n 1)
    if [ -s "$target" ]; then
        local fsize=$(stat -c%s "$target")
        local fsize_h=$(numfmt --to=iec $fsize)
        if [ $fsize -lt $thresh ]; then 
            [ $QUIET_MODE -eq 0 ] && echo -ne "${RED}${label}:${fsize_h}(Small)${NC} " >&2
            CURRENT_FAIL_REASON="${CURRENT_FAIL_REASON}${label}:Small;"
            return 1 
        else 
            [ $QUIET_MODE -eq 0 ] && echo -ne "${GREEN}${label}:${fsize_h}${NC} " >&2; return 0; 
        fi
    else 
        if [ "$allow" == "1" ]; then
            [ $QUIET_MODE -eq 0 ] && echo -ne "${YELLOW}${label}:Miss(OK)${NC} " >&2; return 0
        else
            [ $QUIET_MODE -eq 0 ] && echo -ne "${RED}${label}:Miss${NC} " >&2
            CURRENT_FAIL_REASON="${CURRENT_FAIL_REASON}${label}:Missing;"
            return 1
        fi
    fi
}

# --- Main Scan Loop ---
if [ ! -d "$INPUT_DIR" ]; then
    echo -e "${RED}Error: Directory '$INPUT_DIR' does not exist.${NC}" >&2
    exit 1
fi

find "$INPUT_DIR" -maxdepth 1 -mindepth 1 -type d > /tmp/edta_dirs_$$

while read work_dir; do
    dir_name=$(basename "$work_dir")
    if [[ "$dir_name" == *".EDTA."* || "$dir_name" == *".MAKER."* ]]; then continue; fi
    if [[ "$dir_name" == *"-edta_audit" ]]; then continue; fi

    HAS_LOG=$(find "$work_dir" -maxdepth 1 -name "run.log" | head -n 1)
    HAS_SUM=$(find "$work_dir" -maxdepth 1 -name "*.TEanno.sum" | head -n 1)
    HAS_RAW=$(find "$work_dir" -maxdepth 1 -type d -name "*.mod.EDTA.raw" | head -n 1)
    if [ -z "$HAS_LOG" ] && [ -z "$HAS_SUM" ] && [ -z "$HAS_RAW" ]; then continue; fi
    
    ((CNT_TOTAL++))
    CURRENT_FAIL_REASON=""
    
    SUM_FILE=$(find "$work_dir" -name "*.TEanno.sum" | awk '{ print length, $0 }' | sort -n -s | cut -d" " -f2- | head -n 1)
    RAW_DIR=$(find "$work_dir" -type d -name "*.mod.EDTA.raw" | awk '{ print length, $0 }' | sort -n -s | cut -d" " -f2- | head -n 1)
    
    if [ -f "$work_dir/run.log" ] && ps aux | grep -v grep | grep -q "$work_dir"; then
        log "[${CNT_TOTAL}] ${BLUE}$dir_name${NC} -> ⏳ Running... (Skipped)"; continue
    fi

    [ $QUIET_MODE -eq 0 ] && echo -ne "[${CNT_TOTAL}] ${YELLOW}$dir_name${NC}\n    " >&2
    NEED_FIX=0

    # 1. Check TE %
    TE_VAL=0
    if [ -n "$SUM_FILE" ] && [ -f "$SUM_FILE" ]; then
        TE_CONTENT=$(grep "Total TE:" "$SUM_FILE" | awk '{print $3}' | sed 's/[()%]//g')
        if [ -z "$TE_CONTENT" ]; then TE_CONTENT=$(grep "^Total" "$SUM_FILE" | grep "%" | awk '{for(i=1;i<=NF;i++) if($i~/%/) print $i}' | sed 's/%//g' | head -n 1); fi
        [ -z "$TE_CONTENT" ] && TE_CONTENT=0
        is_low=$(awk -v v="$TE_CONTENT" -v th="$MIN_TE_PERCENT" 'BEGIN {print (v < th) ? 1 : 0}')
        
        if [ "$is_low" -eq 1 ]; then
            [ $QUIET_MODE -eq 0 ] && echo -ne "📊${RED}${TE_CONTENT}%${NC} | " >&2
            NEED_FIX=1
            CURRENT_FAIL_REASON="${CURRENT_FAIL_REASON}LowTE(${TE_CONTENT}%);"
        else
            [ $QUIET_MODE -eq 0 ] && echo -ne "📊${GREEN}${TE_CONTENT}%${NC} | " >&2
        fi
    else
        [ $QUIET_MODE -eq 0 ] && echo -ne "📊${RED}NoResults${NC} | " >&2
        NEED_FIX=1
        CURRENT_FAIL_REASON="${CURRENT_FAIL_REASON}No_TEanno_Sum;"
    fi

    # 2. Check 5 Libraries
    if [ -n "$RAW_DIR" ] && [ -d "$RAW_DIR" ]; then
        check_asset "$RAW_DIR" "*.LTR.raw.fa" "LTR" $LIMIT_LTR 0 || NEED_FIX=1
        check_asset "$RAW_DIR" "*.TIR.intact.raw.fa" "TIR" $LIMIT_TIR 0 || NEED_FIX=1
        check_asset "$RAW_DIR" "*.LINE.raw.fa" "LINE" $LIMIT_LINE 0 || NEED_FIX=1
        check_asset "$RAW_DIR" "*.Helitron.intact.raw.fa" "Heli" $LIMIT_HELI $ALLOW_MISSING_HELI || NEED_FIX=1
        check_asset "$RAW_DIR" "*.SINE.raw.fa" "SINE" $LIMIT_SINE $ALLOW_MISSING_SINE || NEED_FIX=1
    else
        [ $QUIET_MODE -eq 0 ] && echo -ne "${RED}RawDirMissing${NC}" >&2
        NEED_FIX=1
        CURRENT_FAIL_REASON="${CURRENT_FAIL_REASON}RawDirMissing(Try -O);"
    fi
    
    [ $QUIET_MODE -eq 0 ] && echo "" >&2

    # 3. Handle Results
    if [ $NEED_FIX -eq 1 ]; then
        CURRENT_FAIL_REASON=${CURRENT_FAIL_REASON%;}
        echo -e "$dir_name\t$CURRENT_FAIL_REASON" >> "$TEMP_FAIL_LIST"
        ((CNT_FIX++))
        if [ -n "$MOVE_FAIL_DIR" ]; then mv "$work_dir" "$MOVE_FAIL_DIR/"; fi
    else
        echo "$dir_name" >> "$TEMP_PASS_LIST"
        ((CNT_PASS++))
        if [ -n "$MOVE_PASS_DIR" ]; then mv "$work_dir" "$MOVE_PASS_DIR/"; fi
    fi
done < /tmp/edta_dirs_$$; rm -f /tmp/edta_dirs_$$

# --- Report ---
mv "$TEMP_FAIL_LIST" "$FAIL_FILE"
mv "$TEMP_PASS_LIST" "$PASS_FILE"

if [ $QUIET_MODE -eq 0 ]; then
    log "${CYAN}----------------------------------------------------------------${NC}"
    log "Audit Complete! Scanned: $CNT_TOTAL"
    log "Passed : ${GREEN}$CNT_PASS${NC} -> Saved to $PASS_FILE"
    log "Failed : ${RED}$CNT_FIX${NC} -> Saved to $FAIL_FILE"
    if [ -n "$MOVE_PASS_DIR" ]; then log "Moved Passed to: $MOVE_PASS_DIR"; fi
    if [ -n "$MOVE_FAIL_DIR" ]; then log "Moved Failed to: $MOVE_FAIL_DIR"; fi
fi
