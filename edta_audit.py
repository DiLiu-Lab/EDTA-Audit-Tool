#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EDTA-Audit (normalized) - scan sample directories and classify incomplete outputs as TECH_FAIL

Key change:
- Scan "sample directories" (default pattern: GCA_*_genomic or GCF_*_genomic) under root
- For each sample dir, try to locate expected EDTA outputs; if missing -> tech FAIL and still reported
- If tech FAIL: overall_status=FAIL; bio_status=NA (not evaluated); metrics=NA

Outputs: all.tsv, pass.tsv, suspect.tsv, fail.tsv, summary.txt
"""

import argparse
import csv
import os
import re
import sys
from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Dict, List, Tuple, Optional
from collections import Counter


# -------------------------
# utils
# -------------------------
def now_ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def eprint(msg: str, quiet: bool = False):
    if not quiet:
        print(msg, file=sys.stderr)


def is_nonempty_file(path: str) -> bool:
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def is_dir(path: str) -> bool:
    return os.path.isdir(path)


def safe_int(x: str) -> Optional[int]:
    try:
        return int(x.replace(",", ""))
    except Exception:
        return None


def safe_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except Exception:
        return None


def fmt_num(x: float, nd: int = 3) -> str:
    if x != x:  # NaN
        return "NA"
    if nd == 0:
        return str(int(round(x)))
    return f"{x:.{nd}f}"


def to_opt(x: float) -> Optional[float]:
    return None if (x != x) else x


# -------------------------
# TEanno.sum parser
# -------------------------
_SPLIT_COLS = re.compile(r"\s{2,}")


def parse_teanno_sum_repeat_classes(sum_path: str) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    in_table = False
    current_class: Optional[str] = None

    def norm_class(s: str) -> str:
        s = s.strip().upper()
        if s == "NONTIR":
            return "NONTIR"
        return s

    with open(sum_path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if not in_table:
                if re.search(r"^\s*Class\s+Count\s+bpMasked\s+%masked\s*$", line):
                    in_table = True
                    current_class = None
                continue

            if line.strip().startswith("Repeat Stats"):
                break

            s = line.strip()
            if not s:
                continue

            if set(s) <= set("-=") and len(s) >= 5:
                continue

            cols = _SPLIT_COLS.split(s)
            if len(cols) < 4:
                continue

            name, count_s, bp_s, pct_s = cols[0], cols[1], cols[2], cols[3]

            if count_s == "--" and bp_s == "--" and pct_s == "--":
                current_class = norm_class(name)
                continue

            bp = safe_int(bp_s)
            m = re.match(r"^(\d+(?:\.\d+)?)%?$", pct_s)
            pct = safe_float(m.group(1)) if m else None
            if pct is None:
                continue

            if name.lower() == "total":
                out["TOTAL"] = {"bp": float(bp) if bp is not None else float("nan"), "pct": float(pct)}
                current_class = None
                continue

            if name.lower() == "total interspersed":
                out["TOTAL_INTERSPERSED"] = {"bp": float(bp) if bp is not None else float("nan"), "pct": float(pct)}
                continue

            is_family = raw.startswith(" ") or raw.startswith("\t")
            if is_family and current_class:
                key = f"{current_class}/{name.strip()}"
            else:
                key = name.strip().upper()

            out[key] = {"bp": float(bp) if bp is not None else float("nan"), "pct": float(pct)}

    return out


def get_key(te: Dict[str, Dict[str, float]], key: str) -> Tuple[float, float, bool]:
    if key not in te:
        return float("nan"), float("nan"), False
    v = te[key]
    return float(v.get("pct", float("nan"))), float(v.get("bp", float("nan"))), True


def sum_by_prefix(te: Dict[str, Dict[str, float]], prefix: str) -> Tuple[float, float, bool]:
    pct = 0.0
    bp = 0.0
    bp_ok = False
    present = False
    for k, v in te.items():
        if k.startswith(prefix):
            present = True
            pct += float(v.get("pct", 0.0))
            b = float(v.get("bp", float("nan")))
            if b == b:
                bp += b
                bp_ok = True
    if not present:
        return float("nan"), float("nan"), False
    if not bp_ok:
        bp = float("nan")
    return pct, bp, True


def pick_total(te: Dict[str, Dict[str, float]]) -> Tuple[float, float, str, bool]:
    pct, bp, ok = get_key(te, "TOTAL")
    if ok:
        return pct, bp, "TOTAL", True
    pct, bp, ok = get_key(te, "TOTAL_INTERSPERSED")
    if ok:
        return pct, bp, "TOTAL_INTERSPERSED", True
    return float("nan"), float("nan"), "MISSING", False


# -------------------------
# profiles
# -------------------------
@dataclass
class PlantThresholds:
    total_pass: float = 20.0
    total_suspect: float = 10.0
    ltr_pass: float = 5.0
    ltr_suspect: float = 1.0
    tir_pass: float = 2.0
    tir_suspect: float = 0.5
    heli_warn_low: float = 0.1
    heli_warn_mid: float = 0.5
    ltr_known_hint_low: float = 1.0
    ltr_unknown_pct_hint: float = 10.0
    ltr_unknown_ratio_hint: float = 0.60


@dataclass
class AnimalThresholds:
    total_pass: float = 10.0
    total_suspect: float = 3.0
    min_ltr: float = 0.1
    min_tir: float = 0.1
    min_nonltr: float = 0.2


@dataclass
class FungiThresholds:
    total_pass: float = 5.0
    total_suspect: float = 0.5
    min_ltr: float = 0.1
    min_tir: float = 0.1


def tri_level(val: Optional[float], pass_cut: float, suspect_cut: float) -> str:
    if val is None:
        return "FAIL"
    if val >= pass_cut:
        return "PASS"
    if val >= suspect_cut:
        return "SUSPECT"
    return "FAIL"


def classify_plant(total_te, ltr_total, tir, helitron, line_pct, sine_pct, ltr_known, ltr_unknown, thr: PlantThresholds):
    notes, tags = [], []
    st_total = tri_level(total_te, thr.total_pass, thr.total_suspect)
    st_ltr = tri_level(ltr_total, thr.ltr_pass, thr.ltr_suspect)
    st_tir = tri_level(tir, thr.tir_pass, thr.tir_suspect)

    if st_total == "FAIL": tags.append("HARD_FAIL_TOTAL_TE")
    elif st_total == "SUSPECT": tags.append("HARD_SUSPECT_TOTAL_TE")

    if st_ltr == "FAIL": tags.append("HARD_FAIL_LTR_TOTAL")
    elif st_ltr == "SUSPECT": tags.append("HARD_SUSPECT_LTR_TOTAL")

    if st_tir == "FAIL": tags.append("HARD_FAIL_TIR")
    elif st_tir == "SUSPECT": tags.append("HARD_SUSPECT_TIR")

    if helitron is None:
        tags.append("WARN_HELITRON_NOT_REPORTED")
    else:
        if helitron < thr.heli_warn_low: tags.append("WARN_LOW_HELITRON")
        elif helitron < thr.heli_warn_mid: tags.append("WARN_LOW_HELITRON")

    if line_pct is None: tags.append("WARN_LINE_NOT_REPORTED")
    elif line_pct == 0: tags.append("WARN_LINE_ZERO")

    if sine_pct is None: tags.append("WARN_SINE_NOT_REPORTED")
    elif sine_pct == 0: tags.append("WARN_SINE_ZERO")

    if ltr_known is None:
        tags.append("WARN_LTR_KNOWN_NOT_REPORTED")
    else:
        if ltr_known < thr.ltr_known_hint_low:
            tags.append("HINT_LOW_LTR_KNOWN")

    if (ltr_unknown is not None) and (ltr_total is not None) and ltr_total > 0:
        ratio = ltr_unknown / ltr_total
        if (ltr_unknown >= thr.ltr_unknown_pct_hint) or (ratio >= thr.ltr_unknown_ratio_hint):
            tags.append("HINT_HIGH_LTR_UNKNOWN")

    if "FAIL" in (st_total, st_ltr, st_tir):
        return "FAIL", notes, tags
    if "SUSPECT" in (st_total, st_ltr, st_tir):
        return "SUSPECT", notes, tags
    return "PASS", notes, tags


def classify_animal(total_te, ltr_total, tir, line_pct, sine_pct, thr: AnimalThresholds):
    notes, tags = [], []
    st_total = tri_level(total_te, thr.total_pass, thr.total_suspect)
    if st_total == "FAIL": tags.append("HARD_FAIL_TOTAL_TE")
    elif st_total == "SUSPECT": tags.append("HARD_SUSPECT_TOTAL_TE")

    if ltr_total is None or ltr_total < thr.min_ltr: tags.append("WARN_LOW_LTR")
    if tir is None or tir < thr.min_tir: tags.append("WARN_LOW_TIR")

    nonltr = 0.0
    ok = False
    if line_pct is not None: nonltr += line_pct; ok = True
    if sine_pct is not None: nonltr += sine_pct; ok = True
    if (not ok) or (nonltr < thr.min_nonltr): tags.append("WARN_LOW_NONLTR")

    return st_total, notes, tags


def classify_fungi(total_te, ltr_total, tir, helitron, line_pct, sine_pct, thr: FungiThresholds):
    notes, tags = [], []
    st_total = tri_level(total_te, thr.total_pass, thr.total_suspect)
    if st_total == "FAIL": tags.append("HARD_FAIL_TOTAL_TE")
    elif st_total == "SUSPECT": tags.append("HARD_SUSPECT_TOTAL_TE")

    if ltr_total is None or ltr_total < thr.min_ltr: tags.append("WARN_LOW_LTR")
    if tir is None or tir < thr.min_tir: tags.append("WARN_LOW_TIR")

    if helitron is None: tags.append("WARN_HELITRON_NOT_REPORTED")
    if line_pct is None: tags.append("WARN_LINE_NOT_REPORTED")
    if sine_pct is None: tags.append("WARN_SINE_NOT_REPORTED")

    return st_total, notes, tags


# -------------------------
# tech + sample discovery
# -------------------------
@dataclass
class TechConfig:
    require_teanno_sum: bool = True
    require_teanno_gff3: bool = True
    require_telib_fa: bool = True
    require_anno_dir: bool = True


def extract_genome_id(sample_id: str) -> str:
    m = re.match(r"^(GC[AF]_\d+\.\d+)", sample_id)
    return m.group(1) if m else sample_id


def find_sample_dirs(root: str, recursive: bool, max_depth: int, sample_regex: str) -> List[str]:
    r"""
    Find directories matching sample_regex (default: ^GC[AF]_\d+\.\d+.*_genomic$)
    """
    root = os.path.abspath(root)
    root_depth = root.rstrip(os.sep).count(os.sep)
    rx = re.compile(sample_regex)

    skip_markers = [
        ".EDTA.anno" + os.sep,
        ".EDTA.raw" + os.sep,
        ".EDTA.combine" + os.sep,
        ".EDTA.final" + os.sep,
    ]

    samples = []
    for cur, dirs, files in os.walk(root):
        cur_depth = cur.rstrip(os.sep).count(os.sep) - root_depth
        if (not recursive) and cur_depth > 0:
            continue
        if recursive and max_depth >= 0 and cur_depth > max_depth:
            dirs[:] = []
            continue

        cur_abs = os.path.abspath(cur) + os.sep
        if any(m in cur_abs for m in skip_markers):
            dirs[:] = []
            continue

        base = os.path.basename(cur.rstrip(os.sep))
        if rx.match(base):
            samples.append(os.path.abspath(cur))
            # 不继续往下找（一个样本目录作为一个单位）
            dirs[:] = []

    return sorted(set(samples))


def locate_outputs(sample_dir: str) -> Dict[str, str]:
    """
    Locate expected files inside sample_dir:
      - *.fna.mod.EDTA.TEanno.sum (or any *.EDTA.TEanno.sum)
      - corresponding *.TEanno.gff3, *.TElib.fa
      - *.anno directory
    """
    teanno_sum = ""
    for fn in os.listdir(sample_dir):
        if fn.endswith(".EDTA.TEanno.sum") or fn.endswith("EDTA.TEanno.sum"):
            teanno_sum = os.path.join(sample_dir, fn)
            break

    teanno_gff3 = ""
    telib_fa = ""
    anno_dir = ""

    if teanno_sum:
        base = teanno_sum[:-len("TEanno.sum")]
        teanno_gff3 = base + "TEanno.gff3"
        telib_fa = base + "TElib.fa"
        prefix = teanno_sum.replace(".TEanno.sum", "")
        anno_dir = prefix + ".anno"
    else:
        # 没有 TEanno.sum 时，尝试用“最可能的 prefix”猜测（可用于 missing 统计）
        # 找 *.fna.mod 或 *.fna 作为种子
        seed = ""
        for fn in os.listdir(sample_dir):
            if fn.endswith(".fna.mod"):
                seed = os.path.join(sample_dir, fn)
                break
        if not seed:
            for fn in os.listdir(sample_dir):
                if fn.endswith(".fna"):
                    seed = os.path.join(sample_dir, fn)
                    break
        if seed:
            guess_prefix = seed + ".EDTA."
            teanno_gff3 = guess_prefix + "TEanno.gff3"
            telib_fa = guess_prefix + "TElib.fa"
            anno_dir = guess_prefix[:-1] + "anno"  # -> ".EDTA.anno"
        else:
            # 留空
            pass

    return {
        "sample_dir": sample_dir,
        "teanno_sum": teanno_sum,
        "teanno_gff3": teanno_gff3,
        "telib_fa": telib_fa,
        "anno_dir": anno_dir,
    }


def tech_check(paths: Dict[str, str], cfg: TechConfig) -> Tuple[str, List[str]]:
    missing: List[str] = []

    if cfg.require_teanno_sum:
        if (not paths["teanno_sum"]) or (not is_nonempty_file(paths["teanno_sum"])):
            missing.append("TEanno.sum")

    if cfg.require_teanno_gff3:
        if (not paths["teanno_gff3"]) or (not is_nonempty_file(paths["teanno_gff3"])):
            missing.append("TEanno.gff3")

    if cfg.require_telib_fa:
        if (not paths["telib_fa"]) or (not is_nonempty_file(paths["telib_fa"])):
            missing.append("TElib.fa")

    if cfg.require_anno_dir:
        if (not paths["anno_dir"]) or (not is_dir(paths["anno_dir"])):
            missing.append("EDTA.anno/")

    return ("FAIL", missing) if missing else ("PASS", [])


# -------------------------
# record
# -------------------------
@dataclass
class Record:
    profile: str
    genome_id: str
    sample_id: str
    overall_status: str
    tech_status: str
    bio_status: str

    total_te_pct: float
    total_te_bp: float
    total_te_source: str

    ltr_total_pct: float
    ltr_total_bp: float
    ltr_known_pct: float
    ltr_known_bp: float
    ltr_unknown_pct: float
    ltr_unknown_bp: float

    tir_pct: float
    tir_bp: float

    helitron_pct: float
    helitron_bp: float

    line_present: int
    line_pct: float
    line_bp: float

    sine_present: int
    sine_pct: float
    sine_bp: float

    notes: str
    tags: str
    missing: str
    teanno_sum: str
    sample_dir: str


def classify_one(sample_dir: str,
                 profile: str,
                 plant_thr: PlantThresholds,
                 animal_thr: AnimalThresholds,
                 fungi_thr: FungiThresholds,
                 tech_cfg: TechConfig) -> Record:
    sample_dir = os.path.abspath(sample_dir)
    sample_id = os.path.basename(sample_dir)
    genome_id = extract_genome_id(sample_id)

    paths = locate_outputs(sample_dir)
    tech_status, missing_list = tech_check(paths, tech_cfg)

    notes: List[str] = []
    tags: List[str] = []

    # Tech fail: do not evaluate bio
    if tech_status == "FAIL":
        for m in missing_list:
            tags.append(f"TECH_MISSING_{m}")
        return Record(
            profile=profile,
            genome_id=genome_id,
            sample_id=sample_id,
            overall_status="FAIL",
            tech_status="FAIL",
            bio_status="NA",
            total_te_pct=float("nan"),
            total_te_bp=float("nan"),
            total_te_source="NA",
            ltr_total_pct=float("nan"),
            ltr_total_bp=float("nan"),
            ltr_known_pct=float("nan"),
            ltr_known_bp=float("nan"),
            ltr_unknown_pct=float("nan"),
            ltr_unknown_bp=float("nan"),
            tir_pct=float("nan"),
            tir_bp=float("nan"),
            helitron_pct=float("nan"),
            helitron_bp=float("nan"),
            line_present=0,
            line_pct=float("nan"),
            line_bp=float("nan"),
            sine_present=0,
            sine_pct=float("nan"),
            sine_bp=float("nan"),
            notes="TECH_FAIL: incomplete EDTA output",
            tags="; ".join(tags),
            missing="; ".join(missing_list),
            teanno_sum=paths["teanno_sum"] or "NA",
            sample_dir=sample_dir
        )

    # Parse TEanno.sum
    teanno_sum_path = paths["teanno_sum"]
    try:
        te = parse_teanno_sum_repeat_classes(teanno_sum_path)
    except Exception as e:
        tags.append("TECH_PARSE_ERROR")
        return Record(
            profile=profile,
            genome_id=genome_id,
            sample_id=sample_id,
            overall_status="FAIL",
            tech_status="FAIL",
            bio_status="NA",
            total_te_pct=float("nan"),
            total_te_bp=float("nan"),
            total_te_source="NA",
            ltr_total_pct=float("nan"),
            ltr_total_bp=float("nan"),
            ltr_known_pct=float("nan"),
            ltr_known_bp=float("nan"),
            ltr_unknown_pct=float("nan"),
            ltr_unknown_bp=float("nan"),
            tir_pct=float("nan"),
            tir_bp=float("nan"),
            helitron_pct=float("nan"),
            helitron_bp=float("nan"),
            line_present=0,
            line_pct=float("nan"),
            line_bp=float("nan"),
            sine_present=0,
            sine_pct=float("nan"),
            sine_bp=float("nan"),
            notes=f"TECH_FAIL: TEanno.sum parse error: {e}",
            tags="; ".join(tags),
            missing="TEanno.sum(parse_error)",
            teanno_sum=teanno_sum_path,
            sample_dir=sample_dir
        )

    total_pct, total_bp, total_source, _ = pick_total(te)
    tags.append(f"INFO_TOTAL_SOURCE_{total_source}")

    ltr_total_pct, ltr_total_bp, ltr_present = sum_by_prefix(te, "LTR/")
    copia_pct, copia_bp, copia_present = get_key(te, "LTR/Copia")
    gypsy_pct, gypsy_bp, gypsy_present = get_key(te, "LTR/Gypsy")
    ltr_known_pct = float("nan")
    ltr_known_bp = float("nan")
    if copia_present or gypsy_present:
        ltr_known_pct = 0.0
        ltr_known_bp = 0.0
        bp_ok = False
        if copia_present:
            ltr_known_pct += copia_pct
            if copia_bp == copia_bp:
                ltr_known_bp += copia_bp
                bp_ok = True
        if gypsy_present:
            ltr_known_pct += gypsy_pct
            if gypsy_bp == gypsy_bp:
                ltr_known_bp += gypsy_bp
                bp_ok = True
        if not bp_ok:
            ltr_known_bp = float("nan")

    ltr_unk_pct, ltr_unk_bp, ltr_unk_present = get_key(te, "LTR/unknown")
    if not ltr_unk_present:
        ltr_unk_pct, ltr_unk_bp, ltr_unk_present = get_key(te, "LTR/Unknown")
    if not ltr_unk_present:
        ltr_unk_pct, ltr_unk_bp = float("nan"), float("nan")

    tir_pct, tir_bp, _ = sum_by_prefix(te, "TIR/")

    hel_pct, hel_bp, hel_present = get_key(te, "NONTIR/helitron")
    if not hel_present:
        hel_pct, hel_bp, hel_present = get_key(te, "NONTIR/Helitron")
    if not hel_present:
        hel_pct, hel_bp = float("nan"), float("nan")

    line_pct, line_bp, line_present = sum_by_prefix(te, "LINE/")
    sine_pct, sine_bp, sine_present = sum_by_prefix(te, "SINE/")

    total_o = to_opt(total_pct)
    ltr_total_o = to_opt(ltr_total_pct)
    tir_o = to_opt(tir_pct)
    hel_o = to_opt(hel_pct)
    line_o = to_opt(line_pct)
    sine_o = to_opt(sine_pct)
    ltr_known_o = to_opt(ltr_known_pct)
    ltr_unk_o = to_opt(ltr_unk_pct)

    if profile == "plant":
        bio_status, _, bio_tags = classify_plant(total_o, ltr_total_o, tir_o, hel_o, line_o, sine_o, ltr_known_o, ltr_unk_o, plant_thr)
    elif profile == "animal":
        bio_status, _, bio_tags = classify_animal(total_o, ltr_total_o, tir_o, line_o, sine_o, animal_thr)
    else:
        bio_status, _, bio_tags = classify_fungi(total_o, ltr_total_o, tir_o, hel_o, line_o, sine_o, fungi_thr)

    tags.extend(bio_tags)
    overall_status = bio_status  # tech is PASS here

    return Record(
        profile=profile,
        genome_id=genome_id,
        sample_id=sample_id,
        overall_status=overall_status,
        tech_status="PASS",
        bio_status=bio_status,

        total_te_pct=total_pct,
        total_te_bp=total_bp,
        total_te_source=total_source,

        ltr_total_pct=ltr_total_pct,
        ltr_total_bp=ltr_total_bp,
        ltr_known_pct=ltr_known_pct,
        ltr_known_bp=ltr_known_bp,
        ltr_unknown_pct=ltr_unk_pct,
        ltr_unknown_bp=ltr_unk_bp,

        tir_pct=tir_pct,
        tir_bp=tir_bp,

        helitron_pct=hel_pct,
        helitron_bp=hel_bp,

        line_present=1 if line_present else 0,
        line_pct=line_pct,
        line_bp=line_bp,

        sine_present=1 if sine_present else 0,
        sine_pct=sine_pct,
        sine_bp=sine_bp,

        notes="",
        tags="; ".join([t for t in tags if t]),
        missing="",
        teanno_sum=teanno_sum_path,
        sample_dir=sample_dir
    )


# -------------------------
# writers
# -------------------------
def write_tsv(path: str, rows: List[Record]):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fields = [
        "profile","genome_id","sample_id",
        "overall_status","tech_status","bio_status",
        "total_te_pct","total_te_bp","total_te_source",
        "ltr_total_pct","ltr_total_bp","ltr_known_pct","ltr_known_bp","ltr_unknown_pct","ltr_unknown_bp",
        "tir_pct","tir_bp","helitron_pct","helitron_bp",
        "line_present","line_pct","line_bp","sine_present","sine_pct","sine_bp",
        "notes","tags","missing","teanno_sum","sample_dir"
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            d = asdict(r)
            d["total_te_pct"] = fmt_num(r.total_te_pct, 3)
            d["total_te_bp"] = fmt_num(r.total_te_bp, 0)
            d["ltr_total_pct"] = fmt_num(r.ltr_total_pct, 3)
            d["ltr_total_bp"] = fmt_num(r.ltr_total_bp, 0)
            d["ltr_known_pct"] = fmt_num(r.ltr_known_pct, 3)
            d["ltr_known_bp"] = fmt_num(r.ltr_known_bp, 0)
            d["ltr_unknown_pct"] = fmt_num(r.ltr_unknown_pct, 3)
            d["ltr_unknown_bp"] = fmt_num(r.ltr_unknown_bp, 0)
            d["tir_pct"] = fmt_num(r.tir_pct, 3)
            d["tir_bp"] = fmt_num(r.tir_bp, 0)
            d["helitron_pct"] = fmt_num(r.helitron_pct, 3)
            d["helitron_bp"] = fmt_num(r.helitron_bp, 0)
            d["line_pct"] = fmt_num(r.line_pct, 3)
            d["line_bp"] = fmt_num(r.line_bp, 0)
            d["sine_pct"] = fmt_num(r.sine_pct, 3)
            d["sine_bp"] = fmt_num(r.sine_bp, 0)
            w.writerow(d)


def write_summary(path: str, records: List[Record], profile: str):
    total = len(records)
    c_overall = Counter(r.overall_status for r in records)
    c_tech = Counter(r.tech_status for r in records)
    c_bio = Counter(r.bio_status for r in records)

    fail_tags = Counter()
    suspect_tags = Counter()
    tech_missing = Counter()

    def split_tags(s: str) -> List[str]:
        if not s:
            return []
        return [x.strip() for x in s.split(";") if x.strip()]

    for r in records:
        tags = split_tags(r.tags)
        if r.tech_status == "FAIL":
            for t in tags:
                if t.startswith("TECH_MISSING_"):
                    tech_missing[t.replace("TECH_MISSING_", "")] += 1
        if r.overall_status == "FAIL":
            for t in tags:
                if t.startswith("TECH_") or t.startswith("HARD_FAIL_"):
                    fail_tags[t] += 1
        elif r.overall_status == "SUSPECT":
            for t in tags:
                if t.startswith("HARD_SUSPECT_"):
                    suspect_tags[t] += 1

    def dump_counter(title: str, c: Counter):
        lines = [title]
        if not c:
            lines.append("  (none)")
        else:
            for k, v in c.most_common():
                lines.append(f"  {k}\t{v}")
        return "\n".join(lines)

    with open(path, "w", encoding="utf-8") as f:
        f.write("EDTA-Audit summary\n")
        f.write(f"Generated: {now_ts()}\n")
        f.write(f"Profile: {profile}\n\n")
        f.write("Counts\n")
        f.write(f"  TOTAL\t{total}\n")
        f.write(f"  OVERALL_PASS\t{c_overall.get('PASS', 0)}\n")
        f.write(f"  OVERALL_SUSPECT\t{c_overall.get('SUSPECT', 0)}\n")
        f.write(f"  OVERALL_FAIL\t{c_overall.get('FAIL', 0)}\n")
        f.write(f"  TECH_PASS\t{c_tech.get('PASS', 0)}\n")
        f.write(f"  TECH_FAIL\t{c_tech.get('FAIL', 0)}\n")
        f.write(f"  BIO_PASS\t{c_bio.get('PASS', 0)}\n")
        f.write(f"  BIO_SUSPECT\t{c_bio.get('SUSPECT', 0)}\n")
        f.write(f"  BIO_FAIL\t{c_bio.get('FAIL', 0)}\n")
        f.write(f"  BIO_NA\t{c_bio.get('NA', 0)}\n\n")
        f.write(dump_counter("FAIL reason tags (TECH_* + HARD_FAIL_*)", fail_tags) + "\n\n")
        f.write(dump_counter("SUSPECT reason tags (HARD_SUSPECT_*)", suspect_tags) + "\n\n")
        f.write(dump_counter("TECH missing breakdown", tech_missing) + "\n")


# -------------------------
# main
# -------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Audit EDTA outputs by scanning sample directories; incomplete outputs are TECH_FAIL.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("-d", "--dir", required=True, help="Root directory containing sample dirs")
    p.add_argument("-o", "--out", default="audit_out", help="Output directory")
    p.add_argument("-m", "--max-depth", type=int, default=6, help="Max recursion depth (-1 unlimited)")
    p.add_argument("--no-recursive", action="store_true", help="Disable recursive scan")
    p.add_argument("-q", "--quiet", action="store_true")
    p.add_argument("--profile", default="plant", choices=["plant","animal","fungi"])
    p.add_argument("--sample-regex", default=r"^GC[AF]_\d+\.\d+.*_genomic$",
                   help="Regex to detect sample directories (basename match)")

    # tech options
    p.add_argument("--no-require-anno-dir", action="store_true",
                   help="Do not require *.EDTA.anno/ for technical PASS")
    return p


def main():
    args = build_parser().parse_args()
    root = os.path.abspath(args.dir)
    outdir = os.path.abspath(args.out)
    recursive = not args.no_recursive
    quiet = args.quiet

    plant_thr = PlantThresholds()
    animal_thr = AnimalThresholds()
    fungi_thr = FungiThresholds()

    tech_cfg = TechConfig(require_anno_dir=(not args.no_require_anno_dir))

    eprint(f"[{now_ts()}] Scanning sample dirs: {root} (recursive={1 if recursive else 0}, max_depth={args.max_depth})", quiet)
    eprint(f"[{now_ts()}] Profile: {args.profile}", quiet)
    eprint(f"[{now_ts()}] sample_regex: {args.sample_regex}", quiet)

    sample_dirs = find_sample_dirs(root, recursive=recursive, max_depth=args.max_depth, sample_regex=args.sample_regex)
    eprint(f"[{now_ts()}] Found sample dirs: {len(sample_dirs)}", quiet)

    if not sample_dirs:
        eprint(f"[{now_ts()}] ERROR: No sample dirs matched regex under root.", quiet=False)
        sys.exit(2)

    records: List[Record] = []
    for sd in sample_dirs:
        records.append(classify_one(sd, args.profile, plant_thr, animal_thr, fungi_thr, tech_cfg))

    pass_rows = [r for r in records if r.overall_status == "PASS"]
    suspect_rows = [r for r in records if r.overall_status == "SUSPECT"]
    fail_rows = [r for r in records if r.overall_status == "FAIL"]

    os.makedirs(outdir, exist_ok=True)
    write_tsv(os.path.join(outdir, "all.tsv"), records)
    write_tsv(os.path.join(outdir, "pass.tsv"), pass_rows)
    write_tsv(os.path.join(outdir, "suspect.tsv"), suspect_rows)
    write_tsv(os.path.join(outdir, "fail.tsv"), fail_rows)
    write_summary(os.path.join(outdir, "summary.txt"), records, args.profile)

    eprint(f"[{now_ts()}] Done.", quiet)
    eprint(f"[{now_ts()}] Total samples: {len(records)} | PASS: {len(pass_rows)} | SUSPECT: {len(suspect_rows)} | FAIL: {len(fail_rows)}", quiet)


if __name__ == "__main__":
    main()

