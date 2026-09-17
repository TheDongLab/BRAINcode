#!/usr/bin/env bash
#SBATCH --job-name=qtl_meta
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH -p day
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
set -euo pipefail

# Usage:
#   bash QTL.pipeline.sh [username]
#   bash QTL.pipeline.sh [username] --dry-run
#
# DAG per QTL:
#   prep_<QTL>.sh <tissue>
#       ├── run_<QTL>.sh <tissue>               [standard]
#       └── run_<QTL>.sh <tissue> interaction   [interaction]
#
#   all STANDARD QTL runs complete
#       ├── run_<QTL>_MR.sh -> run_<QTL>_SuSiE_coloc.sh -> run_coloc_plot.sh <QTL>
#       ├── run_<QTL>_SMR_HEIDI.sh -> run_smr_heidi_plot.sh <QTL>
#       └── run_read_RPM_normalized_plots_for_QTL.sh <QTL>
#
# Child scripts retain their own #SBATCH resources and existing log files.
# This wrapper creates only .done/.jobid workflow-state files.

# ── Arguments / user ──────────────────────────────────────────────────────────
DRY_RUN=false
TARGET_USER=""
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        -h|--help)
            echo "Usage: bash $(basename "$0") [username] [--dry-run]"
            exit 0
            ;;
        --*) echo "ERROR: Unknown option: $arg" >&2; exit 2 ;;
        *)
            [[ -z "$TARGET_USER" ]] || { echo "ERROR: Only one username may be supplied." >&2; exit 2; }
            TARGET_USER="$arg"
            ;;
    esac
done
TARGET_USER="${TARGET_USER:-$USER}"
TARGET_HOME="$(getent passwd "$TARGET_USER" 2>/dev/null | awk -F: 'NR==1{print $6}')"
TARGET_HOME="${TARGET_HOME:-/home/${TARGET_USER}}"
LAB_ROOT="$TARGET_HOME/donglab"
DATA_ROOT="$LAB_ROOT/data/target_ALS"
QTL_DIR="$LAB_ROOT/pipelines/scripts/QTL"
MR_DIR="$LAB_ROOT/pipelines/scripts/MR"
STATUS_DIR="$DATA_ROOT/QTL/.workflow_status"

# ── Configuration ─────────────────────────────────────────────────────────────
QTL_TYPES=(eQTL sQTL cQTL)
TISSUES=(
    Motor_Cortex
    Frontal_Cortex
    Cervical_Spinal_Cord
    Lumbar_Spinal_Cord
    Cerebellum
)

# Interaction mode. Standard mode always runs for all five tissues.
INTERACTION_TISSUES_eQTL=(
    Motor_Cortex Frontal_Cortex Cervical_Spinal_Cord Lumbar_Spinal_Cord Cerebellum
)
INTERACTION_TISSUES_sQTL=(
    Motor_Cortex Frontal_Cortex Cervical_Spinal_Cord Lumbar_Spinal_Cord Cerebellum
)
INTERACTION_TISSUES_cQTL=(
    Motor_Cortex Frontal_Cortex Cervical_Spinal_Cord Lumbar_Spinal_Cord Cerebellum
)

# ── Helpers ───────────────────────────────────────────────────────────────────
status_file() { echo "$STATUS_DIR/$1.$2.done"; }
jobid_file()  { echo "$STATUS_DIR/$1.$2.jobid"; }
join_colon()  { local IFS=:; echo "$*"; }

interaction_tissues() {
    case "$1" in
        eQTL) printf '%s\n' "${INTERACTION_TISSUES_eQTL[@]}" ;;
        sQTL) printf '%s\n' "${INTERACTION_TISSUES_sQTL[@]}" ;;
        cQTL) printf '%s\n' "${INTERACTION_TISSUES_cQTL[@]}" ;;
    esac
}

contains() {
    local needle="$1" x; shift
    for x in "$@"; do [[ "$x" == "$needle" ]] && return 0; done
    return 1
}

active_marker() {
    local jf="$1" jid reason
    [[ -s "$jf" ]] || return 1
    jid="$(cat "$jf")"
    [[ -n "$jid" ]] || return 1
    if squeue -h -j "$jid" -o '%T' 2>/dev/null | grep -q .; then
        reason="$(squeue -h -j "$jid" -o '%r' 2>/dev/null | head -1 || true)"
        if [[ "$reason" == *DependencyNeverSatisfied* ]]; then
            scancel "$jid" >/dev/null 2>&1 || true
            rm -f "$jf"
            return 1
        fi
        echo "$jid"
        return 0
    fi
    return 1
}

submit_marker() {
    local mode="$1" sf="$2" jf="$3" deps="$4" depkind marker
    [[ "$mode" == soft ]] && depkind=afterany || depkind=afterok
    marker="$(sbatch --parsable \
        --dependency="${depkind}:${deps}" \
        --output=/dev/null --error=/dev/null \
        --wrap="mkdir -p '$STATUS_DIR'; touch '$sf'" | cut -d';' -f1)"
    echo "$marker" > "$jf"
    echo "$marker"
}

# submit_stage QTL STAGE hard|soft DEPENDENCY SCRIPT [ARGS...]
# hard: .done only after successful child job
# soft: plotting branch is terminal/best-effort; attempt is considered complete
submit_stage() {
    local qtl="$1" stage="$2" mode="$3" dep="$4" script="$5"; shift 5
    local sf jf old child marker
    sf="$(status_file "$qtl" "$stage")"
    jf="$(jobid_file "$qtl" "$stage")"

    if [[ -f "$sf" ]]; then
        echo "[DONE]   $qtl $stage" >&2
        echo ""
        return
    fi
    if ! $DRY_RUN && old="$(active_marker "$jf" 2>/dev/null)"; then
        echo "[ACTIVE] $qtl $stage -> $old" >&2
        echo "$old"
        return
    fi
    if $DRY_RUN; then
        printf '[PLAN]   %s %s' "$qtl" "$stage" >&2
        [[ -n "$dep" ]] && printf ' afterok:%s' "$dep" >&2
        printf '\n         sbatch --chdir=%q %q' "$(dirname "$script")" "$script" >&2
        [[ $# -gt 0 ]] && printf ' %q' "$@" >&2
        printf '\n' >&2
        echo "DRYRUN_${qtl}_${stage}"
        return
    fi

    if [[ -n "$dep" ]]; then
        child="$(sbatch --parsable --chdir="$(dirname "$script")" \
            --dependency="afterok:${dep}" "$script" "$@" | cut -d';' -f1)"
    else
        child="$(sbatch --parsable --chdir="$(dirname "$script")" \
            "$script" "$@" | cut -d';' -f1)"
    fi
    marker="$(submit_marker "$mode" "$sf" "$jf" "$child")"
    echo "[SUBMIT] $qtl $stage -> job $child; status $marker" >&2
    echo "$marker"
}

submit_barrier() {
    local qtl="$1" stage="$2"; shift 2
    local sf jf old deps marker
    sf="$(status_file "$qtl" "$stage")"
    jf="$(jobid_file "$qtl" "$stage")"
    if [[ -f "$sf" ]]; then echo "[DONE]   $qtl $stage" >&2; echo ""; return; fi
    if ! $DRY_RUN && old="$(active_marker "$jf" 2>/dev/null)"; then
        echo "[ACTIVE] $qtl $stage -> $old" >&2
        echo "$old"
        return
    fi
    [[ $# -gt 0 ]] || {
        if $DRY_RUN; then echo "[PLAN]   restore $qtl $stage aggregate status" >&2
        else touch "$sf"
        fi
        echo ""
        return
    }
    deps="$(join_colon "$@")"
    if $DRY_RUN; then
        echo "[PLAN]   $qtl $stage barrier afterok:$deps" >&2
        echo "DRYRUN_${qtl}_${stage}"
        return
    fi
    marker="$(submit_marker hard "$sf" "$jf" "$deps")"
    echo "[SUBMIT] $qtl $stage barrier -> $marker" >&2
    echo "$marker"
}

# ── Preflight ─────────────────────────────────────────────────────────────────
preflight() {
    local failed=0 qtl t s cmd
    local -a scripts=()
    echo "============================================================"
    echo "QTL PIPELINE PREFLIGHT"
    echo "============================================================"
    echo "User      : $TARGET_USER"
    echo "Home      : $TARGET_HOME"
    echo "Data      : $DATA_ROOT"
    echo "QTL dir   : $QTL_DIR"
    echo "MR dir    : $MR_DIR"
    echo "Dry run   : $DRY_RUN"

    for cmd in sbatch squeue scancel getent awk grep cut; do
        command -v "$cmd" >/dev/null 2>&1 || { echo "[MISSING CMD] $cmd"; failed=1; }
    done
    for d in "$DATA_ROOT" "$DATA_ROOT/QTL" "$QTL_DIR" "$MR_DIR"; do
        [[ -d "$d" ]] || { echo "[MISSING DIR] $d"; failed=1; }
    done
    for t in "${TISSUES[@]}"; do
        [[ -d "$DATA_ROOT/$t" ]] && echo "[OK TISSUE] $t" || { echo "[MISSING TISSUE] $t"; failed=1; }
    done
    for qtl in "${QTL_TYPES[@]}"; do
        scripts+=(
            "$QTL_DIR/prep_${qtl}.sh"
            "$QTL_DIR/run_${qtl}.sh"
            "$MR_DIR/run_${qtl}_MR.sh"
            "$MR_DIR/run_${qtl}_SMR_HEIDI.sh"
            "$MR_DIR/run_${qtl}_SuSiE_coloc.sh"
        )
    done
    scripts+=(
        "$QTL_DIR/run_smr_heidi_plot.sh"
        "$QTL_DIR/run_coloc_plot.sh"
        "$QTL_DIR/run_read_RPM_normalized_plots_for_QTL.sh"
    )
    for s in "${scripts[@]}"; do
        [[ -f "$s" && -r "$s" ]] && echo "[OK SCRIPT] $s" || { echo "[MISSING SCRIPT] $s"; failed=1; }
    done
    for qtl in "${QTL_TYPES[@]}"; do
        s="$QTL_DIR/run_${qtl}.sh"
        [[ ! -f "$s" ]] || grep -qi interaction "$s" || echo "[WARNING] No literal 'interaction' found in $s"
    done
    [[ "$failed" -eq 0 ]] || { echo "PREFLIGHT FAILED — no jobs submitted."; exit 1; }
    echo "PREFLIGHT PASSED"
    echo "============================================================"
}

# ── QTL branch ────────────────────────────────────────────────────────────────
submit_qtl() {
    local qtl="$1" t prep run p s i standard_ready interaction_ready
    local standard_sf interaction_sf
    local -a standard_jobs=() interaction_jobs=() ints=()
    prep="$QTL_DIR/prep_${qtl}.sh"
    run="$QTL_DIR/run_${qtl}.sh"
    mapfile -t ints < <(interaction_tissues "$qtl")
    standard_sf="$(status_file "$qtl" QTL_standard)"
    interaction_sf="$(status_file "$qtl" QTL_interaction)"

    if [[ -f "$standard_sf" && -f "$interaction_sf" ]]; then
        echo "[DONE]   $qtl QTL standard + interaction" >&2
        echo ""
        return
    fi

    for t in "${TISSUES[@]}"; do
        p="$(submit_stage "$qtl" "prep_${t}" hard "" "$prep" "$t")"
        if [[ ! -f "$standard_sf" ]]; then
            s="$(submit_stage "$qtl" "run_${t}_standard" hard "$p" "$run" "$t")"
            [[ -n "$s" ]] && standard_jobs+=("$s")
        fi
        if [[ ! -f "$interaction_sf" ]] && contains "$t" "${ints[@]}"; then
            i="$(submit_stage "$qtl" "run_${t}_interaction" hard "$p" "$run" "$t" interaction)"
            [[ -n "$i" ]] && interaction_jobs+=("$i")
        fi
    done

    if [[ -f "$standard_sf" ]]; then standard_ready=""
    else standard_ready="$(submit_barrier "$qtl" QTL_standard "${standard_jobs[@]}")"
    fi
    if [[ -f "$interaction_sf" ]]; then interaction_ready=""
    else interaction_ready="$(submit_barrier "$qtl" QTL_interaction "${interaction_jobs[@]}")"
    fi
    : "$interaction_ready"
    echo "$standard_ready"
}

# ── Main ──────────────────────────────────────────────────────────────────────
preflight
$DRY_RUN || mkdir -p "$STATUS_DIR"

echo ""
echo "============================================================"
echo "SUBMITTING QTL WORKFLOW"
echo "============================================================"
for qtl in "${QTL_TYPES[@]}"; do
    echo ""
    echo "-------------------- $qtl --------------------"
    qtl_ready="$(submit_qtl "$qtl")"

    # These three branches need only the completed STANDARD raw QTL results.
    mr_ready="$(submit_stage "$qtl" MR hard "$qtl_ready" "$MR_DIR/run_${qtl}_MR.sh")"
    smr_ready="$(submit_stage "$qtl" SMR_HEIDI hard "$qtl_ready" "$MR_DIR/run_${qtl}_SMR_HEIDI.sh")"
    submit_stage "$qtl" RPM_plot soft "$qtl_ready" \
        "$QTL_DIR/run_read_RPM_normalized_plots_for_QTL.sh" "$qtl" >/dev/null

    # Required ordering: MR -> SuSiE -> coloc plot.
    susie_ready="$(submit_stage "$qtl" SuSiE_coloc hard "$mr_ready" "$MR_DIR/run_${qtl}_SuSiE_coloc.sh")"

    # Terminal plot branches. Missing/no-significant-result plots never gate the pipeline.
    submit_stage "$qtl" SMR_HEIDI_plot soft "$smr_ready" "$QTL_DIR/run_smr_heidi_plot.sh" "$qtl" >/dev/null
    submit_stage "$qtl" coloc_plot soft "$susie_ready" "$QTL_DIR/run_coloc_plot.sh" "$qtl" >/dev/null
done

echo ""
echo "============================================================"
if $DRY_RUN; then
    echo "DRY RUN COMPLETE — no jobs submitted and no status files written."
else
    echo "WORKFLOW SUBMITTED"
    echo "Status: $STATUS_DIR"
    echo "Monitor: squeue -u $USER"
fi
echo "============================================================"
