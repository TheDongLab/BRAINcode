from pathlib import Path

def summarize_gangstr_vs_trgt(base=None):
    if base is None:
        base = Path.home() / "donglab/data/target_ALS/WGS_LR/repeat_comparison_gangSTR_vs_TRGT"
    else:
        base = Path(base)

    infile = base / "comparison/all_calls.tsv"

    # Callability
    master = gang_called = trgt_called = both_called = 0
    gang_only = trgt_only = neither = 0

    # Agreement groups
    groups = {
        "exact": {"N": 0, "gang_q_sum": 0, "q_lt_09": 0, "q_lt_05": 0,
                  "trgt_purity_sum": 0, "purity_lt_09": 0, "weak_sd": 0},
        "within_1": {"N": 0, "gang_q_sum": 0, "q_lt_09": 0, "q_lt_05": 0,
                     "trgt_purity_sum": 0, "purity_lt_09": 0, "weak_sd": 0},
        "within_2": {"N": 0, "gang_q_sum": 0, "q_lt_09": 0, "q_lt_05": 0,
                     "trgt_purity_sum": 0, "purity_lt_09": 0, "weak_sd": 0},
        "gt_2": {"N": 0, "gang_q_sum": 0, "q_lt_09": 0, "q_lt_05": 0,
                 "trgt_purity_sum": 0, "purity_lt_09": 0, "weak_sd": 0},
    }

    # >2-repeat, high-confidence TRGT subset
    trgt_high_gt2 = 0
    trgt_high_gang_low = 0
    trgt_high_gang_mid = 0
    both_high = 0

    def called(gt):
        return gt not in ("", ".", "./.", ".|.")

    def pair(x):
        try:
            v = sorted(float(y) for y in x.split(","))
            return v if len(v) == 2 else None
        except (ValueError, AttributeError):
            return None

    with open(infile) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if len(x) < 19:
                continue

            master += 1

            # all_calls.tsv:
            # 9  GangSTR GT
            # 10 GangSTR REPCN
            # 13 GangSTR Q
            # 14 TRGT GT
            # 15 TRGT MC
            # 18 TRGT SD
            # 19 TRGT AP
            ggt = x[8]
            grepcn = x[9]
            gq = x[12]
            tgt = x[13]
            tmc = x[14]
            tsd = x[17]
            tap = x[18]

            gc = called(ggt)
            tc = called(tgt)

            if gc:
                gang_called += 1
            if tc:
                trgt_called += 1

            if gc and tc:
                both_called += 1
            elif gc:
                gang_only += 1
            elif tc:
                trgt_only += 1
            else:
                neither += 1

            if not (gc and tc):
                continue

            g = pair(grepcn)
            t = pair(tmc)
            p = pair(tap)
            sd = pair(tsd)

            if g is None or t is None or p is None or sd is None:
                continue

            d1 = abs(g[0] - t[0])
            d2 = abs(g[1] - t[1])
            d = max(d1, d2)

            if d == 0:
                group = "exact"
            elif d <= 1:
                group = "within_1"
            elif d <= 2:
                group = "within_2"
            else:
                group = "gt_2"

            q = float(gq)
            min_purity = min(p)
            min_sd = min(sd)

            z = groups[group]
            z["N"] += 1
            z["gang_q_sum"] += q
            z["trgt_purity_sum"] += min_purity

            if q < 0.9:
                z["q_lt_09"] += 1
            if q < 0.5:
                z["q_lt_05"] += 1
            if min_purity < 0.9:
                z["purity_lt_09"] += 1
            if min_sd < 3:
                z["weak_sd"] += 1

            # High-confidence TRGT >2-repeat disagreement:
            # purity >=0.9 on BOTH alleles
            # >=3 HiFi reads supporting BOTH alleles
            if d > 2 and min_purity >= 0.9 and min_sd >= 3:
                trgt_high_gt2 += 1

                if q < 0.5:
                    trgt_high_gang_low += 1
                elif q < 0.9:
                    trgt_high_gang_mid += 1
                else:
                    both_high += 1

    comparable = sum(z["N"] for z in groups.values())

    print("\n=== CATALOG / CALLABILITY ===")
    print(f"MASTER              {master:,}")
    print(f"GangSTR called      {gang_called:,} ({gang_called/master*100:.2f}%)")
    print(f"TRGT called         {trgt_called:,} ({trgt_called/master*100:.2f}%)")
    print(f"Both called         {both_called:,} ({both_called/master*100:.2f}%)")
    print(f"GangSTR only        {gang_only:,}")
    print(f"TRGT only           {trgt_only:,}")
    print(f"Neither called      {neither:,}")

    print("\n=== AGREEMENT AMONG BOTH-CALLED LOCI ===")
    labels = {
        "exact": "Exact",
        "within_1": "Differ by 1 repeat",
        "within_2": "Differ by 2 repeats",
        "gt_2": "Differ by >2 repeats",
    }

    for key in ["exact", "within_1", "within_2", "gt_2"]:
        n = groups[key]["N"]
        print(f"{labels[key]:22s} {n:>9,} ({n/comparable*100:5.2f}%)")

    print("\n=== CONFIDENCE BY AGREEMENT GROUP ===")
    header = (
        f"{'Group':22s} {'N':>9s} {'Mean G-Q':>10s} "
        f"{'G-Q<0.9':>16s} {'G-Q<0.5':>16s} "
        f"{'Mean T-purity':>14s} {'T-purity<0.9':>18s} "
        f"{'T weak reads':>16s}"
    )
    print(header)

    for key in ["exact", "within_1", "within_2", "gt_2"]:
        z = groups[key]
        n = z["N"]

        mean_q = z["gang_q_sum"] / n
        mean_p = z["trgt_purity_sum"] / n

        print(
            f"{labels[key]:22s} "
            f"{n:9,d} "
            f"{mean_q:10.3f} "
            f"{z['q_lt_09']:7,d} ({z['q_lt_09']/n*100:5.1f}%) "
            f"{z['q_lt_05']:7,d} ({z['q_lt_05']/n*100:5.1f}%) "
            f"{mean_p:14.3f} "
            f"{z['purity_lt_09']:7,d} ({z['purity_lt_09']/n*100:5.1f}%) "
            f"{z['weak_sd']:7,d} ({z['weak_sd']/n*100:5.1f}%)"
        )

    print("\n=== >2-REPEAT DISAGREEMENTS WITH HIGH-CONFIDENCE TRGT ===")
    print(
        "High-confidence TRGT = min purity >=0.9 "
        "AND >=3 HiFi reads supporting each allele"
    )
    print(f"Total                    {trgt_high_gt2:,}")

    if trgt_high_gt2:
        print(
            f"GangSTR Q <0.5           {trgt_high_gang_low:,} "
            f"({trgt_high_gang_low/trgt_high_gt2*100:.1f}%)"
        )
        print(
            f"GangSTR 0.5<=Q<0.9       {trgt_high_gang_mid:,} "
            f"({trgt_high_gang_mid/trgt_high_gt2*100:.1f}%)"
        )
        print(
            f"GangSTR Q >=0.9          {both_high:,} "
            f"({both_high/trgt_high_gt2*100:.1f}%)"
        )

    return {
        "master": master,
        "gangstr_called": gang_called,
        "trgt_called": trgt_called,
        "both_called": both_called,
        "gangstr_only": gang_only,
        "trgt_only": trgt_only,
        "neither_called": neither,
        "groups": groups,
        "high_confidence_trgt_gt2": trgt_high_gt2,
        "high_trgt_low_gangstr": trgt_high_gang_low,
        "high_trgt_mid_gangstr": trgt_high_gang_mid,
        "both_high_confidence": both_high,
    }
