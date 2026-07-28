import argparse
import numpy as np
import pandas as pd

def main(csv_path: str, out_csv: str):
    df = pd.read_csv(csv_path)

    # Mean mapping (order) for instance
    m = (df.groupby(["instance_name", "order"])[["depth_out", "num_swaps"]]
            .mean()
            .reset_index())

    # Best for instance (minimum)
    m["depth_best"] = m.groupby("instance_name")["depth_out"].transform("min")
    m["swaps_best"] = m.groupby("instance_name")["num_swaps"].transform("min")

    # Variation respect the best (division by 0 return nan)
    m["depth_var_pct"] = np.where(
        m["depth_best"] != 0,
        (m["depth_out"] - m["depth_best"]) / m["depth_best"] * 100.0,
        np.nan
    )
    m["num_swaps_var_pct"] = np.where(
        m["swaps_best"] != 0,
        (m["num_swaps"] - m["swaps_best"]) / m["swaps_best"] * 100.0,
        np.nan
    )

    # Ranks (1 = best) for immediately confront 
    m["rank_depth"] = m.groupby("instance_name")["depth_out"].rank(method="min")
    m["rank_swaps"] = m.groupby("instance_name")["num_swaps"].rank(method="min")

    # order
    out = (m[[
            "instance_name", "order",
            "depth_out", "depth_var_pct",
            "num_swaps", "num_swaps_var_pct",
            "rank_depth", "rank_swaps"
        ]]
        .rename(columns={
            "depth_out": "depth_mean",
            "num_swaps": "num_swaps_mean"
        })
        .sort_values(
            by=["instance_name", "depth_var_pct", "num_swaps_var_pct", "order"],
            ascending=[True, True, True, True]
        ))

    out["depth_mean"] = out["depth_mean"].round(3)
    
    # depth_mean - min(depth_mean) / min(depth_mean) * 100
    out["depth_var_pct"] = out["depth_var_pct"].round(2)

    out["num_swaps_mean"] = out["num_swaps_mean"].round(3)
    
    # num_swaps_mean - min(num_swaps_mean) / min(num_swaps) * 100
    out["num_swaps_var_pct"] = out["num_swaps_var_pct"].round(2)

    out.to_csv(out_csv, index=False)
    
    print(f"Save: {out_csv}")
    #Show first line 
    print(out.head(12).to_string(index=False))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="csv path (es: esp_angelo_1.csv)")
    ap.add_argument("--out", default="mapping_variations_by_order.csv", help="CSV output")
    args = ap.parse_args()
    main(args.csv, args.out)
