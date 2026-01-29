# build_summary_2025.py

import json
import pandas as pd
import pathlib
from snakemake.script import snakemake

DK2_PREFIX = "DK1" # DK2 / Zealand corresponds to DK1 buses


###
# Helpers
###

def is_dk2(bus):
    return bus.startswith(DK2_PREFIX)


###
# MAIN
###

if __name__ == "__main__":

    rows = []
    rows_dk2 = []

    for metrics_file in snakemake.input.metrics:
        scenario = pathlib.Path(metrics_file).parts[-2]

        with open(metrics_file, "r") as f:
            d = json.load(f)

        ###
        # FULL SYSTEM SUMMARY ROW
        ###
        row = {
            "scenario": scenario,
            
            # RESILIENCE (DK2-based)
            
            "RI_DK2": d.get("ri_a_dk2", 0.0),
            "dk2_hours_with_shed_outage": d.get("dk2_duration", {}).get(
                "hours_with_shed_outage", 0
            ),
            "dk2_longest_shed_streak": d.get("dk2_duration", {}).get(
                "longest_streak_outage", 0
            ),

            
            # SYSTEM ELECTRICITY METRICS
            
            "elec_price_system_outage_avg": d["electricity_prices"]["system_outage_avg"],
            "elec_ls_system_outage_MWh": (
                d["electricity_load_shed"]["system"].get("outage_MWh", 0.0)
            ),
        }

        ###
        # ELECTRICITY PRICES — per bus
        ###
        elec_prices = d["electricity_prices"]["per_bus"]

        for bus, info in elec_prices.items():
            row[f"elec_price_outage_avg__{bus}"] = info["outage_avg"]
            row[f"elec_price_outage_max__{bus}"] = info["outage_max"]

            for block, vals in info.get("blocks", {}).items():
                row[f"elec_price_{block}_avg__{bus}"] = vals["avg"]
                row[f"elec_price_{block}_max__{bus}"] = vals["max"]

        ###
        # ELECTRICITY LOAD SHEDDING — per bus + blocks
        ###
        elec_ls = d["electricity_load_shed"]

        for bus, val in elec_ls.get("per_bus", {}).items():
            row[f"elec_ls_outage_MWh__{bus}"] = val.get("outage_MWh", 0.0)

        for block, val in elec_ls.get("system", {}).get("blocks", {}).items():
            row[f"elec_ls_{block}_system_MWh"] = val

        ###
        # NET FLOWS (disruption period)
        ###
        for bus, vals in d["net_flows"].items():
            row[f"import_MWh__{bus}"] = vals["import_MWh"]
            row[f"export_MWh__{bus}"] = vals["export_MWh"]
            row[f"net_import_MWh__{bus}"] = vals["net_import_MWh"]

        rows.append(row)

        ###
        # DK2-ONLY SUMMARY ROW
        ###
        dk2_row = {
            "scenario": scenario,
            "RI_DK2": d.get("ri_a_dk2", 0.0),
            "dk2_hours_with_shed_outage": d.get("dk2_duration", {}).get(
                "hours_with_shed_outage", 0
            ),
            "dk2_longest_shed_streak": d.get("dk2_duration", {}).get(
                "longest_streak_outage", 0
            ),
        }

        ###
        # DK2 ELECTRICITY PRICE
        ###
        dk2_prices = [
            info["outage_avg"]
            for bus, info in elec_prices.items()
            if is_dk2(bus)
        ]

        dk2_row["dk2_elec_price_outage_avg"] = (
            float(sum(dk2_prices) / len(dk2_prices)) if dk2_prices else 0.0
        )

        ###
        # DK2 LOAD SHEDDING
        ###
        dk2_row["dk2_elec_ls_outage_MWh"] = sum(
            val.get("outage_MWh", 0.0)
            for bus, val in elec_ls.get("per_bus", {}).items()
            if is_dk2(bus)
        )

        ###
        # DK2 NET FLOWS
        ###
        dk2_flows = [
            d["net_flows"][bus]
            for bus in d["net_flows"]
            if is_dk2(bus)
        ]

        dk2_row["dk2_import_MWh"] = sum(f["import_MWh"] for f in dk2_flows)
        dk2_row["dk2_export_MWh"] = sum(f["export_MWh"] for f in dk2_flows)
        dk2_row["dk2_net_import_MWh"] = sum(f["net_import_MWh"] for f in dk2_flows)

        rows_dk2.append(dk2_row)

    ###
    # SAVE TABLES
    ###
    df = pd.DataFrame(rows)
    df.to_csv(snakemake.output[0], index=False)
    print("Summary table written:", snakemake.output[0])

    dk2_path = pathlib.Path(snakemake.output[0]).with_name("summary_table_DK2.csv")
    df_dk2 = pd.DataFrame(rows_dk2)
    df_dk2.to_csv(dk2_path, index=False)
    print("DK2-only summary written:", dk2_path)