# build_summary_2035.py

import json
import pandas as pd
import pathlib
from snakemake.script import snakemake


BUS_PRIORITY = ["DK0 0", "DK1 0"] # Preferred bus order


def sort_buses(buses):
    """Ensure DK1 / DK0 appear first, rest alphabetical"""
    return sorted(
        buses,
        key=lambda b: (BUS_PRIORITY.index(b)
                       if b in BUS_PRIORITY else 99, b)
    )


if __name__ == "__main__":

    rows = []
    rows_dk2 = []

    for metrics_file in snakemake.input.metrics:
        scenario = pathlib.Path(metrics_file).parts[-2]

        with open(metrics_file, "r") as f:
            d = json.load(f)

        row = {
            "scenario": scenario,
        }

        ###
        # RESILIENCE INDEX
        ###
        row["RI_DK2"] = d.get("ri_a_dk2", 0.0)

        row["dk2_hours_with_shed_outage"] = d.get("dk2_duration", {}
        ).get("hours_with_shed_outage", 0)

        row["dk2_longest_shed_streak"] = d.get("dk2_duration", {}
        ).get("longest_streak_outage", 0)

        ###
        # SYSTEM-LEVEL METRICS
        ###
        row["elec_price_system_outage_avg"] = d["electricity_prices"]["system_outage_avg"]

        row["elec_ls_system_outage_MWh"] = (
            d["electricity_load_shed"]["system"].get("outage_MWh", 0.0)
        )

        row["heat_ls_system_outage_MWh"] = (
            d["heat_load_shed"]["system"].get("outage_MWh", 0.0)
        )

        row["bio_ls_system_outage_MWh"] = (d.get("biomass_load_shed", {}).get("system", {}).get("outage_MWh", 0.0)
        )


        ###
        # ELECTRICITY PRICES — per bus
        ###
        elec_prices = d["electricity_prices"]["per_bus"]

        for bus in sort_buses(elec_prices.keys()):
            info = elec_prices[bus]

            row[f"elec_price_outage_avg__{bus}"] = info["outage_avg"]
            row[f"elec_price_outage_max__{bus}"] = info["outage_max"]

            # 30-day blocks
            for block, vals in info.get("blocks", {}).items():
                row[f"elec_price_{block}_avg__{bus}"] = vals["avg"]
                row[f"elec_price_{block}_max__{bus}"] = vals["max"]

        ###
        # HEAT PRICES — Pper bus
        ###
        heat_prices = d["heat_prices"]["per_bus"]

        for bus in sort_buses(heat_prices.keys()):
            info = heat_prices[bus]

            row[f"heat_price_outage_avg__{bus}"] = info["outage_avg"]
            row[f"heat_price_outage_max__{bus}"] = info["outage_max"]

            for block, vals in info.get("blocks", {}).items():
                row[f"heat_price_{block}_avg__{bus}"] = vals["avg"]
                row[f"heat_price_{block}_max__{bus}"] = vals["max"]

        ###
        # ELECTRICITY LOAD SHEDDING — per bus + blocks
        ###
        elec_ls = d["electricity_load_shed"]

        for bus, val in elec_ls.get("per_bus", {}).items():
            row[f"elec_ls_outage_MWh__{bus}"] = val.get("outage_MWh", 0.0)

        for block, val in elec_ls.get("system", {}).get("blocks", {}).items():
            row[f"elec_ls_{block}_system_MWh"] = val

        ###
        # BIOMASS LOAD SHEDDING — per bus + blocks
        ###
        bio_ls = d.get("biomass_load_shed", {})

        for bus, val in bio_ls.get("per_bus", {}).items():
            row[f"bio_ls_outage_MWh__{bus}"] = val.get("outage_MWh", 0.0)
        
        ###
        # NET IMPORT / EXPORT (outage periods)
        ###
        flows = d["net_flows"]

        for bus in sort_buses(flows.keys()):
            row[f"import_MWh__{bus}"] = flows[bus]["import_MWh"]
            row[f"export_MWh__{bus}"] = flows[bus]["export_MWh"]
            row[f"net_import_MWh__{bus}"] = flows[bus]["net_import_MWh"]

        rows.append(row)

        ###
        # DK2-ONLY SUMMARY ROW
        ###
        dk2_row = {
            "scenario": scenario,
            "RI_DK2": d.get("ri_a_dk2", 0.0),
            "dk2_hours_with_shed_outage": d.get("dk2_duration", {}).get("hours_with_shed_outage", 0),
            "dk2_longest_shed_streak": d.get("dk2_duration", {}).get("longest_streak_outage", 0),
        }

        ###
        # DK2 ELECTRICITY PRICES
        ###
        elec_prices = d["electricity_prices"]["per_bus"]

        ac_price = elec_prices.get("DK1 0", {}).get("outage_avg", 0.0)
        lv_price = elec_prices.get("DK1 0 low voltage", {}).get("outage_avg", 0.0)

        dk2_row["dk2_elec_price_outage_avg_ac"] = ac_price
        dk2_row["dk2_elec_price_outage_avg_lv"] = lv_price
        dk2_row["dk2_elec_price_outage_avg_mean"] = (
            (ac_price + lv_price) / 2 if (ac_price and lv_price) else max(ac_price, lv_price)
        )

        ###
        # DK2 HEAT PRICES (mean)
        ###
        heat_prices = d["heat_prices"]["per_bus"]
        dk2_heat_vals = [
            v["outage_avg"]
            for bus, v in heat_prices.items()
            if bus.startswith("DK1")
        ]

        dk2_row["dk2_heat_price_outage_avg"] = (
            float(sum(dk2_heat_vals) / len(dk2_heat_vals)) if dk2_heat_vals else 0.0
        )
        
        ###
        # DK2 LOAD SHEDDING
        ###
        dk2_row["dk2_elec_ls_outage_MWh"] = (
            d["electricity_load_shed"]
            .get("per_bus", {})
            .get("DK1 0", {})
            .get("outage_MWh", 0.0)
        )

        dk2_row["dk2_heat_ls_outage_MWh"] = sum(
            v.get("outage_MWh", 0.0)
            for bus, v in d["heat_load_shed"].get("per_bus", {}).items()
            if bus.startswith("DK1")
        )

        dk2_row["dk2_bio_ls_outage_MWh"] = sum(
            v.get("outage_MWh", 0.0)
            for bus, v in d.get("biomass_load_shed", {}).get("per_bus", {}).items()
            if bus.startswith("DK1")
        )

        ###
        # DK2 NET FLOWS
        ###
        flows = d["net_flows"].get("DK1 0", {})
        dk2_row["dk2_import_MWh"] = flows.get("import_MWh", 0.0)
        dk2_row["dk2_export_MWh"] = flows.get("export_MWh", 0.0)
        dk2_row["dk2_net_import_MWh"] = flows.get("net_import_MWh", 0.0)

        rows_dk2.append(dk2_row)

    ###
    # SAVE
    ###
    df = pd.DataFrame(rows)
    df.to_csv(snakemake.output[0], index=False)
    print("Summary table written:", snakemake.output[0])

    # DK2-only summary
    dk2_path = pathlib.Path(snakemake.output[0]).with_name("summary_table_DK2.csv")
    df_dk2 = pd.DataFrame(rows_dk2)
    df_dk2.to_csv(dk2_path, index=False)

    print("DK2 summary written:", dk2_path)

    