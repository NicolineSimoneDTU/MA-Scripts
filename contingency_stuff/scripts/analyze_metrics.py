# analyse_metrics.py

import pypsa
import pandas as pd
import numpy as np
import json
from snakemake.script import snakemake


# Helpers
def get_30d_blocks(outage_mask, snapshots):
    blocks = []
    outage_snaps = snapshots[outage_mask]

    if outage_snaps.empty:
        return blocks

    start = outage_snaps[0]
    final = outage_snaps[-1]

    while start < final:
        end = min(start + pd.Timedelta(days=30), final)
        blocks.append((start, end))
        start = end

    return blocks


def longest_run(mask):
    if mask.empty:
        return 0
    groups = (mask == False).cumsum()
    return int(mask.groupby(groups).sum().max() or 0)


####
# MAIN ANALYSIS
###

def analyze_metrics_network(n, outage_mask):

    snapshots = n.snapshots
    blocks = get_30d_blocks(outage_mask, snapshots)

    # buses (electricity only)
    buses = n.buses
    elec_buses = buses.index
    prices = n.buses_t.marginal_price

    ###
    # ELECTRICITY PRICES
    ###
    elec_price_results = {}

    for bus in elec_buses:
        series = prices[bus]

        elec_price_results[bus] = {
            "outage_avg": float(series.loc[outage_mask].mean()),
            "outage_max": float(series.loc[outage_mask].max()),
            "blocks": {},
        }

        for i, (start, end) in enumerate(blocks, 1):
            mask = (snapshots >= start) & (snapshots <= end)
            elec_price_results[bus]["blocks"][f"block_{i}"] = {
                "avg": float(series.loc[mask].mean()),
                "max": float(series.loc[mask].max()),
            }

    system_elec_price_outage = (
        float(np.mean([v["outage_avg"] for v in elec_price_results.values()]))
        if elec_price_results else 0.0
    )

    ###
    # ELECTRICITY LOAD SHEDDING
    ###
    elec_ls = n.generators[n.generators.carrier == "elec load shed"]
    elec_ls_results = {"system": {}, "per_bus": {}}

    dispatch = pd.DataFrame(0.0, index=snapshots, columns=[])

    if not elec_ls.empty:
        dispatch = n.generators_t.p[elec_ls.index]
        ls_ts = dispatch.sum(axis=1)

        elec_ls_results["system"]["outage_MWh"] = float(
            ls_ts.loc[outage_mask].sum()
        )

        elec_ls_results["system"]["blocks"] = {}
        for i, (start, end) in enumerate(blocks, 1):
            mask = (snapshots >= start) & (snapshots <= end)
            elec_ls_results["system"]["blocks"][f"block_{i}"] = float(
                ls_ts.loc[mask].sum()
            )

        for gen, bus in zip(elec_ls.index, elec_ls.bus):
            elec_ls_results["per_bus"].setdefault(bus, {})
            elec_ls_results["per_bus"][bus]["outage_MWh"] = float(
                dispatch[gen].loc[outage_mask].sum()
            )

    ###
    # DK2 (Zealand) METRICS
    ###
    dk2_buses = buses.index[buses.index.str.startswith("DK1")] # DK1 buses correspond to Zealand / DK2 region
    dk2_ls_gens = elec_ls[elec_ls.bus.isin(dk2_buses)].index

    if len(dk2_ls_gens) > 0:
        ls_ts_dk2 = dispatch[dk2_ls_gens].sum(axis=1)
        eens_dk2 = float(ls_ts_dk2.loc[outage_mask].sum())
    else:
        ls_ts_dk2 = pd.Series(0.0, index=snapshots)
        eens_dk2 = 0.0

    shed_dk2 = ls_ts_dk2 > 1e-6
    shed_dk2_outage = shed_dk2 & outage_mask

    dk2_duration = {
        "hours_with_shed_outage": int(shed_dk2_outage.sum()),
        "longest_streak_outage": longest_run(shed_dk2_outage),
    }

    ###
    # NET IMPORT / EXPORT (outage period)
    ###
    flow_results = {}

    for bus in elec_buses:
        net_flow = pd.Series(0.0, index=snapshots)

        for line, row in n.lines.iterrows():
            p0 = n.lines_t.p0[line]
            if row.bus0 == bus:
                net_flow += p0
            elif row.bus1 == bus:
                net_flow -= p0

        for link, row in n.links.iterrows():
            p0 = n.links_t.p0[link]
            if row.bus0 == bus:
                net_flow += p0
            elif row.bus1 == bus:
                net_flow -= p0

        nf = net_flow.loc[outage_mask]
        import_MWh = float((-nf[nf < 0]).sum())
        export_MWh = float((nf[nf > 0]).sum())

        flow_results[bus] = {
            "import_MWh": import_MWh,
            "export_MWh": export_MWh,
            "net_import_MWh": import_MWh - export_MWh,
        }

    ###
    # RESILIENCE INDEX (RI-A, DK2)
    ###
    outage_hours = outage_mask.sum()
    ri_a_dk2 = eens_dk2 / outage_hours if outage_hours > 0 else 0.0

    ###
    # RETURN
    ###
    return {
        "ri_a_dk2": ri_a_dk2,
        "dk2_duration": dk2_duration,
        "electricity_prices": {
            "system_outage_avg": system_elec_price_outage,
            "per_bus": elec_price_results,
        },
        "electricity_load_shed": elec_ls_results,
        "net_flows": flow_results,
    }


###
# SNAKEMAKE ENTRY POINT
###

if __name__ == "__main__":

    n = pypsa.Network(snakemake.input.net)

    start = pd.to_datetime(snakemake.params.start)
    end = start + pd.Timedelta(days=snakemake.params.duration)
    outage_mask = (n.snapshots >= start) & (n.snapshots <= end)

    results = analyze_metrics_network(n, outage_mask)

    with open(snakemake.output[0], "w") as f:
        json.dump(results, f, indent=4)

    print("2025 analysis metrics written successfully")