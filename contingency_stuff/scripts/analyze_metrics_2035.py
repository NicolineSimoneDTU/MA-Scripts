# analyse_metrics_2035.py

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


### MAIN METRICS

def analyze_metrics_network(n, outage_mask):

    snapshots = n.snapshots
    weights = n.snapshot_weightings.generators 
    blocks = get_30d_blocks(outage_mask, snapshots)

    # Bus groups
    buses = n.buses

    ac_buses = buses.index[buses.carrier == "AC"]
    lv_buses = buses.index[buses.carrier == "low voltage"]
    elec_buses = ac_buses.union(lv_buses)

    heat_buses = buses.index[
        buses.carrier.isin(
            ["rural heat", "urban central heat", "urban decentral heat"]
        )
    ]

    prices = n.buses_t.marginal_price

    ### 
    # ELECTRICITY PRICES (per bus)
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
    # HEAT PRICES
    ### 
    heat_price_results = {}

    for bus in heat_buses:
        series = prices[bus]

        heat_price_results[bus] = {
            "outage_avg": float(series.loc[outage_mask].mean()),
            "outage_max": float(series.loc[outage_mask].max()),
            "blocks": {},
        }

        for i, (start, end) in enumerate(blocks, 1):
            mask = (snapshots >= start) & (snapshots <= end)
            heat_price_results[bus]["blocks"][f"block_{i}"] = {
                "avg": float(series.loc[mask].mean()),
                "max": float(series.loc[mask].max()),
            }

    ###
    # ELECTRICITY LOAD SHEDDING
    ###
    elec_ls = n.generators[n.generators.carrier == "elec load shed"]
    elec_ls_results = {"system": {}, "per_bus": {}}

    dispatch_elec = pd.DataFrame(0.0, index=snapshots, columns=[])

    if not elec_ls.empty:
        dispatch_elec = n.generators_t.p[elec_ls.index]

        ls_ts = dispatch_elec.sum(axis=1)
        elec_ls_results["system"]["outage_MWh"] = float((ls_ts.loc[outage_mask] * weights.loc[outage_mask]).sum())


        elec_ls_results["system"]["blocks"] = {}
        for i, (start, end) in enumerate(blocks, 1):
            mask = (snapshots >= start) & (snapshots <= end)
            elec_ls_results["system"]["blocks"][f"block_{i}"] = float(
                ls_ts.loc[mask].sum()
            )

        for gen, bus in zip(elec_ls.index, elec_ls.bus):
            elec_ls_results["per_bus"].setdefault(bus, {})
            elec_ls_results["per_bus"][bus]["outage_MWh"] = float(
                (dispatch_elec[gen].loc[outage_mask] * weights.loc[outage_mask]).sum()
            )

    ###
    # DK2 (Zealand) ONLY METRICS
    ###

    # DK1 buses correspond to Zealand / DK2 region
    dk2_buses = buses.index[buses.index.str.startswith("DK1")]
    dk2_ls_gens = elec_ls[elec_ls.bus.isin(dk2_buses)].index

    if len(dk2_ls_gens) > 0:
        ls_ts_dk2 = dispatch_elec[dk2_ls_gens].sum(axis=1)
        elec_ls_dk2_outage = float(
            (ls_ts_dk2.loc[outage_mask] * weights.loc[outage_mask]).sum()
        )
    else:
        elec_ls_dk2_outage = 0.0


    ###
    # HEAT LOAD SHEDDING
    ###
    heat_ls = n.generators[n.generators.carrier == "heat shed"]
    heat_ls_results = {"system": {}, "per_bus": {}}

    dispatch_heat = pd.DataFrame(0.0, index=snapshots, columns=[])

    if not heat_ls.empty:
        dispatch_heat = n.generators_t.p[heat_ls.index]

        heat_ts = dispatch_heat.sum(axis=1)
        heat_ls_results["system"]["outage_MWh"] = float(
            (heat_ts.loc[outage_mask] * weights.loc[outage_mask]).sum()
        )

        # blocks
        heat_ls_results["system"]["blocks"] = {}
        for i, (start, end) in enumerate(blocks, 1):
            mask = (snapshots >= start) & (snapshots <= end)
            heat_ls_results["system"]["blocks"][f"block_{i}"] = float(heat_ts.loc[mask].sum())

        # per bus
        for gen, bus in zip(heat_ls.index, heat_ls.bus):
            heat_ls_results["per_bus"].setdefault(bus, {})
            heat_ls_results["per_bus"][bus]["outage_MWh"] = float(
                (dispatch_heat[gen].loc[outage_mask] * weights.loc[outage_mask]).sum()
            )
    
    # DK2 only 
    dk2_heat_ls = heat_ls.index[heat_ls.bus.isin(dk2_buses)]

    if len(dk2_heat_ls) > 0:
        heat_ls_dk2_outage = float(
            (dispatch_heat.loc[outage_mask, dk2_heat_ls]
            .mul(weights.loc[outage_mask], axis=0)
            .sum()
            .sum())
        )
    else:
        heat_ls_dk2_outage = 0.0


    
    ###
    # BIOMASS LOAD SHEDDING
    ###
    bio_ls = n.generators[n.generators.carrier == "biomass shed"]
    bio_ls_results = {"system": {}, "per_bus": {}}

    dispatch_bio = pd.DataFrame(0.0, index=snapshots, columns=[])

    if not bio_ls.empty:
        dispatch_bio = n.generators_t.p[bio_ls.index]
        ls_ts_bio = dispatch_bio.sum(axis=1)
        bio_ls_results["system"]["outage_MWh"] = float(
            (ls_ts_bio.loc[outage_mask] * weights.loc[outage_mask]).sum()
        )

        bio_ls_results["system"]["blocks"] = {}
        for i, (start, end) in enumerate(blocks, 1):
            mask = (snapshots >= start) & (snapshots <= end)
            bio_ls_results["system"]["blocks"][f"block_{i}"] = float(
                ls_ts_bio.loc[mask].sum()
            )

        # per bus
        for gen, bus in zip(bio_ls.index, bio_ls.bus):
            bio_ls_results["per_bus"].setdefault(bus, {})
            bio_ls_results["per_bus"][bus]["outage_MWh"] = float(
                (dispatch_bio[gen].loc[outage_mask] * weights.loc[outage_mask]).sum()
            )

    
    # DK2 only biomass ls
    dk2_bio_ls = bio_ls[bio_ls.bus.isin(dk2_buses)].index

    if len(dk2_bio_ls) > 0:
        bio_ls_dk2_outage = float(
            (dispatch_bio.loc[outage_mask, dk2_bio_ls]
            .mul(weights.loc[outage_mask], axis=0)
            .sum()
            .sum())
        )
    else:
        bio_ls_dk2_outage = 0.0


    ###
    # LOAD SHEDDING DURATION METRICS
    ###
    shed_hours = ls_ts_dk2 > 1e-6
    shed_hours_outage = shed_hours & outage_mask

    dk2_duration = {
        "hours_with_shed_outage": int(shed_hours_outage.sum()),
        "longest_streak_outage": longest_run(shed_hours_outage),
    }

    ###
    # NET IMPORT / EXPORT (outage only)
    ###

    flow_results = {}

    for bus in elec_buses:
        net_flow = pd.Series(0.0, index=snapshots)

        # AC lines 
        for line, row in n.lines.iterrows():
            p0 = n.lines_t.p0[line]

            if row.bus0 == bus:
                # positive p0 = export
                net_flow += p0
            elif row.bus1 == bus:
                # positive p0 = import
                net_flow -= p0

        # DC links
        for link, row in n.links.iterrows():
            p0 = n.links_t.p0[link]

            if row.bus0 == bus:
                net_flow += p0
            elif row.bus1 == bus:
                net_flow -= p0

        # Outage aggregation 
        net_flow_outage = net_flow.loc[outage_mask]

        import_MWh = float(
            (-net_flow_outage.clip(upper=0))
            .mul(weights)
            .sum()
        )

        export_MWh = float(
            (net_flow_outage.clip(lower=0))
            .mul(weights)
            .sum()
        )

        net_import_MWh = import_MWh - export_MWh

        flow_results[bus] = {
            "import_MWh": import_MWh,
            "export_MWh": export_MWh,
            "net_import_MWh": net_import_MWh,
        }


    

    ###
    # RESILIENCE INDEX (RI-A, DK2)
    ###
    outage_hours = outage_mask.sum()
    eens_total_dk2 = elec_ls_dk2_outage + heat_ls_dk2_outage + bio_ls_dk2_outage

    ri_a_dk2 = eens_total_dk2 / outage_hours if outage_hours > 0 else 0.0

    return {
        ###
        # CORE RESILIENCE METRIC
        ###
        "ri_a_dk2": ri_a_dk2,

        ###
        # LOAD SHEDDING
        ###
        "electricity_load_shed": elec_ls_results,
        "heat_load_shed": heat_ls_results,
        "biomass_load_shed": bio_ls_results,

        # DK2 summary 
        "dk2_eens": {
            "electricity_MWh": elec_ls_dk2_outage,
            "heat_MWh": heat_ls_dk2_outage,
            "biomass_MWh": bio_ls_dk2_outage,
            "total_MWh": eens_total_dk2,
        },

        "dk2_duration": dk2_duration,

        ###
        # PRICES
        ###
        "electricity_prices": {
            "system_outage_avg": system_elec_price_outage,
            "per_bus": elec_price_results,
        },
        "heat_prices": {
            "per_bus": heat_price_results,
        },

        ###
        # FLOWS
        ###
        "net_flows": flow_results,
    }    


# ###
# SNAKEMAKE ENTRY POINT
# ###

if __name__ == "__main__":

    n = pypsa.Network(snakemake.input.net)

    start = pd.to_datetime(snakemake.params.start)
    end = start + pd.Timedelta(days=snakemake.params.duration)

    outage_mask = (n.snapshots >= start) & (n.snapshots <= end)

    results = analyze_metrics_network(n, outage_mask)

    with open(snakemake.output[0], "w") as f:
        json.dump(results, f, indent=4)

    print("Analysis metrics written successfully")