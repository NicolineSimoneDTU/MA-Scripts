# scripts/add_load_shedding.py

import pypsa
import pandas as pd

if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.net)

    
    # Parameters from Snakemake / config
    elec_ls_cost = float(snakemake.params.elec_ls_cost)      # €/MWh
    heat_ls_factor = float(snakemake.params.heat_ls_factor)  
    add_heat_ls = bool(snakemake.params.add_heat_ls)         
    use_lv_buses = bool(snakemake.params.use_lv_buses)       # load shed should be added on low voltage buses for 2035
    add_bio_ls = bool(snakemake.params.add_bio_ls)           
    bio_ls_cost = elec_ls_cost * 0.25  

    
    # ELECTRICITY LOAD SHEDDING
    if use_lv_buses and "low voltage" in n.buses.carrier.values:
        # sector-coupled - LS only on LV buses
        elec_buses = n.buses.query('carrier == "low voltage"').index 
    else:
        # only-electricity - LS on electricity buses
        elec_buses = n.buses.query('carrier == "AC" or carrier == "DC"').index

    if len(elec_buses) == 0:
        print("No electricity buses found for load shedding.")
    else:
        if "elec load shed" not in n.carriers.index:
            n.add("Carrier", "elec load shed",
                  nice_name="Electricity load shedding",
                  color="#ff0000")

        n.add(
            "Generator",
            elec_buses,
            " load shed",
            bus=elec_buses,
            carrier="elec load shed",
            marginal_cost=elec_ls_cost,  # €/MWh
            p_nom=1e9,                   # MW 
        )

    # HEAT SHEDDING (only for sector-coupled)
    if add_heat_ls and any(n.buses.carrier.str.contains("heat")):
        heat_ls_cost = heat_ls_factor * elec_ls_cost

        if "heat shed" not in n.carriers.index:
            n.add("Carrier", "heat shed",
                  nice_name="Heat shedding",
                  color="#4F000B")

        for heat_carrier in ["rural heat", "urban central heat", "urban decentral heat"]:
            buses_i = n.buses.query('carrier == @heat_carrier').index
            if len(buses_i) == 0:
                continue

            n.add(
                "Generator",
                buses_i,
                " shedding",
                bus=buses_i,
                carrier="heat shed",
                marginal_cost=heat_ls_cost,
                p_nom=1e6,  # MW
            )
    
    # BIOMASS LOAD SHEDDING (only for sector-coupled)
    if add_bio_ls and any(n.buses.carrier.str.contains("bio")):
        
        if "biomass shed" not in n.carriers.index:
            n.add("Carrier", "biomass shed",
                  nice_name="Biomass shedding",
                  color="#8B4513")
    
        for bio_carrier in ["biogas", "solid biomass", "solid biomass for industry", "unsustainable bioliquids"]:
            bio_buses_i = n.buses.query('carrier == @bio_carrier').index
            if len(bio_buses_i) == 0:
                continue

            n.add(
                "Generator",
                bio_buses_i,
                " shedding",
                bus=bio_buses_i,
                carrier="biomass shed",
                marginal_cost=bio_ls_cost,
                p_nom=1e6,  # MW
            )

    # Save new network
    n.export_to_netcdf(snakemake.output[0])