import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# load network function
def load_network(file_path):

    """
    Load a PyPSA network from a given file path. 
    """

    n_new = pypsa.Network(file_path)

    return n_new

# solving dispatch optimisation (with fixed capacities and infeasible flag)
def solve_dispatch(n):
    n.optimize.fix_optimal_capacities()
    n.optimize(
        snapshots=n.snapshots,
        solver_name="gurobi",
        solver_options={
            "OutputFlag": 0,
        },
    )

    # flag if not optimal
    status = getattr(n.model, "status", None)
    term = getattr(n.model, "termination_condition", None)
    print("Solver status:", status, "termination:", term)

    if str(term) not in ["optimal", "optimal_inaccurate"]:
        raise RuntimeError(f"Dispatch solve failed: {term}")

    return n

# interconnector disruption function
def apply_line_link_outages(n, line_list=None, link_list=None, outage_mask=None):
    """
    Apply line and link outages during a given time window.
    Restores normal availability outside the outage period.
    """

    if outage_mask is None:
        raise ValueError("You must provide an outage mask.")

    # 
    # LINES
    # 
    if line_list is not None:
        for line in line_list:

            # Create column if missing
            if line not in n.lines_t.s_max_pu.columns:
                # fill entire year with 0.7 (fully available)
                n.lines_t.s_max_pu[line] = 0.7 

            # Cut line during outage only
            n.lines_t.s_max_pu.loc[outage_mask, line] = 0.0

    # 
    # LINKS
    # 
    if link_list is not None:
        for link in link_list:

            # Create columns if missing — fill whole year with 1/-1
            if link not in n.links_t.p_max_pu.columns:
                n.links_t.p_max_pu[link] = 1.0
            if link not in n.links_t.p_min_pu.columns:
                n.links_t.p_min_pu[link] = -1.0

            # Apply outage only in time window
            n.links_t.p_max_pu.loc[outage_mask, link] = 0.0
            n.links_t.p_min_pu.loc[outage_mask, link] = 0.0

    return n


# cyber attack on offshore windfarms function
def apply_wind_attack_by_bus(n, buses, scaling_factor, outage_mask):
    """
    Apply cyber attack on wind farms located at specific buses.
    Restores normal availability outside the outage period.

    """

    # Find wind generators connected to the selected buses
    wind_gens = n.generators[
        (n.generators.carrier.str.contains("offwind", case=False)) &
        (n.generators.bus.isin(buses))
    ].index

    if len(wind_gens) == 0:
        print(f"No wind generators found at the selected buses {buses}.")
        return n

    # Reduce offshore availability using scaling factor
    for gen in wind_gens:
        n.generators_t.p_max_pu.loc[outage_mask, gen] *= scaling_factor

    print(f"Applied wind attack at buses {buses}. Affected generators: {list(wind_gens)}")
    print(f"Affected wind capacity (MW):",
      n.generators.loc[wind_gens, "p_nom"].sum())

    return n

# biomass supply shortage function 2025
def apply_biomass_shortage(n, buses, scaling_factor, outage_mask):

    """
    Apply a biomass shortage by scaling down the biomass generation capacity.
    Restores normal availability outside the outage period.
    
    """
    bio_mask = n.generators.carrier.str.contains("bio", case=False, na=False)
    bio_mask &= n.generators.bus.isin(buses)

    biomass_gens = n.generators.index[bio_mask]

    if biomass_gens.empty:
        print("No biomass generators found - skipping biomass shortage.")
        return n

    print("Biomass gens affected:", list(biomass_gens))

    for gen in biomass_gens:
        if gen not in n.generators_t.p_max_pu.columns:
            n.generators_t.p_max_pu[gen] = 1.0
        if gen not in n.generators_t.p_min_pu.columns:
            n.generators_t.p_min_pu[gen] = 0.0

        n.generators_t.p_max_pu.loc[outage_mask, gen] *= scaling_factor


    return n

# biomass supply shortage function 2035
def apply_biomass_shortage_2035(n, areas, scaling_factor, outage_mask):

    """
    Apply a biomass shortage by scaling down the biomass link capacity.
    Restores normal availability outside the outage period.
    
    """
    
    # Convert area to list
    if areas is None:
        area_list = None
    elif isinstance(areas, str):
        area_list = [areas]
    else:
        area_list = list(areas)

    # Find all biomass links 
    bio_links_df = n.links[n.links.carrier.str.contains("biomass", case=False, na=False)]
    all_bio_links = bio_links_df.index.tolist()

    print("Total biomass links:", len(all_bio_links))

    # Ensure all bio links have p_max_pu columns 
    for link in all_bio_links:
        if link not in n.links_t.p_max_pu.columns:
            n.links_t.p_max_pu[link] = 1.0
        if link not in n.links_t.p_min_pu.columns:
            n.links_t.p_min_pu[link] = 0.0

    # Filter by area 
    if area_list is not None:
        mask_area = bio_links_df.apply(
            lambda row: any(row.bus0.startswith(a) or row.bus1.startswith(a)
                            for a in area_list),
            axis=1
        )
        bio_links = bio_links_df[mask_area].index.tolist()
    else:
        bio_links = all_bio_links

    print("Biomass links to scale:", len(bio_links))

    # Apply scaling/shortage
    for link in bio_links:
        before = n.links_t.p_max_pu.loc[outage_mask, link].unique()
        n.links_t.p_max_pu.loc[outage_mask, link] *= scaling_factor
        after = n.links_t.p_max_pu.loc[outage_mask, link].unique()
        print(f"{link}: {before} → {after}")

    return n
