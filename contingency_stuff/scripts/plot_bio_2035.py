# plot_bio_2035.py

import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os


# COLOR HANDLING

tech_colors = {
     'onwind': '#235ebc',
    'offwind-ac': '#6895dd',
    'offwind-dc': '#74c6f2',
    'offwind-float': '#b5e2fa',
    'solar': '#f9d002',
    'solar-hsat': '#fdb915',
    'OCGT': '#e0986c',
    'CCGT': '#a85522',
    'coal': '#545454',
    'lignite': '#826837',
    'nuclear': '#ff8c00',
    'biomass': '#baa741',
    'oil': '#c9c9c9',
    'hydro': '#298c81',
    'ror': '#3dbfb0',
    'unsustainable solid biomass': '#998622',
    'unsustainable bioliquids': '#32CD32',
    'unsustainable biogas': '#e3d37d',
    'solar rooftop': '#ffea80', 
    'solid biomass': '#baa741',
    'biogas': '#e3d37d',
    'gas': '#e05b09',
    'oil primary': '#d2d2d2',
    'BEV charger': '#baf238',
    'V2G': '#e5ffa8',
    'water pits charger': "#b36a5e",
    'water pits discharger': "#b37468",
    'gas pipeline': '#ebbca0',
    'urban central gas CHP': '#8d5e56',
    'gas for industry': '#853403',
    'DC': "#8a1caf",
    'SMR CC': '#4f1745',
    'solid biomass for industry': '#7a6d26',
    'urban central air heat pump': '#6cfb6b',
    'co2 sequestered': '#f2682f',
    'urban central aquifer thermal energy storage': "#6d00fc", 
    'non-sequestered HVC': '#8f79b5',
    'co2': '#f29dae',
    'methanol': '#FF7B00',
    'EV battery': '#baf238',
    'PHS': '#51dbcc',
    'urban central solar thermal': '#d7b24c',
    'rural solar thermal': '#f1c069', 
    'urban decentral solar thermal': '#e5bc5a',
    'urban central heat vent': '#a74747',
    'urban decentral heat vent': '#a33c3c', 
    'rural heat vent': '#ff5c5c', 
    
    'urban central water tanks': '#e9977d',
    'urban central water tanks charger': '#b57a67',
    'urban central water tanks discharger': '#b9816e',

    'urban central water pits': "#d96f4c",
    'urban central water pits charger': "#a85d47",
    'urban central water pits discharger': "#b36452",

    'rural water tanks': '#f7b7a3',

    'urban central water tanks': '#b9816e', 
    'urban decentral water tanks': '#baac9e',
    'home battery': '#80c944',
    'battery': '#ace37f',
    'elec load shed': '#ff0000',
    'heat shed': '#ff0009',
    'biomass shed': "#e600ff"
    
}

def get_color(carrier):
    return tech_colors.get(carrier, "#777777")


############################################################
# PLOTTING FUNCTIONS


def get_biomass_links(n, area=None):
    # All biomass links
    bio_links = n.links[
        n.links.carrier.str.contains("bio", case=False, na=False)
    ].index.tolist()

    # No area filtering
    if area is None:
        return bio_links

    # Turn single string into list
    if isinstance(area, str):
        area_list = [area]
    else:
        area_list = list(area)

    # Filter by area list
    filtered = []
    for l in bio_links:
        bus0 = n.links.at[l, "bus0"]
        bus1 = n.links.at[l, "bus1"]
        if any(bus0.startswith(a) or bus1.startswith(a) for a in area_list):
            filtered.append(l)

    return filtered


def plot_biomass_link_capacity(n, area=None, outage_start=None, outage_end=None):

    bio_links = get_biomass_links(n, area=area)

    if not bio_links:
        print("No biomass links found for this selection.")
        return None

    for link in bio_links:
        if link not in n.links_t.p_max_pu.columns:
            n.links_t.p_max_pu[link] = 1.0

    if outage_start is not None and outage_end is not None:
        outage_start = pd.to_datetime(outage_start)
        outage_end   = pd.to_datetime(outage_end)

    if isinstance(area, list):
        area_label = ", ".join(area)
    else:
        area_label = area if area else ""

    title_suffix = f" ({area_label})" if area_label else ""

    # Plot 
    fig, ax = plt.subplots(figsize=(14, 4))

    for lk in bio_links:
        ax.plot(
            n.links_t.p_max_pu.index,
            n.links_t.p_max_pu[lk],
            label=lk
        )

    ax.set_title(f"Biomass link availability (p_max_pu){title_suffix}")
    ax.set_ylabel("p_max_pu")
    ax.set_xlabel("Time")

    if outage_start is not None and outage_end is not None:
        ax.axvspan(outage_start, outage_end, color="red", alpha=0.15)

    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.grid(True, alpha=0.3)

    ax.legend(title="Biomass link", bbox_to_anchor=(1.01, 1), loc="upper left")

    plt.tight_layout()



#### dispatch plots
def plot_dispatch_by_carrier(n, outage_start, outage_end):
    dispatch_g = n.generators_t.p.copy()
    carriers_g = n.generators.carrier

    dispatch_s = n.storage_units_t.p.copy().clip(lower=0)
    carriers_s = n.storage_units.carrier

    dispatch = pd.concat([dispatch_g, dispatch_s], axis=1)
    carriers = pd.concat([carriers_g, carriers_s])

    dispatch_by_carrier = dispatch.groupby(carriers, axis=1).sum()

    outage_start = pd.to_datetime(outage_start)
    outage_end = pd.to_datetime(outage_end)
    dispatch_outage = dispatch_by_carrier.loc[outage_start:outage_end]

    fig, axes = plt.subplots(1, 2, figsize=(20, 6), sharey=True)

    colors = [get_color(c) for c in dispatch_by_carrier.columns]

    dispatch_by_carrier.plot.area(ax=axes[0], color=colors, linewidth=0, legend=False)
    axes[0].set_title("Full-Year Dispatch by Carrier")

    dispatch_outage.plot.area(ax=axes[1], color=colors, linewidth=0, legend=False)
    axes[1].set_title("Outage-Period Dispatch by Carrier")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5)

    plt.tight_layout(rect=[0, 0.08, 1, 1])

def plot_bus_dispatch_by_carrier(n, bus, outage_start, outage_end):

    gens_at_bus = n.generators[n.generators.bus == bus]
    stores_at_bus = n.storage_units[n.storage_units.bus == bus]

    gen_names = gens_at_bus.index.tolist()
    store_names = stores_at_bus.index.tolist()

    if len(gen_names) + len(store_names) == 0:
        print(f"No generators or storage found at bus '{bus}'.")
        return

    # extract time series
    dispatch_g = n.generators_t.p[gen_names] if gen_names else pd.DataFrame(index=n.snapshots)
    dispatch_s = n.storage_units_t.p[store_names] if store_names else pd.DataFrame(index=n.snapshots)

    if not dispatch_s.empty:
        dispatch_s = dispatch_s.clip(lower=0)

    dispatch = pd.concat([dispatch_g, dispatch_s], axis=1)

    carriers = pd.concat([gens_at_bus.carrier, stores_at_bus.carrier])
    dispatch_by_carrier = dispatch.groupby(carriers, axis=1).sum()

    outage_start = pd.to_datetime(outage_start)
    outage_end = pd.to_datetime(outage_end)
    dispatch_outage = dispatch_by_carrier.loc[outage_start:outage_end]

    colors = [get_color(c) for c in dispatch_by_carrier.columns]

    fig, axes = plt.subplots(1, 2, figsize=(22, 6), sharey=True)

    # Full period
    dispatch_by_carrier.plot.area(ax=axes[0], color=colors, linewidth=0)
    axes[0].set_title(f"Dispatch at {bus} - Full Year")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("MW")
    axes[0].grid(alpha=0.3)

    # Outage period
    dispatch_outage.plot.area(ax=axes[1], color=colors, linewidth=0)
    axes[1].set_title(f"Dispatch at {bus} - {outage_start.date()} → {outage_end.date()}")
    axes[1].set_xlabel("Time")
    axes[1].grid(alpha=0.3)

    # Remove subplot legends
    axes[0].legend().remove()
    axes[1].legend().remove()

    # Shared legend at bottom
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5)

    plt.tight_layout(rect=[0, 0.08, 1, 1])

### dispatch 2 - per area
def plot_dispatch_for_area(n, area_prefix, outage_start, outage_end):
 

    # Find all buses belonging to this area 
    buses = n.buses.index[n.buses.index.str.startswith(area_prefix)]
    if len(buses) == 0:
        print(f"No buses found starting with '{area_prefix}'")
        return
    
    print(f"Found {len(buses)} buses in area '{area_prefix}'")

    # Select all generators + stores located on these buses
    gens = n.generators[n.generators.bus.isin(buses)]
    stores = n.storage_units[n.storage_units.bus.isin(buses)]

    gen_names = gens.index.tolist()
    store_names = stores.index.tolist()

    print(f"Generators in area: {len(gen_names)}")
    print(f"Storage units in area: {len(store_names)}")

    # Extract dispatch time series 
    dispatch_g = n.generators_t.p[gen_names] if gen_names else pd.DataFrame(index=n.snapshots)
    dispatch_s = n.storage_units_t.p[store_names] if store_names else pd.DataFrame(index=n.snapshots)

    # Storage discharging only 
    if not dispatch_s.empty:
        dispatch_s = dispatch_s.clip(lower=0)

    dispatch = pd.concat([dispatch_g, dispatch_s], axis=1)

    # Group by carrier 
    carriers = pd.concat([gens.carrier, stores.carrier])
    dispatch_by_carrier = dispatch.groupby(carriers, axis=1).sum()

    # Extract outage period 
    outage_start = pd.to_datetime(outage_start)
    outage_end = pd.to_datetime(outage_end)
    dispatch_outage = dispatch_by_carrier.loc[outage_start:outage_end]

    # Plot 
    colors = [get_color(c) for c in dispatch_by_carrier.columns]

    fig, axes = plt.subplots(1, 2, figsize=(22, 6), sharey=True)

    # Full year
    dispatch_by_carrier.plot.area(ax=axes[0], color=colors, linewidth=0)
    axes[0].set_title(f"Dispatch - {area_prefix} - Full Period")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("MW")
    axes[0].grid(alpha=0.3)

    # Outage
    dispatch_outage.plot.area(ax=axes[1], color=colors, linewidth=0)
    axes[1].set_title(f"Dispatch - {area_prefix} - Outage Period")
    axes[1].set_xlabel("Time")
    axes[1].grid(alpha=0.3)

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5)

    plt.tight_layout(rect=[0, 0.08, 1, 1])


### More dispatch per area (incl. stores)
def plot_dispatch_for_area_2(n, area_prefix, outage_start, outage_end):

    buses = n.buses.index[n.buses.index.str.startswith(area_prefix)]
    if buses.empty:
        print(f"No buses found starting with '{area_prefix}'")
        return

    # Components 
    gens = n.generators[n.generators.bus.isin(buses)]
    sus = n.storage_units[n.storage_units.bus.isin(buses)]
    stores = n.stores[n.stores.bus.isin(buses)]

    print(f"Generators: {len(gens)}")
    print(f"Storage units: {len(sus)}")
    print(f"Stores: {len(stores)}")

    # Generators
    gen_p = n.generators_t.p[gens.index] if not gens.empty else pd.DataFrame(index=n.snapshots)
    gen_carriers = gens.carrier

    # Storage units
    if not sus.empty:
        su_p = n.storage_units_t.p[sus.index]
        su_dis = su_p.clip(lower=0)
        su_ch = su_p.clip(upper=0)

        su_dis.columns = [f"{c} discharge" for c in sus.carrier]
        su_ch.columns = [f"{c} charge" for c in sus.carrier]
    else:
        su_dis = su_ch = pd.DataFrame(index=n.snapshots)

    # Stores
    if not stores.empty:
        st_p = n.stores_t.p[stores.index]
        st_dis = st_p.clip(lower=0)
        st_ch = st_p.clip(upper=0)

        st_dis.columns = [f"{c} discharge" for c in stores.carrier]
        st_ch.columns = [f"{c} charge" for c in stores.carrier]
    else:
        st_dis = st_ch = pd.DataFrame(index=n.snapshots)

    # Combine
    dispatch = pd.concat([gen_p, su_dis, su_ch, st_dis, st_ch], axis=1)

    carriers = pd.concat([
        gen_carriers,
        pd.Series([f"{c} discharge" for c in sus.carrier], index=su_dis.columns),
        pd.Series([f"{c} charge" for c in sus.carrier], index=su_ch.columns),
        pd.Series([f"{c} discharge" for c in stores.carrier], index=st_dis.columns),
        pd.Series([f"{c} charge" for c in stores.carrier], index=st_ch.columns),
    ])

    dispatch_by_carrier = dispatch.groupby(carriers, axis=1).sum()

    # Disruption window
    outage_start = pd.to_datetime(outage_start)
    outage_end = pd.to_datetime(outage_end)
    dispatch_outage = dispatch_by_carrier.loc[outage_start:outage_end]

    # Plotting
    colors = [
        get_color(c.replace(" discharge", "").replace(" charge", ""))
        for c in dispatch_by_carrier.columns
    ]

    fig, axes = plt.subplots(1, 2, figsize=(24, 9), sharey=True)

    dispatch_by_carrier.plot.area(
        ax=axes[0], color=colors, linewidth=0, legend=False
    )
    axes[0].set_title(f"Dispatch - {area_prefix} - Full Period")
    axes[0].set_ylabel("MW")
    axes[0].grid(alpha=0.3)

    dispatch_outage.plot.area(
        ax=axes[1], color=colors, linewidth=0, legend=False
    )
    axes[1].set_title(f"Dispatch - {area_prefix} - Outage Period")
    axes[1].grid(alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=10, frameon=False)

    plt.tight_layout(rect=[0, 0.18, 1, 1])


### plot load shedding (electricity)
def plot_load_shedding(n, start, end):
    load_gens = n.generators[n.generators.carrier == "elec load shed"].index
    if load_gens.empty:
        print("No load shedding generators found.")
        return

    total_shed = n.generators_t.p[load_gens].sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(18, 4), sharey=True)

    axes[0].plot(total_shed.index, total_shed)
    axes[0].set_title("Load Shedding - Full Year")

    mask = (total_shed.index >= start) & (total_shed.index <= end)
    axes[1].plot(total_shed.index[mask], total_shed[mask])
    axes[1].set_title("Load Shedding - Outage Window")

    plt.tight_layout()

### plot load shedding per bus (electricity)
def plot_load_shedding_per_bus(n, start, end):
    
    load_gens = n.generators[n.generators.carrier == "elec load shed"]

    if load_gens.empty:
        print("No load shedding generators found.")
        return

    gen_names = load_gens.index
    buses = load_gens.bus
    load_shed = n.generators_t.p[gen_names]

    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    # Full year
    for gen, bus in zip(gen_names, buses):
        axes[0].plot(
            n.snapshots,
            load_shed[gen],
            label=f"{bus}",
        )

    axes[0].set_title("Load Shedding per Bus - Full Year")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("MW")
    axes[0].grid(True, alpha=0.3)
    axes[0].axvspan(start, end, color="gray", alpha=0.2)
    axes[0].legend(loc="upper left", fontsize=8)

    # Zoom 
    for gen, bus in zip(gen_names, buses):
        axes[1].plot(
            n.snapshots[mask],
            load_shed[gen][mask],
            label=f"{bus}"
        )

    axes[1].set_title(f"Load Shedding per Bus - {start.date()} - {end.date()}")
    axes[1].set_xlabel("Time")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()

### plot biomass load shedding
def plot_biomass_shedding_per_bus(n, start, end, areas):

    # Filter biomass shedding generators by area
    bio_gens = n.generators[
        (n.generators.carrier == "biomass shed")
        & (n.generators.bus.str[:3].isin(areas))
    ]

    if bio_gens.empty:
        print("No biomass shedding generators found for selected areas.")
        return

    gen_names = bio_gens.index
    buses = bio_gens.bus
    bio_shed = n.generators_t.p[gen_names]

    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    fig, axes = plt.subplots(1, 2, figsize=(20, 6), sharey=True)

    # Full period
    for gen, bus in zip(gen_names, buses):
        axes[0].plot(
            n.snapshots,
            bio_shed[gen],
            label=bus,
        )

    axes[0].set_title("Biomass Shedding per Bus - Full Period")
    axes[0].set_ylabel("MW")
    axes[0].grid(True, alpha=0.3)
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # Outage window
    for gen, bus in zip(gen_names, buses):
        axes[1].plot(
            n.snapshots[mask],
            bio_shed[gen][mask],
            label=bus
        )

    axes[1].set_title(f"Biomass Shedding per Bus - {start.date()} - {end.date()}")
    axes[1].grid(True, alpha=0.3)

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        fontsize=8
    )

    plt.tight_layout(rect=[0, 0.15, 1, 1])

### plot heat shedding per bus
def heat_shedding_per_bus(n, start, end, areas):

    # Filter heat shedding generators by area
    heat_gens = n.generators[
        (n.generators.carrier == "heat shed")
        & (n.generators.bus.str[:3].isin(areas))
    ]

    if heat_gens.empty:
        print("No heat shedding generators found for selected areas.")
        return

    gen_names = heat_gens.index
    buses = heat_gens.bus
    heat_shed = n.generators_t.p[gen_names]

    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    fig, axes = plt.subplots(1, 2, figsize=(20, 6), sharey=True)

    # Full period
    for gen, bus in zip(gen_names, buses):
        axes[0].plot(
            n.snapshots,
            heat_shed[gen],
            label=bus,
        )

    axes[0].set_title("Heat Shedding per Bus - Full Period")
    axes[0].set_ylabel("MW")
    axes[0].grid(True, alpha=0.3)
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # Outage window
    for gen, bus in zip(gen_names, buses):
        axes[1].plot(
            n.snapshots[mask],
            heat_shed[gen][mask],
            label=bus
        )

    axes[1].set_title(f"Heat Shedding per Bus - {start.date()} - {end.date()}")
    axes[1].grid(True, alpha=0.3)

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        fontsize=8
    )

    plt.tight_layout(rect=[0, 0.15, 1, 1])


######################## ELECTRICITY PRICES #########################

def plot_elec_prices_full_and_zoom(n, start, end, outage_mask=None):

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    # Select electricity buses
    elec_carriers = ["AC", "DC", "low voltage"]
    elec_buses = n.buses.index[n.buses.carrier.isin(elec_carriers)]

    prices = n.buses_t.marginal_price[elec_buses]

    fig, axes = plt.subplots(1, 2, figsize=(18, 5), sharey=True)

    # ---- Full year ----
    prices.plot(ax=axes[0], linewidth=0.8)
    axes[0].set_title("Electricity prices - Full Year")
    axes[0].set_ylabel("€/MWh")
    axes[0].grid(alpha=0.3)

    # Shade outage window
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # ---- Zoom ----
    zoom = prices.loc[start:end]
    zoom.plot(ax=axes[1], linewidth=0.8)
    axes[1].set_title(f"Electricity prices - {start.date()} - {end.date()}")
    axes[1].set_ylabel("€/MWh")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()


def plot_heat_prices_full_and_zoom(n, start, end, outage_mask=None):

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    # Select heat buses
    heat_carriers = ["rural heat", "urban central heat", "urban decentral heat"]
    heat_buses = n.buses.index[n.buses.carrier.isin(heat_carriers)]

    if len(heat_buses) == 0:
        print("No heat buses found - skipping heat price plot.")
        return

    prices = n.buses_t.marginal_price[heat_buses]

    fig, axes = plt.subplots(1, 2, figsize=(18, 5), sharey=True)

    # Full year 
    prices.plot(ax=axes[0], linewidth=0.8)
    axes[0].set_title("Heat prices - Full Year")
    axes[0].set_ylabel("€/MWh")
    axes[0].grid(alpha=0.3)

    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # Zoom 
    zoom = prices.loc[start:end]
    zoom.plot(ax=axes[1], linewidth=0.8)
    axes[1].set_title(f"Heat prices - {start.date()} - {end.date()}")
    axes[1].set_ylabel("€/MWh")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()


### plot electricity price per country
def plot_price_elec_by_country(n, country, start, end):

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    # Select electricity buses in this country
    elec_carriers = ["AC", "low voltage"]
    mask_buses = n.buses.carrier.isin(elec_carriers) & (n.buses.country == country)
    buses = n.buses.index[mask_buses]

    if len(buses) == 0:
        print(f"No electricity buses found in country {country}.")
        return

    prices = n.buses_t.marginal_price[buses]
    if prices.empty:
        print(f"No marginal prices for electricity buses in {country}.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    # Full year
    prices.plot(ax=axes[0], legend=True)
    axes[0].set_title(f"Electricity prices - {country} - Full Year")
    axes[0].set_ylabel("€/MWh")
    axes[0].axvspan(start, end, color="gray", alpha=0.2)
    axes[0].grid(alpha=0.3)

    # Zoom
    zoom = prices.loc[start:end]
    zoom.plot(ax=axes[1], legend=False)
    axes[1].set_title(f"Electricity prices - {country} - {start.date()} - {end.date()}")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()

def plot_price_heat_by_country(n, country, start, end):

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    heat_carriers = ["rural heat", "urban central heat", "urban decentral heat"]
    mask_buses = n.buses.carrier.isin(heat_carriers) & (n.buses.country == country)
    buses = n.buses.index[mask_buses]

    if len(buses) == 0:
        print(f"No heat buses found in country {country}.")
        return

    prices = n.buses_t.marginal_price[buses]
    if prices.empty:
        print(f"No marginal prices for heat buses in {country}.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    # Full year
    prices.plot(ax=axes[0], legend=True)
    axes[0].set_title(f"Heat prices - {country} - Full Year")
    axes[0].set_ylabel("€/MWh")
    axes[0].axvspan(start, end, color="gray", alpha=0.2)
    axes[0].grid(alpha=0.3)

    # Zoom
    zoom = prices.loc[start:end]
    zoom.plot(ax=axes[1], legend=False)
    axes[1].set_title(f"Heat prices - {country} - {start.date()} - {end.date()}")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()


########################## POWER FLOWS #########################
def get_flows_at_bus(n, bus):
    records = []

    # AC LINES 
    for line, row in n.lines.iterrows():
        if row.bus0 == bus:     # flow OUT of bus
            flow = n.lines_t.p0[line]
            direction = 1       # export
            other = row.bus1
            cap  = row.s_nom  # MW
        elif row.bus1 == bus:   # flow INTO bus
            flow = n.lines_t.p1[line]
            direction = -1      # import
            other = row.bus0
            cap  = row.s_nom  # MW
        else:
            continue

        name = f"{bus} ↔ {other} (AC {line}, {cap:.0f} MW)"
        records.append(pd.DataFrame({
            "timestamp": flow.index,
            "element": name,
            "type": "AC",
            "flow": flow * direction     # positive = export, negative = import
        }))

    # DC LINKS 
    for link, row in n.links.iterrows():
        if row.bus0 == bus:     # flow OUT from bus0
            flow = n.links_t.p0[link]
            direction = 1
            other = row.bus1
            cap  = row.p_nom  # MW
        elif row.bus1 == bus:   # flow INTO bus1
            flow = n.links_t.p1[link]
            direction = -1
            other = row.bus0
            cap  = row.p_nom  # MW
        else:
            continue

        name = f"{bus} ↔ {other} (DC {link}, {cap:.0f} MW)"
        records.append(pd.DataFrame({
            "timestamp": flow.index,
            "element": name,
            "type": "DC",
            "flow": flow * direction
        }))

    if not records:
        return pd.DataFrame()

    df = pd.concat(records).set_index("timestamp")
    return df.sort_index()

def plot_flows_at_bus(n, bus, start=None, end=None):
    df = get_flows_at_bus(n, bus)

    if start:
        df = df.loc[pd.to_datetime(start):]
    if end:
        df = df.loc[:pd.to_datetime(end)]

    plt.figure(figsize=(18, 6))

    for element in df["element"].unique():
        sub = df[df.element == element]
        plt.plot(sub.index, sub["flow"], label=element)

    plt.axhline(0, color="black", linewidth=0.8)
    plt.title(f"Power Flows at {bus}")
    plt.ylabel("Flow (MW)")
    plt.legend(loc="upper left")
    plt.grid(alpha=0.3)
    plt.tight_layout()

### plot biomass dispatvh per area
def plot_bio_dispatch_for_area(n, area_prefix, outage_start, outage_end):

    # Buses in this area
    buses = n.buses.index[n.buses.index.str.startswith(area_prefix)]
    if len(buses) == 0:
        print(f"No buses found starting with '{area_prefix}'")
        return

    # Biomass-related generators & stores on those buses
    bio_keywords = ["bio", "biomass", "biogas"]

    gen_mask_area = n.generators.bus.isin(buses)
    gen_mask_bio = n.generators.carrier.str.contains(
        "|".join(bio_keywords), case=False, na=False
    )
    gens = n.generators[gen_mask_area & gen_mask_bio]

    store_mask_area = n.storage_units.bus.isin(buses)
    store_mask_bio = n.storage_units.carrier.str.contains(
        "|".join(bio_keywords), case=False, na=False
    )
    stores = n.storage_units[store_mask_area & store_mask_bio]

    gen_names = gens.index.tolist()
    store_names = stores.index.tolist()

    if len(gen_names) + len(store_names) == 0:
        print(f"No biomass-related generators or storage in area '{area_prefix}'.")
        return

    print(f"[{area_prefix}] Biomass gens: {len(gen_names)}, stores: {len(store_names)}")

    # Extract dispatch time series 
    dispatch_g = (
        n.generators_t.p[gen_names]
        if gen_names
        else pd.DataFrame(index=n.snapshots)
    )
    dispatch_s = (
        n.storage_units_t.p[store_names]
        if store_names
        else pd.DataFrame(index=n.snapshots)
    )

    # Only count discharging from storage
    if not dispatch_s.empty:
        dispatch_s = dispatch_s.clip(lower=0)

    dispatch = pd.concat([dispatch_g, dispatch_s], axis=1)

    carriers = pd.concat([gens.carrier, stores.carrier])
    dispatch_by_carrier = dispatch.groupby(carriers, axis=1).sum()

    # outage window 
    outage_start = pd.to_datetime(outage_start)
    outage_end = pd.to_datetime(outage_end)
    dispatch_outage = dispatch_by_carrier.loc[outage_start:outage_end]

    # Plot
    colors = [get_color(c) for c in dispatch_by_carrier.columns]

    fig, axes = plt.subplots(1, 2, figsize=(22, 6), sharey=True)

    # Full period
    dispatch_by_carrier.plot.area(ax=axes[0], color=colors, linewidth=0)
    axes[0].set_title(f"Biomass dispatch - {area_prefix} - Full Period")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("MW")
    axes[0].grid(alpha=0.3)

    # Outage period
    dispatch_outage.plot.area(ax=axes[1], color=colors, linewidth=0)
    axes[1].set_title(
        f"Biomass dispatch - {area_prefix} - {outage_start.date()} → {outage_end.date()}"
    )
    axes[1].set_xlabel("Time")
    axes[1].grid(alpha=0.3)

    # Shared legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4)

    plt.tight_layout(rect=[0, 0.08, 1, 1])


def plot_biomass_link_flows(n, areas=None, start=None, end=None):


    if isinstance(areas, str):
        areas = [areas]

    # Select biomass links
    bio_links = n.links[n.links.carrier.str.contains("bio", case=False, na=False)]

    if areas is not None:
        mask = bio_links.apply(
            lambda row: any(row.bus0.startswith(a) or row.bus1.startswith(a) for a in areas),
            axis=1
        )
        bio_links = bio_links[mask]

    if bio_links.empty:
        print("No biomass-related links found.")
        return None

    # Extract flows 
    flows = {}
    for lk in bio_links.index:
        # Use positive flow into bus1 (absolute biomass use)
        flows[lk] = n.links_t.p0[lk]

    df = pd.DataFrame(flows)

    # Time filtering
    if start:
        df = df.loc[pd.to_datetime(start):]
    if end:
        df = df.loc[:pd.to_datetime(end)]

    # Plot 
    fig, ax = plt.subplots(figsize=(18, 6))
    df.plot(ax=ax, linewidth=0.8)

    title_area = ", ".join(areas) if areas else "All areas"
    ax.set_title(f"Biomass link flows - {title_area}")
    ax.set_ylabel("MW (flow p0)")
    ax.grid(alpha=0.3)
    plt.tight_layout()


def plot_stacked_biomass_shedding_per_bus(n, start, end, areas):

    bio_gens = n.generators[
        (n.generators.carrier == "biomass shed")
        & (n.generators.bus.str[:3].isin(areas))
    ]

    if bio_gens.empty:
        print("No biomass shedding generators found.")
        return

    # Power dispatch (MW)
    p = n.generators_t.p[bio_gens.index]

    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    p = p.loc[start:end]

    w = n.snapshot_weightings.generators.loc[p.index]  
    energy_MWh_per_gen = p.mul(w, axis=0).sum()

    df = (
        pd.DataFrame({
            "bus": bio_gens.bus.values,
            "MWh shed": energy_MWh_per_gen.values
        })
        .groupby("bus")
        .sum()
        .sort_values("MWh shed", ascending=False)
    )

    # Plot
    fig, ax = plt.subplots(figsize=(10, 5))
    bars = df["MWh shed"].plot.bar(ax=ax, color="darkolivegreen")

    ax.set_ylabel("MWh")
    ax.set_title("Total Biomass Shedding per Bus")
    ax.grid(axis="y", alpha=0.3)

    for bar in bars.patches:
        height = bar.get_height()
        if height > 0:
            ax.annotate(
                f"{height:.1f}",
                (bar.get_x() + bar.get_width() / 2, height),
                ha="center",
                va="bottom",
                fontsize=9,
                xytext=(0, 3),
                textcoords="offset points",
            )

    plt.tight_layout()

# heat shedding stacked
def plot_stacked_heat_shedding_per_bus(n, start, end, areas):

    heat_gens = n.generators[
        (n.generators.carrier == "heat shed")
        & (n.generators.bus.str[:3].isin(areas))
    ]

    if heat_gens.empty:
        print("No heat shedding generators found.")
        return

    p = n.generators_t.p[heat_gens.index]

    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    p = p.loc[start:end]

    w = n.snapshot_weightings.generators.loc[p.index]
    energy_MWh_per_gen = p.mul(w, axis=0).sum()

    df = (
        pd.DataFrame({
            "bus": heat_gens.bus.values,
            "MWh shed": energy_MWh_per_gen.values
        })
        .groupby("bus")
        .sum()
        .sort_values("MWh shed", ascending=False)
    )

    fig, ax = plt.subplots(figsize=(10, 5))
    bars = df["MWh shed"].plot.bar(ax=ax, color="darkgreen")

    ax.set_ylabel("MWh")
    ax.set_title(f"Total Heat Shedding per Bus ({', '.join(areas)})")
    ax.grid(axis="y", alpha=0.3)

    for bar in bars.patches:
        height = bar.get_height()
        if height > 0:
            ax.annotate(
                f"{height:.1f}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # offset above bar
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9
            )

    plt.tight_layout()


############################################################
# SNAKEMAKE ENTRY POINT

if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.net)
    scenario = snakemake.wildcards.scenario

    outdir = f"results_contingencies/2035/bio/{scenario}/plots"
    os.makedirs(outdir, exist_ok=True)

    start = pd.to_datetime(snakemake.params.start)
    duration = pd.Timedelta(days=snakemake.params.duration)
    end = start + duration
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    areas=snakemake.params.areas
    scaling_factor=snakemake.params.scaling_factor

    print(f"Generating plots for scenario {scenario}…")

    # Bio links under shortage
    plot_biomass_link_capacity(n, areas, start, end)
    plt.savefig(f"{outdir}/bio_links_status.png")
    plt.close()

    # bio dispatch
    areas_to_plot = ["DK0", "DK1", "SE1", "DE0", "GB2", "GB3"]
    for area in areas_to_plot:
        plot_bio_dispatch_for_area(n, area, start, end)
        plt.savefig(f"{outdir}/bio_dispatch_{area}.png")
        plt.close()
    
    # bio link flows
    areas_to_plot = ["DK0", "DK1", "SE1", "DE0"]
    for area in areas_to_plot:
        plot_biomass_link_flows(n, areas=[area], start=start, end=end)
        plt.savefig(f"{outdir}/bio_link_flows_{area}.png")
        plt.close()

    # Dispatch by carrier
    plot_dispatch_by_carrier(n, start, end)
    plt.savefig(f"{outdir}/dispatch_by_carrier.png")
    plt.close()

    # Dispatch by carrier at specific buses
    buses_to_plot = ["DK0 0", "DK1 0", "SE1 0", "DE0 0"]
    for bus in buses_to_plot:
        plot_bus_dispatch_by_carrier(n, bus, start, end)
        plt.savefig(f"{outdir}/dispatch_by_carrier_{bus.replace(' ','_')}.png")
        plt.close()
    
    # Dispatch by carrier for specific areas 
    areas_to_plot = ["DK0", "DK1", "SE1", "DE0"]
    for area in areas_to_plot:
        plot_dispatch_for_area(n, area, start, end)
        plt.savefig(f"{outdir}/dispatch_by_carrier_area_{area}.png")
        plt.close()
    
    # Dispatch by carrier for specific areas
    areas_to_plot = ["DK0 0", "DK1 0", "SE1 0", "DE0 0", "GB2 0", "GB3 0", "FI1 0"]
    for area in areas_to_plot:
        plot_dispatch_for_area_2(n, area, start, end)
        plt.savefig(f"{outdir}/dispatch_by_carrier_area_2_{area}.png")
        plt.close()

    # Load shedding per bus
    plot_load_shedding_per_bus(n, start, end)
    plt.savefig(f"{outdir}/load_shed_per_bus.png")
    plt.close()

    # Biomass shedding per bus
    plot_biomass_shedding_per_bus(n, start, end, areas)
    plt.savefig(f"{outdir}/biomass_shed_per_bus.png")
    plt.close()

    # Biomass shedding per bus STACKED
    plot_stacked_biomass_shedding_per_bus(n, start, end, areas)
    plt.savefig(f"{outdir}/biomass_shed_stacked_per_bus.png")
    plt.close()

    #Heat shedding per bus
    heat_shedding_per_bus(n, start, end, areas)
    plt.savefig(f"{outdir}/heat_shed_per_bus.png")
    plt.close()

    #Heat shedding per bus STACKED
    plot_stacked_heat_shedding_per_bus(n, start, end, areas)
    plt.savefig(f"{outdir}/heat_shed_stacked_per_bus.png")
    plt.close()

    # Electricity Prices
    plot_elec_prices_full_and_zoom(n, start, end, outage_mask=mask)
    plt.savefig(f"{outdir}/prices_elec_plot.png")
    plt.close()

    #Heat prices
    plot_heat_prices_full_and_zoom(n, start, end, outage_mask=mask)
    plt.savefig(f"{outdir}/prices_heat_plot.png")
    plt.close()

    #Electricity prices per bus
    countries_to_plot = ["DK", "SE", "DE"]
    for c in countries_to_plot:
        plot_price_elec_by_country(n, c, start, end)
        plt.savefig(f"{outdir}/electricity_price_{c}.png")
        plt.close()

    #Heat prices per bus
    countries_to_plot = ["DK", "SE", "DE"]
    for c in countries_to_plot:
        plot_price_heat_by_country(n, c, start, end)
        plt.savefig(f"{outdir}/heat_price_{c}.png")
        plt.close()

    # Flows for selected buses
    buses_to_plot = ["DK0 0", "DK1 0", "SE1 0", "DE0 0"]
    for bus in buses_to_plot:
        plot_flows_at_bus(n, bus, start, end)
        plt.savefig(f"{outdir}/flows_{bus.replace(' ','_')}.png")
        plt.close()

    with open(snakemake.output[0], "w") as f:
        f.write("OK\n")

    print("All plots created successfully!")