# plot_trans_2035.py
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
    'heat shed': '#ff0009'
    
}

def get_color(carrier):
    return tech_colors.get(carrier, "#777777")


############################################################
# PLOTTING FUNCTIONS

def plot_outage_status(n, lines=None, links=None, start=None, end=None, outdir="plots"):

    os.makedirs(outdir, exist_ok=True)

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    # Helper for month formatting
    def format_month_axis(ax):
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.grid(True, alpha=0.3)

    # AC LINES 
    if lines:
        plt.figure(figsize=(14, 4))

        for ln in lines:
            if ln not in n.lines.index:
                print(f"⚠ Line '{ln}' not found.")
                continue

            series = n.lines_t.s_max_pu[ln]
            plt.plot(series.index, series.values, label=ln)

        plt.title("AC Line availability (s_max_pu)")
        plt.ylabel("s_max_pu")
        plt.axvspan(start, end, color="red", alpha=0.15)
        plt.legend(title="Line")
        ax = plt.gca()
        format_month_axis(ax)

        plt.tight_layout()
        plt.savefig(f"{outdir}/ac_lines_status.png")
        plt.close()

    # DC LINKS 
    if links:
        plt.figure(figsize=(14, 4))

        for lk in links:
            if lk not in n.links.index:
                print(f"⚠ Link '{lk}' not found.")
                continue

            series = n.links_t.p_max_pu[lk]
            plt.plot(series.index, series.values, label=lk)

        plt.title("DC Link availability (p_max_pu)")
        plt.ylabel("p_max_pu")
        plt.axvspan(start, end, color="red", alpha=0.15)
        plt.legend(title="Link")
        ax = plt.gca()
        format_month_axis(ax)

        plt.tight_layout()
        plt.savefig(f"{outdir}/dc_links_status.png")
        plt.close()

### dispatch plots 
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
    dispatch_by_carrier.plot.area(ax=axes[0], color=colors, linewidth=0, legend=False)
    axes[0].set_title(f"Dispatch at {bus} – Full Year")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("MW")
    axes[0].grid(alpha=0.3)

    # Outage period
    dispatch_outage.plot.area(ax=axes[1], color=colors, linewidth=0, legend=False)
    axes[1].set_title(f"Dispatch at {bus} – {outage_start.date()} → {outage_end.date()}")
    axes[1].set_xlabel("Time")
    axes[1].grid(alpha=0.3)

    # Remove subplot legends
    axes[0].legend().remove()
    axes[1].legend().remove()

    # Shared legend at bottom
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5)

    plt.tight_layout(rect=[0, 0.08, 1, 1])

### dispatch per area
def plot_dispatch_for_area(n, area_prefix, outage_start, outage_end):
  

    # buses belonging to this area 
    buses = n.buses.index[n.buses.index.str.startswith(area_prefix)]
    if len(buses) == 0:
        print(f"No buses found starting with '{area_prefix}'")
        return
    
    print(f"Found {len(buses)} buses in area '{area_prefix}'")

    # all generators + stores located on these buses
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

### dispatch plot per area (incl stores)
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

    # Disruption period
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
    fig.legend(handles, labels, loc="lower center", ncol=6)

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

# plot load shedding per bus (electricity)
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

    axes[1].set_title(f"Load Shedding per Bus - {start.date()} → {end.date()}")
    axes[1].set_xlabel("Time")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()


### plot load shedding DK2 2035
def plot_load_shedding_DK2_2035(n, start, end):

    # Identify DK2 buses
    dk2_buses = n.buses.index[n.buses.index.str.startswith("DK1")]

    # Electricity load shedding
    elec_ls = n.generators[
        (n.generators.carrier == "elec load shed") &
        (n.generators.bus.isin(dk2_buses))
    ]

    elec_ts = (
        n.generators_t.p[elec_ls.index].sum(axis=1)
        if not elec_ls.empty
        else pd.Series(0.0, index=n.snapshots)
    )

    # Heat load shedding
    heat_ls = n.generators[
        (n.generators.carrier == "heat shed") &
        (n.generators.bus.isin(dk2_buses))
    ]

    heat_ts = (
        n.generators_t.p[heat_ls.index].sum(axis=1)
        if not heat_ls.empty
        else pd.Series(0.0, index=n.snapshots)
    )

    # Time window
    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(18, 4), sharey=True)

    # Full period 
    axes[0].plot(elec_ts.index, elec_ts, label="Electric load shed", color="red")
    axes[0].plot(heat_ts.index, heat_ts, label="Heat load shed", color="darkorange")
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    axes[0].set_title("DK2 Load Shedding - Full Period")
    axes[0].set_ylabel("MW")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    # Outage window
    axes[1].plot(
        elec_ts.loc[mask].index,
        elec_ts.loc[mask],
        label="Electric load shed",
        color="red",
    )
    axes[1].plot(
        heat_ts.loc[mask].index,
        heat_ts.loc[mask],
        label="Heat load shed",
        color="darkorange",
    )

    axes[1].set_title("DK2 Load Shedding - Outage Window")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    plt.tight_layout()

######################## ELECTRICITY PRICES #########################

def plot_elec_prices_full_and_zoom(n, start, end, outage_mask=None):

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    # electricity buses
    elec_carriers = ["AC", "DC", "low voltage"]
    elec_buses = n.buses.index[n.buses.carrier.isin(elec_carriers)]

    prices = n.buses_t.marginal_price[elec_buses]

    fig, axes = plt.subplots(1, 2, figsize=(18, 5), sharey=True)

    # Full year 
    prices.plot(ax=axes[0], linewidth=0.8)
    axes[0].set_title("Electricity prices - Full Year")
    axes[0].set_ylabel("€/MWh")
    axes[0].grid(alpha=0.3)

    # Shade outage window
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # Zoom 
    zoom = prices.loc[start:end]
    zoom.plot(ax=axes[1], linewidth=0.8)
    axes[1].set_title(f"Electricity prices - {start.date()} - {end.date()}")
    axes[1].set_ylabel("€/MWh")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()

## plot heat prices
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
    axes[1].set_title(f"Heat prices - {start.date()} → {end.date()}")
    axes[1].set_ylabel("€/MWh")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()



### plot electricity prices by country
def plot_price_elec_by_country(n, country, start, end):

    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)

    #  electricity buses in this country
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

### plot heat prices by country
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

################################## POWER FLOWS ##############################
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




############################################################
# SNAKEMAKE ENTRY POINT

if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.net)
    scenario = snakemake.wildcards.scenario

    outdir = f"results_contingencies/2035/transmission/{scenario}/plots"
    os.makedirs(outdir, exist_ok=True)

    start = pd.to_datetime(snakemake.params.start)
    duration = pd.Timedelta(days=snakemake.params.duration)
    end = start + duration
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    lines=snakemake.params.lines
    links=snakemake.params.links

    print(f"Generating plots for scenario {scenario}…")

    # Links and lines under outage
    plot_outage_status(n,lines,links,start=start,end=end, outdir=outdir)
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
    areas_to_plot = ["DK0", "DK1", "SE1", "DE0"]
    for area in areas_to_plot:
        plot_dispatch_for_area_2(n, area, start, end)
        plt.savefig(f"{outdir}/dispatch_by_carrier_area_2_{area}.png")
        plt.close()

    # load shed DK2 2035
    plot_load_shedding_DK2_2035(n, start, end)
    plt.savefig(f"{outdir}/load_shedding_DK2_2035.png")
    plt.close()

    # Load shedding per bus
    plot_load_shedding_per_bus(n, start, end)
    plt.savefig(f"{outdir}/load_shed_per_bus.png")
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