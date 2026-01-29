# plot_bio_2025.py
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
    'PHS': '#51dbcc',
    'battery': '#ace37f',
    
    'elec load shed': '#ff0000',
    'heat shed': '#4F000B'
    
}

def get_color(carrier):
    return tech_colors.get(carrier, "#777777")


############################################################
# PLOTTING FUNCTIONS

### plot bio availability 
def plot_affected_biomass_generators(n, buses, start, end):

    # Locate biomass generators
    bio_mask = n.generators.carrier.str.contains("bio", case=False, na=False)
    bio_mask &= n.generators.bus.isin(buses)

    biomass_gens = n.generators.index[bio_mask]

    if biomass_gens.empty:
        print(f"No biomass generators at buses {buses}.")
        return

    print("Biomass generators at selected buses:")
    for g in biomass_gens:
        print(f"  - {g}")

    # Ensure p_max_pu columns exist for plotting
    for gen in biomass_gens:
        if gen not in n.generators_t.p_max_pu.columns:
            n.generators_t.p_max_pu[gen] = 1.0
            n.generators_t.p_min_pu[gen] = 0.0

    # Time
    start = pd.to_datetime(start)
    end   = pd.to_datetime(end)
    mask  = (n.snapshots >= start) & (n.snapshots <= end)

    # Extract p_max_pu time series
    pmax = n.generators_t.p_max_pu[biomass_gens]

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    # FULL YEAR
    pmax.plot(ax=axes[0], linewidth=1)
    axes[0].set_title("Biomass Generator Availability - Full Year")
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel("p_max_pu")
    axes[0].grid(alpha=0.3)

    # Highlight outage period
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # zoomed outage period 
    pmax.loc[mask].plot(ax=axes[1], linewidth=1)
    axes[1].set_title(f"Zoom: {start.date()} → {end.date()}")
    axes[1].set_xlabel("Time")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()



### dispatch by carrier
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

    dispatch_by_carrier.plot.area(ax=axes[0], color=colors, linewidth=0)
    axes[0].set_title("Full-Year Dispatch by Carrier")

    dispatch_outage.plot.area(ax=axes[1], color=colors, linewidth=0)
    axes[1].set_title("Outage-Period Dispatch by Carrier")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5)

    plt.tight_layout(rect=[0, 0.08, 1, 1])

### plot dispatch by carrier per bus
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
    axes[1].set_title(f"Dispatch at {bus} - {outage_start.date()} - {outage_end.date()}")
    axes[1].set_xlabel("Time")
    axes[1].grid(alpha=0.3)

    # Remove subplot legends
    axes[0].legend().remove()
    axes[1].legend().remove()

    # Shared legend at bottom
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5)

    plt.tight_layout(rect=[0, 0.08, 1, 1])

## plot load shedding (electricity)
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

### plot load shedding in DK2 (electricity)
def plot_load_shedding_DK2(n, start, end, bus="DK1 0"):

    # Find electric load shedding generators at this bus
    ls_gens = n.generators[
        (n.generators.carrier == "elec load shed") &
        (n.generators.bus == bus)
    ]

    if ls_gens.empty:
        print(f"No electricity load shedding at bus {bus}.")
        return

    # Aggregate load shedding time series
    ls_ts = n.generators_t.p[ls_gens.index].sum(axis=1)

    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    mask = (ls_ts.index >= start) & (ls_ts.index <= end)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(18, 4), sharey=True)

    # Full year
    axes[0].plot(ls_ts.index, ls_ts, color="red", linewidth=1.5)
    axes[0].set_title(f"Electric Load Shedding at {bus} - Full Year")
    axes[0].set_ylabel("MW")
    axes[0].grid(True, alpha=0.3)
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # Outage window
    axes[1].plot(ls_ts.index[mask], ls_ts[mask], color="red", linewidth=1.5)
    axes[1].set_title(f"Electric Load Shedding at {bus} - Outage Window")
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()

######################## ELECTRICITY PRICES #########################
def plot_prices_full_and_zoom(n, start, end, outage_mask=None):
    prices = n.buses_t.marginal_price
    fig, axes = plt.subplots(1, 2, figsize=(18, 5), sharey=True)

    prices.plot(ax=axes[0])
    axes[0].set_title("Prices - Full Year")

    zoom = prices.loc[start:end]
    zoom.plot(ax=axes[1])
    axes[1].set_title(f"Prices {start.date()} → {end.date()}")

    if outage_mask is not None:
        axes[0].axvspan(start, end, color="gray", alpha=0.2)

    plt.tight_layout()

def plot_price_per_bus(n, bus, start, end):
    start = pd.to_datetime(start)
    end = pd.to_datetime(end)

    prices = n.buses_t.marginal_price


    if bus not in prices.columns:
        print(f"Bus {bus} not found in network prices.")
        return

    series = prices[bus]
    mask = (series.index >= start) & (series.index <= end)

    fig, axes = plt.subplots(1, 2, figsize=(18, 6), sharey=True)

    # LEFT: full-year price
    axes[0].plot(series.index, series.values, label=bus)
    axes[0].set_title(f"Electricity Price at {bus} - Full Year")
    axes[0].set_ylabel("€/MWh")
    axes[0].grid(alpha=0.3)

    # Shade outage window
    axes[0].axvspan(start, end, color="gray", alpha=0.2)

    # RIGHT: outage zoom
    axes[1].plot(series.index[mask], series[mask], label=bus)
    axes[1].set_title(f"{bus} – {start.date()} → {end.date()}")
    axes[1].grid(alpha=0.3)

    plt.tight_layout()

############################### POWER FLOWS ################################
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

    outdir = f"results_contingencies/2025/bio_2025/{scenario}/plots"
    os.makedirs(outdir, exist_ok=True)

    start = pd.to_datetime(snakemake.params.start)
    duration = pd.Timedelta(days=snakemake.params.duration)
    end = start + duration
    mask = (n.snapshots >= start) & (n.snapshots <= end)

    scaling_factor = snakemake.params.scaling_factor
    buses = snakemake.params.buses

    print(f"Generating plots for scenario {scenario}…")

    # CP status
    plot_affected_biomass_generators(n, buses, start=start,end=end)
    plt.savefig(f"{outdir}/bio_status.png")
    plt.close()

    # Dispatch by carrier
    plot_dispatch_by_carrier(n, start, end)
    plt.savefig(f"{outdir}/dispatch_by_carrier.png")
    plt.close()

    # Dispatch by carrier at specific buses
    buses_to_plot = ["DK0 0", "DK1 0", "SE1 0"]
    for bus in buses_to_plot:
        plot_bus_dispatch_by_carrier(n, bus, start, end)
        plt.savefig(f"{outdir}/dispatch_by_carrier_{bus.replace(' ','_')}.png")
        plt.close()

    # Load shedding DK2
    plot_load_shedding_DK2(n, start, end, bus="DK1 0")
    plt.savefig(f"{outdir}/load_shedding_DK2.png")
    plt.close()

    # Load shedding per bus
    plot_load_shedding_per_bus(n, start, end)
    plt.savefig(f"{outdir}/load_shed_per_bus.png")
    plt.close()

    # Electricity Prices
    plot_prices_full_and_zoom(n, start, end, outage_mask=mask)
    plt.savefig(f"{outdir}/price_plot.png")
    plt.close()

    #Electricity prices per bus
    buses_to_plot = ["DK0 0", "DK1 0", "SE1 0", "DE0 0"]
    for bus in buses_to_plot:
        plot_price_per_bus(n, bus, start, end)
        plt.savefig(f"{outdir}/electricity_price_{bus.replace(' ','_')}.png")
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