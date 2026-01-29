# apply_bio_2035.py

from snakemake.script import snakemake
import pypsa
import pandas as pd

# Load inputs
network_path = snakemake.input[0]
out_path = snakemake.output[0]

areas = snakemake.params.areas
scaling_factor = snakemake.params.scaling_factor

start = pd.to_datetime(snakemake.params.start)
duration = pd.Timedelta(days=snakemake.params.duration)

n = pypsa.Network(network_path)

# Create mask
if snakemake.params.duration == 0:
    mask = pd.Series(False, index=n.snapshots)
else:
    end = start + duration
    mask = (n.snapshots >= start) & (n.snapshots <= end)

# Removing EU bus
if "EU" in n.buses.index:
    print("Removing empty bus: EU")
    n.remove("Bus", "EU")

# Apply disruption
from my_functions import apply_biomass_shortage_2035, solve_dispatch

apply_biomass_shortage_2035(n, areas, scaling_factor, mask)
solve_dispatch(n)

# Save modified network
n.export_to_netcdf(out_path)