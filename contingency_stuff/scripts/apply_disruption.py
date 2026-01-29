from snakemake.script import snakemake
import pypsa
import pandas as pd

# Load inputs
network_path = snakemake.input[0]
out_path = snakemake.output[0]

start = pd.to_datetime(snakemake.params.start)
duration = pd.Timedelta(days=snakemake.params.duration)

lines = snakemake.params.lines
links = snakemake.params.links

n = pypsa.Network(network_path)

# Create mask
if snakemake.params.duration == 0:
    mask = pd.Series(False, index=n.snapshots)
else:
    end = start + duration
    mask = (n.snapshots >= start) & (n.snapshots <= end)

# Apply disruptions
from my_functions import apply_line_link_outages, solve_dispatch

apply_line_link_outages(n, lines, links, mask)
solve_dispatch(n)

# Save modified network
n.export_to_netcdf(out_path)