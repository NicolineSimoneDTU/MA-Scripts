configfile: "contingency_stuff/c_wind_2035.yaml"

SCENARIOS = config["scenarios"]

rule all:
    input:
        # For each scenario
        expand("results_contingencies/2035/wind/{scenario}/finished.txt",
               scenario=SCENARIOS.keys()),
        expand("results_contingencies/2035/wind/{scenario}/plots_done.txt",
               scenario=SCENARIOS.keys()),
        
        # One global summary file
        "results_contingencies/2035/wind/summary_table.csv"

rule add_load_shedding:
    input:
        net=config["network_file"]
    output:
        "results_contingencies/2035/wind/{scenario}/network_ls.nc"
    params:
        elec_ls_cost=lambda w: config["load_shedding"]["elec_cost"],
        heat_ls_factor=lambda w: config["load_shedding"]["heat_factor"],
        add_heat_ls=lambda w: config["load_shedding"]["add_heat"],
        use_lv_buses=lambda w: config["load_shedding"]["use_lv_buses"],
        add_bio_ls=lambda w: config["load_shedding"].get("add_bio", False),
    script:
        "scripts/add_load_shedding.py"


rule apply_outage:
    input:
        net="results_contingencies/2035/wind/{scenario}/network_ls.nc"
    output:
        "results_contingencies/2035/wind/{scenario}/network.nc"
    params:
        start = lambda w: SCENARIOS[w.scenario]["start"],
        duration = lambda w: SCENARIOS[w.scenario]["duration_days"],
        buses = lambda w: SCENARIOS[w.scenario]["buses"],
        scaling_factor = lambda w: SCENARIOS[w.scenario]["scaling_factor"]
    script:
        "scripts/apply_wind_2035.py"

rule analyze:
    input:
        net="results_contingencies/2035/wind/{scenario}/network.nc"
    output:
        "results_contingencies/2035/wind/{scenario}/metrics.json"
    params:
        start=lambda w: SCENARIOS[w.scenario]["start"],
        duration=lambda w: SCENARIOS[w.scenario]["duration_days"]
    script:
        "scripts/analyze_metrics_2035.py"

rule plot_all:
    input:
        net="results_contingencies/2035/wind/{scenario}/network.nc"
    output:
        "results_contingencies/2035/wind/{scenario}/plots_done.txt"
    params:
        start=lambda w: SCENARIOS[w.scenario]["start"],
        duration=lambda w: SCENARIOS[w.scenario]["duration_days"],
        buses=lambda w: SCENARIOS[w.scenario]["buses"],
        scaling_factor=lambda w: SCENARIOS[w.scenario]["scaling_factor"]
    script:
        "scripts/plot_wind_2035.py"

rule summary:
    input:
        metrics=expand("results_contingencies/2035/wind/{scenario}/metrics.json",
                       scenario=SCENARIOS.keys())
    output:
        "results_contingencies/2035/wind/summary_table.csv"
    script:
        "scripts/build_summary_2035.py"

rule done:
    input:
        "results_contingencies/2035/wind/{scenario}/metrics.json"
    output:
        "results_contingencies/2035/wind/{scenario}/finished.txt"
    shell:
        "echo 'OK' > {output}"