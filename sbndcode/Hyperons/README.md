# SBND Hyperon Analysis Pipeline

This repository contains the LArSoft modules, FHiCL configurations, and Snakemake workflow for the $\Lambda^0$ hyperon analysis on SBND.

The primary focus of this analysis is identifying $\Lambda^0$ hyperons generated via **associated production**. In Cabibbo-favored neutrino interactions ($\Delta S = 0$), strangeness is conserved in the strong/electromagnetic vertex, meaning the $\Lambda^0$ (strangeness $S=-1$) must be produced alongside another strange particle, typically a Kaon (e.g., $K^+$, with $S=+1$). Therefore, our complete signal signature requires finding a detached $\Lambda^0$ decay topology in coincidence with an associated Kaon track originating directly from the primary neutrino vertex.

The pipeline is designed to process `reco2` LArSoft files, flatten the complex data into standard ROOT `TTrees`, and apply topological and kinematic selections to isolate hyperon decays. This is a work in progress, so there are many apects of the selection to expand on. I will attempt to keep it up to date :)

## 1. Target Decay Chains & Topologies

The primary signal channel driven by the Charged Current (CC) neutrino interaction follows this full production and decay chain:

### A. Primary Associated Production Interaction
$$\nu_{\mu} + N \rightarrow \mu^{-} + \Lambda^0 + K^{+} + N'$$

*   **$\mu^{-}$ (Prompt Track):** A long, minimum-ionizing particle (MIP) track originating from the primary neutrino vertex.
*   **$K^{+}$ (Prompt Track):** A heavily ionizing track originating from the primary neutrino vertex, acting as the explicit tag for strangeness conservation ($\Delta S = 0$).
*   **$\Lambda^0$ (Invisible Gap):** A neutral hyperon that leaves no ionization trace in the Liquid Argon TPC as it travels away from the primary vertex.

### B. Secondary Hyperon Decay (The V0 Channel)
The neutral $\Lambda^0$ travels a short distance before decaying via the weak interaction. We target the charged decay mode (Branching Ratio $\approx 63.9\%$):
$$\Lambda^0 \rightarrow p + \pi^{-}$$

*   **$p$ (Proton Track):** A highly ionizing track originating at a secondary vertex (V0), typically stopping in the active volume and exhibiting a distinct Bragg peak.
*   **$\pi^{-}$ (Pion Track):** A lower-ionizing, MIP-like track originating from the same secondary V0 vertex.

## 2. The Analyser Stage (`AnalyzeHyperon_module.cc`)
The core of the analysis is the `AnalyzeHyperon` LArSoft analyzer module.

**Key operations performed in the Analyser:**
* **Data Unpacking:** Reads in reconstructed objects from the `reco2` stage, specifically `recob::Track`, `recob::Vertex`, `recob::Shower`, and PID objects.
* **Truth-Matching:** Utilises LArSoft BackTracker services to match reconstructed tracks to their true simulated `simb::MCParticle` counterparts (crucial for evaluating signal efficiency and background misidentification).
* **Tree Flattening:** Extracts event-level variables, track kinematics (start/end points, direction, momentum), and calorimetry metrics (dE/dx), and writes them into a flat `TTree` using the `TFileService`.

## 3. Event Selection Strategy
The physics selection is designed to identify the classic detached "V0" topology characteristic of hyperon decays, specifically $\Lambda^0 \rightarrow p + \pi^-$. 

The selection is applied through the following hierarchical cuts:

### A. Pre-Selection & Quality
* **Beam Timing:** Events must fall within the BNB hardware trigger window.
* **Fiducial Volume:** The reconstructed primary vertex (PV) and the detached V0 vertex must be fully contained within the defined SBND active TPC volume.

### B. Slice Selection & Cosmic Rejection (Pandora)
Before looking for specific hyperon topologies, the analyser must isolate the actual neutrino interaction from the surrounding cosmic ray background. During the `reco2` stage, the Pandora pattern recognition algorithms group clusters of hits into distinct, isolated interactions known as **"slices."**

To select the correct slice, the analyser performs the following:
* **Topological Scoring:** Evaluates the `Pandora NuScore` (a multivariate score distinguishing neutrino-like topologies from cosmic-like topologies) for each `recob::Slice`. 
* **Neutrino Slice Identification:** The slice with the highest NuScore is explicitly tagged as the primary neutrino interaction.
* **PFParticle Isolation:** All subsequent analysis (V0 finding, PID, and track kinematics) is strictly restricted to the `recob::PFParticle` objects (tracks and showers) contained entirely within this selected neutrino slice, ignoring all other reconstructed objects in the event. i.e. Remove reconstructed objects that are not contained within the neutrino tagged slice.

### C. Topological Selection (The V0 Finder) **TODO**
Because the $\Lambda^0$ is neutral, it leaves no track in the LAr. Can select events based on its decay products:
* **Two-Track Signature:** Look for a secondary vertex (V0) with exactly two outgoing reconstructed tracks.
* **Detachment Cut:** The distance between the neutrino Primary Vertex and the V0 vertex must be greater than `[X] cm` to reject tracks originating directly from the primary interaction.
* **Collinearity/Pointing:** The momentum vector sum of the two V0 tracks must point back to the Primary Vertex.


## 4. The Workflow (Snakemake & Grid)
Processing thousands of ROOT files is handled via a containerised Snakemake workflow.

* **Environment:** The pipeline uses Apptainer to mount the official Fermilab `fnal-dev-sl7` image, ensuring complete OS compatibility.
* **Setup Scripts:** `snakemake_setup.sh` is used to give paths for inputs and local larsoft developer areas. The `setup_env.sh` automatically configures the base `sbndcode v10_14_02_03` release, sources the local compiled products, and wires the `$FHICL_FILE_PATH` to the working directory.
* **Execution:** Snakemake runs `lar -c run_analyze_hyperon_events.fcl` over the input datasets (both signal and background). It will also perform various other rules that are yet to be expanded on. 

***

### How to Run
To process the latest batch of files locally or interactively:
```bash
# Clean previous shadows
snakemake --cleanup-shadow

# Execute pipeline (adjust cores as needed)
snakemake --cores 1

# I also like to use SLURM. Something like the following can be eaily used to run snakemake on slurm.
snakemake --cluster "sbatch --time={resources.time_min} --mem={resources.mem_mb} --cpus-per-task={resources.cpus}" --jobs 1000 --latency-wait 60
```