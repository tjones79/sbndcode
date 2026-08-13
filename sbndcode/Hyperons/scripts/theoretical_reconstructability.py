#!/usr/bin/env python3
import ROOT
import argparse
import math

def get_hierarchy_status(true_id, tree):
    """Checks if Pandora successfully built an object for this truth ID."""
    if true_id == -1: return "Not Applicable"
        
    t_true_ids = getattr(tree, "track_TrueGeantID", [])
    t_is_prim = getattr(tree, "track_IsPrimary", [])
    t_scores = getattr(tree, "track_TrackScore", [])
    t_lengths = getattr(tree, "track_Length", [])
    
    s_true_ids = getattr(tree, "shower_TrueGeantID", [])
    s_is_prim = getattr(tree, "shower_IsPrimary", [])
    s_scores = getattr(tree, "shower_TrackScore", [])
    
    found = False
    best_is_primary = False
    max_length = -1.0
    
    for i, tid in enumerate(t_true_ids):
        if tid == true_id and t_scores[i] > 0.5:
            found = True
            if t_lengths[i] > max_length:
                max_length = t_lengths[i]
                best_is_primary = (t_is_prim[i] == 1)
                
    if not found:
        for i, tid in enumerate(s_true_ids):
            if tid == true_id and 0.0 <= s_scores[i] <= 0.5:
                found = True
                best_is_primary = (s_is_prim[i] == 1)
                if best_is_primary: break 
                    
    if not found: return "Not Reco"
    elif best_is_primary: return "Primary"
    else: return "Secondary"

def run_exact_reconstructability(input_filename, output_tag="exact_reconstructability"):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0) 
    
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("ana/tree") or input_file.Get("tree")
    if not tree:
        print("Error: Could not find tree in file.")
        return

    hierarchy_particles = ["Muon (nu)", "Kaon (nu)", "Proton (Lambda)", "Pion (Lambda)"]
    
    mom_ranges = {
        "Muon (nu)": (40, 0.0, 2.0),
        "Kaon (nu)": (40, 0.0, 1.0),
        "Proton (Lambda)": (40, 0.0, 1.5),
        "Pion (Lambda)": (40, 0.0, 0.5)
    }
    
    h_total = {}
    h_theory_pass = {}
    h_pandora_pass = {}

    for p in hierarchy_particles:
        n_bins, x_min, x_max = mom_ranges[p]
        h_total[p] = ROOT.TH1F(f"h_tot_{p}", f"{p};True Momentum [GeV/c];Events", n_bins, x_min, x_max)
        h_theory_pass[p] = ROOT.TH1F(f"h_theory_{p}", f"{p};True Momentum [GeV/c];Events", n_bins, x_min, x_max)
        h_pandora_pass[p] = ROOT.TH1F(f"h_pandora_{p}", f"{p};True Momentum [GeV/c];Events", n_bins, x_min, x_max)

    num_entries = tree.GetEntries()
    for entry in range(num_entries):
        tree.GetEntry(entry)
        
        # Filter: ONLY look at true signal events in the FV
        if getattr(tree, "IsAssocLambdaKPlusWithProtonPiMinus", 0) != 1 or getattr(tree, "IsInsideFV", 0) != 1: 
            continue

        # Extract Geant4 indices, IDs, and Momenta
        particles_data = {
            "Muon (nu)": {"idx": -1, "id": -1, "mom": -1.0},
            "Kaon (nu)": {"idx": -1, "id": -1, "mom": -1.0},
            "Proton (Lambda)": {"idx": -1, "id": -1, "mom": -1.0},
            "Pion (Lambda)": {"idx": -1, "id": -1, "mom": -1.0}
        }
        
        lambda_id = -1
        
        for i, pdg in enumerate(tree.geant_PDG):
            mid = tree.geant_MotherID[i]
            tid = tree.geant_TrackID[i]
            p_mag = math.sqrt(tree.geant_Px[i]**2 + tree.geant_Py[i]**2 + tree.geant_Pz[i]**2)
            
            if pdg == 13 and mid in [0, 10000000]: 
                particles_data["Muon (nu)"] = {"idx": i, "id": tid, "mom": p_mag}
            elif pdg == 321 and mid in [0, 10000000]: 
                particles_data["Kaon (nu)"] = {"idx": i, "id": tid, "mom": p_mag}
            elif pdg == 3122 and mid in [0, 10000000]: 
                lambda_id = tid
        # find lambda daughters, could make this a helper function for all scripts really
        if lambda_id != -1:
            for i, pdg in enumerate(tree.geant_PDG):
                if tree.geant_MotherID[i] == lambda_id:
                    p_mag = math.sqrt(tree.geant_Px[i]**2 + tree.geant_Py[i]**2 + tree.geant_Pz[i]**2)
                    
                    if pdg == 2212: 
                        particles_data["Proton (Lambda)"] = {"idx": i, "id": tree.geant_TrackID[i], "mom": p_mag}
                    elif pdg == -211: 
                        particles_data["Pion (Lambda)"] = {"idx": i, "id": tree.geant_TrackID[i], "mom": p_mag}

        # Apply the exact Pandora View hit thresholds
        for p_name, data in particles_data.items():
            if data["idx"] != -1 and data["mom"] >= 0:
                h_total[p_name].Fill(data["mom"])
                
                # Check exact hits
                u = tree.geant_HitsU[data["idx"]]
                v = tree.geant_HitsV[data["idx"]]
                w = tree.geant_HitsW[data["idx"]]
                total_hits = u + v + w

                # pandora criteria
                is_reconstructable = (total_hits > 15) and (
                    (u > 5 and v > 5) or 
                    (u > 5 and w > 5) or 
                    (v > 5 and w > 5)
                )

                if is_reconstructable:
                    h_theory_pass[p_name].Fill(data["mom"])
                
                # Did pandora successfully reconstruct this particle?
                status = get_hierarchy_status(data["id"], tree)
                if status in ["Primary", "Secondary"]:
                    h_pandora_pass[p_name].Fill(data["mom"])

    c3 = ROOT.TCanvas("c3", "Exact Theoretical vs Actual Efficiency", 1200, 1000)
    c3.Divide(2, 2)
    
    efficiencies = [] 
    draw_objects = [] 
    
    for i, p_name in enumerate(hierarchy_particles):
        c3.cd(i+1)
        ROOT.gPad.SetGridy(1)
        ROOT.gPad.SetGridx(1)
        
        # Theoretical Efficiency (Green) vs Actual Pandora Efficiency (Blue)
        if ROOT.TEfficiency.CheckConsistency(h_theory_pass[p_name], h_total[p_name]):
            eff_theory = ROOT.TEfficiency(h_theory_pass[p_name], h_total[p_name])
            eff_theory.SetTitle(f"{p_name} Efficiency (Is Reconstructable? vs Actual Pandora);True Momentum [GeV/c];Efficiency")
            eff_theory.SetMarkerStyle(21)
            eff_theory.SetMarkerColor(ROOT.kGreen+2)
            eff_theory.SetLineColor(ROOT.kGreen+2)
            eff_theory.Draw("AP")
            efficiencies.append(eff_theory)
            
            ROOT.gPad.Update()
            graph = eff_theory.GetPaintedGraph()
            if graph: graph.GetYaxis().SetRangeUser(0.0, 1.05)

        if ROOT.TEfficiency.CheckConsistency(h_pandora_pass[p_name], h_total[p_name]):
            eff_pandora = ROOT.TEfficiency(h_pandora_pass[p_name], h_total[p_name])
            eff_pandora.SetMarkerStyle(20)
            eff_pandora.SetMarkerColor(ROOT.kAzure+2)
            eff_pandora.SetLineColor(ROOT.kAzure+2)
            eff_pandora.Draw("P SAME")
            efficiencies.append(eff_pandora)
            
        leg = ROOT.TLegend(0.40, 0.25, 0.88, 0.45)
        leg.SetBorderSize(0)
        leg.AddEntry(eff_theory, "Is Reconstructable", "pe")
        leg.AddEntry(eff_pandora, "Actual (Pandora PFP Found)", "pe")
        leg.Draw()
        draw_objects.append(leg)

    c3.SaveAs(f"{output_tag}_overlay.pdf")
    print(f"\nOverlay plots saved to {output_tag}_overlay.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    args = parser.parse_args()
    run_exact_reconstructability(args.input)