#!/usr/bin/env python3
import ROOT
import argparse
from collections import defaultdict

def get_hierarchy_status(true_id, tree):
    """
    Looks up a True Geant4 TrackID in the reco arrays.
    Enforces TrackScore logic and uses track_Length to ignore delta-rays.
    Returns "Primary", "Secondary", or "Not Reco".
    """
    if true_id == -1:
        return "Not Applicable"
        
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
            # Only evaluate the hierarchy of the LONGEST track segment.
            # This prevents tiny primary delta-rays 
            if t_lengths[i] > max_length:
                max_length = t_lengths[i]
                best_is_primary = (t_is_prim[i] == 1)
                
    # Check Showers
    if not found:
        for i, tid in enumerate(s_true_ids):
            if tid == true_id and 0.0 <= s_scores[i] <= 0.5:
                found = True
                best_is_primary = (s_is_prim[i] == 1)
                if best_is_primary:
                    break # Assume the primary shower is the main one
                    
    if not found:
        return "Not Reco"
    elif best_is_primary:
        return "Primary"
    else:
        return "Secondary"

def run_hierarchy_diagnostics(input_filename, output_tag="hierarchy"):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    
    f = ROOT.TFile.Open(input_filename, "READ")
    tree = f.Get("ana/tree") or f.Get("tree")
    if not tree:
        print(f"Error: Could not find tree in {input_filename}")
        return

    particles = ["Muon (nu)", "Kaon (nu)", "Proton (Lambda)", "Pion (Lambda)", "Muon (Kaon Decay)"]
    
    # Track counts: { particle_name : { "Primary": 0, "Secondary": 0, "Not Reco": 0, "Total": 0 } }
    results = { p: {"Primary": 0, "Secondary": 0, "Not Reco": 0, "Total": 0} for p in particles }

    print("Scanning Geant4 Truth to Pandora Hierarchy mapping...")

    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)

        # True signal in FV
        if getattr(tree, "IsAssocLambdaKPlusWithProtonPiMinus", 0) != 1: continue
        if getattr(tree, "IsInsideFV", 0) != 1: continue

        true_mu_id, true_k_id, lambda_id = -1, -1, -1
        true_p_id, true_pi_id, true_k_mu_id = -1, -1, -1
        
        geant_pdgs = getattr(tree, "geant_PDG", [])
        geant_tids = getattr(tree, "geant_TrackID", [])
        geant_mids = getattr(tree, "geant_MotherID", [])

        # Neutrino daughters fnid
        for i, pdg in enumerate(geant_pdgs):
            mid = geant_mids[i]
            tid = geant_tids[i]
            if mid in [0, 10000000]:
                if pdg == 13: true_mu_id = tid
                if pdg == 321: true_k_id = tid
                if pdg == 3122: lambda_id = tid

        # Lambda daughters find
        for i, pdg in enumerate(geant_pdgs):
            mid = geant_mids[i]
            tid = geant_tids[i]
            if lambda_id != -1 and mid == lambda_id:
                if pdg == 2212: true_p_id = tid
                if pdg == -211: true_pi_id = tid
            if true_k_id != -1 and mid == true_k_id:
                # The Kaon decays to a muon ~63% of the time
                if pdg == -13: true_k_mu_id = tid

        # Map to TrackID
        id_map = {
            "Muon (nu)": true_mu_id,
            "Kaon (nu)": true_k_id,
            "Proton (Lambda)": true_p_id,
            "Pion (Lambda)": true_pi_id,
            "Muon (Kaon Decay)": true_k_mu_id
        }

        for p_name, t_id in id_map.items():
            if t_id != -1: # Only evaluate if the particle physically existed in the event
                status = get_hierarchy_status(t_id, tree)
                results[p_name][status] += 1
                results[p_name]["Total"] += 1

    f.Close()

    # Draw stacked histogram
    c1 = ROOT.TCanvas("c1", "PFP Hierarchy Breakdown", 1000, 600)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetRightMargin(0.2)
    
    n_bins = len(particles)
    
    h_prim = ROOT.TH1F("h_prim", "Pandora PFP Hierarchy Mapping;;Percentage of Events", n_bins, 0, n_bins)
    h_sec  = ROOT.TH1F("h_sec", "Secondary", n_bins, 0, n_bins)
    h_miss = ROOT.TH1F("h_miss", "Not Reco", n_bins, 0, n_bins)

    h_prim.SetFillColor(ROOT.kAzure+2)
    h_sec.SetFillColor(ROOT.kOrange+1)
    h_miss.SetFillColor(ROOT.kGray+1)

    for i, p_name in enumerate(particles):
        tot = results[p_name]["Total"]
        if tot > 0:
            pct_prim = (results[p_name]["Primary"] / tot) * 100.0
            pct_sec  = (results[p_name]["Secondary"] / tot) * 100.0
            pct_miss = (results[p_name]["Not Reco"] / tot) * 100.0
            
            h_prim.SetBinContent(i+1, pct_prim)
            h_sec.SetBinContent(i+1, pct_sec)
            h_miss.SetBinContent(i+1, pct_miss)
            
        h_prim.GetXaxis().SetBinLabel(i+1, p_name)

    hs = ROOT.THStack("hs", "Pandora PFP Hierarchy Mapping;;Percentage (%)")
    # To stop silly inversion rubbish
    hs.Add(h_miss)
    hs.Add(h_sec)
    hs.Add(h_prim)


    hs.Draw("BAR")
    
    hs.SetMaximum(100.0)
    hs.GetXaxis().SetLabelSize(0.05)
    
    leg = ROOT.TLegend(0.81, 0.4, 0.99, 0.8)
    leg.SetBorderSize(0)
    leg.AddEntry(h_prim, "Primary PFP", "f")
    leg.AddEntry(h_sec, "Secondary PFP", "f")
    leg.AddEntry(h_miss, "Not Reco'd", "f")
    leg.Draw()

    c1.SaveAs(f"{output_tag}_pfp_hierarchy.pdf")

    print("\n" + "="*85)
    print(" PANDORA PFP HIERARCHY BREAKDOWN (TRUE SIGNAL EVENTS IN FV)")
    print("="*85)
    print(f" {'Particle':<20} | {'Total Events':<12} | {'Primary %':<12} | {'Secondary %':<12} | {'Not Reco %'}")
    print("-" * 85)
    
    for p_name in particles:
        tot = results[p_name]["Total"]
        if tot == 0: continue
        
        pct_prim = (results[p_name]["Primary"] / tot) * 100.0
        pct_sec  = (results[p_name]["Secondary"] / tot) * 100.0
        pct_miss = (results[p_name]["Not Reco"] / tot) * 100.0
        
        print(f" {p_name:<20} | {tot:<12} | {pct_prim:>8.1f} %   | {pct_sec:>8.1f} %   | {pct_miss:>8.1f} %")
        
    print("="*85 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--signal", required=True, help="Filtered Hyperon ROOT file")
    parser.add_argument("-o", "--output_tag", default="hierarchy")
    args = parser.parse_args()
    
    run_hierarchy_diagnostics(args.signal, args.output_tag)