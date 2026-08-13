#!/usr/bin/env python3
import ROOT
import argparse
import math
from collections import defaultdict

# Slicing
def get_slices(true_id, tree):
    slices = set()
    if true_id == -1: return slices
    
    t_true_ids = getattr(tree, "track_TrueGeantID", [])
    t_slice_ids = getattr(tree, "track_SliceID", [])
    t_scores = getattr(tree, "track_TrackScore", [])
    for i, tid in enumerate(t_true_ids):
        if tid == true_id and t_slice_ids[i] != -1 and t_scores[i] > 0.5: 
            slices.add(t_slice_ids[i])
    
    s_true_ids = getattr(tree, "shower_TrueGeantID", [])
    s_slice_ids = getattr(tree, "shower_SliceID", [])
    s_scores = getattr(tree, "shower_TrackScore", [])
    for i, tid in enumerate(s_true_ids):
        if tid == true_id and s_slice_ids[i] != -1 and 0.0 <= s_scores[i] <= 0.5: 
            slices.add(s_slice_ids[i])
        
    return slices

# Tracking Counts (Tracks vs Showers)
def get_reco_counts(true_id, tree):
    if true_id == -1: return 0, 0
    n_t, n_s = 0, 0
    
    t_true_ids = getattr(tree, "track_TrueGeantID", [])
    t_scores = getattr(tree, "track_TrackScore", [])
    for i, tid in enumerate(t_true_ids):
        if tid == true_id and t_scores[i] > 0.5: n_t += 1

    s_true_ids = getattr(tree, "shower_TrueGeantID", [])
    s_scores = getattr(tree, "shower_TrackScore", [])
    for i, tid in enumerate(s_true_ids):
        if tid == true_id and 0.0 <= s_scores[i] <= 0.5: n_s += 1
            
    return n_t, n_s

# Hierarchy Status (Primary, Secondary, Not Reco)
def get_hierarchy_status(true_id, tree):
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


def run_diagnostics(input_filename, output_tag="reco"):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0) 
    
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("ana/tree")
    if not tree:
        print("Error: Could not find tree in file.")
        return

    h_slices = ROOT.TH1F("h_slices", "Signal Slices per Event;Number of Slices (TrueNuHits > 5);Events", 6, -0.5, 5.5)
    h_slices.SetFillColorAlpha(ROOT.kAzure+1, 0.7)

    h2_slices_vs_dist = ROOT.TH2F("h2_slices_vs_dist", "Fracturing vs. #Lambda Flight Distance;True #Lambda Flight Distance [cm];Number of Signal Slices", 20, 0, 20, 5, 0.5, 5.5)
    p_slices_vs_dist = ROOT.TProfile("p_slices_vs_dist", "Average Fracturing vs. #Lambda Flight Distance;True #Lambda Flight Distance [cm];Average Signal Slices", 20, 0, 20)
    p_slices_vs_dist.SetLineColor(ROOT.kRed); p_slices_vs_dist.SetLineWidth(2); p_slices_vs_dist.SetMarkerStyle(20); p_slices_vs_dist.SetMarkerColor(ROOT.kRed)

    h_fractures = ROOT.TH1F("h_fractures", "How the Event Fractures (Excludes Perfectly Contained);;Events", 8, 0, 8)

    categories = ["Perfectly Contained", "Muon Track Chopped", "Proton/Pion Chopped", "Intact V-Vertex Split", "Shattered V-Vertex", "Only Proton Split", "Only Pion Split", "Kaon Split", "Junk / Unreconstructed"]
    colors = { "Perfectly Contained": ROOT.kTeal+1, "Muon Track Chopped": ROOT.kRed-4, "Proton/Pion Chopped": ROOT.kOrange+1, "Intact V-Vertex Split": ROOT.kMagenta+1, "Shattered V-Vertex": ROOT.kMagenta+3, "Only Proton Split": ROOT.kBlue-4, "Only Pion Split": ROOT.kCyan+1, "Kaon Split": ROOT.kGreen+2, "Junk / Unreconstructed": ROOT.kGray+1 }
    h_comp = {cat: ROOT.TH1F(f"h_comp_{cat.replace(' ', '_').replace('/', '_')}", "Primary Slice Completeness", 40, 0.0, 1.05) for cat in categories}
    for cat in categories: h_comp[cat].SetFillColor(colors[cat]); h_comp[cat].SetLineColor(ROOT.kBlack); h_comp[cat].SetLineWidth(1)

    daughter_cats = ["Two Distinct Tracks", "Pion Lost/Merged", "Proton Lost/Merged", "V-Vertex Showered", "Neither Reconstructed"]
    h_daughters = ROOT.TH1F("h_daughters", "Lambda Daughter Tracking;;Events", 5, 0, 5)
    h_daughters.SetFillColorAlpha(ROOT.kSpring-5, 0.8)

    hierarchy_particles = ["Muon (nu)", "Kaon (nu)", "Proton (Lambda)", "Pion (Lambda)"]
    hier_statuses = ["Primary", "Secondary", "Not Reco"]
    hier_colors = {"Primary": ROOT.kAzure+2, "Secondary": ROOT.kOrange+1, "Not Reco": ROOT.kGray+1}
    
    mom_ranges = {
        "Muon (nu)": (40, 0.0, 2.0),
        "Kaon (nu)": (40, 0.0, 1.0),
        "Proton (Lambda)": (40, 0.0, 1.5),
        "Pion (Lambda)": (40, 0.0, 0.5)
    }
    
    h_hier_mom = {p: {} for p in hierarchy_particles}
    for p in hierarchy_particles:
        n_bins, x_min, x_max = mom_ranges[p]
        for status in hier_statuses:
            h_name = f"h_hier_mom_{p.replace(' ', '_').replace('(', '').replace(')', '')}_{status.replace(' ', '_')}"
            h_hier_mom[p][status] = ROOT.TH1F(h_name, f"{p} Hierarchy vs True Momentum;True Momentum [GeV/c];Events", n_bins, x_min, x_max)
            h_hier_mom[p][status].SetFillColor(hier_colors[status])
            h_hier_mom[p][status].SetLineColor(ROOT.kBlack)
            h_hier_mom[p][status].SetLineWidth(1)

    total_events, contained_slices, fractured_slices, pileup_events = 0, 0, 0, 0
    fracture_modes = {k: 0 for k in categories if k != "Perfectly Contained"}
    daughter_modes = {k: 0 for k in daughter_cats}
    
    hierarchy_results = { p: {"Primary": 0, "Secondary": 0, "Not Reco": 0, "Total": 0} for p in hierarchy_particles }

    num_entries = tree.GetEntries()
    for entry in range(num_entries):
        tree.GetEntry(entry)
        if entry > 0 and entry % 5000 == 0: print(f"{entry}/{num_entries} events processed")
        
        # Filter: ONLY look at true signal events in the FV
        if getattr(tree, "IsAssocLambdaKPlusWithProtonPiMinus", 0) != 1 or getattr(tree, "IsInsideFV", 0) != 1: continue
            
        total_events += 1
        if len(tree.nu_pdg) > 1: pileup_events += 1
        
        # Slicing
        signal_slice_count, max_slice_hits = 0, 0
        for hits in tree.slice_TrueNuHits:
            if hits > 5: signal_slice_count += 1
            if hits > max_slice_hits: max_slice_hits = hits
        h_slices.Fill(signal_slice_count)
        
        completeness = -1.0
        event_total_hits = getattr(tree, "event_TotalTrueNuHits", 0)
        if event_total_hits > 0: completeness = float(max_slice_hits) / float(event_total_hits)

        lambda_flight_dist = -1.0
        if hasattr(tree, "geant_PDG"):
            for l_idx, pdg in enumerate(tree.geant_PDG):
                if pdg == 3122:
                    l_dx = tree.geant_endX[l_idx] - tree.geant_startX[l_idx]
                    l_dy = tree.geant_endY[l_idx] - tree.geant_startY[l_idx]
                    l_dz = tree.geant_endZ[l_idx] - tree.geant_startZ[l_idx]
                    lambda_flight_dist = math.sqrt(l_dx**2 + l_dy**2 + l_dz**2)
                    break 
        
        if lambda_flight_dist >= 0:
            h2_slices_vs_dist.Fill(lambda_flight_dist, signal_slice_count)
            p_slices_vs_dist.Fill(lambda_flight_dist, signal_slice_count)

        # G4 momenta and IDs
        true_mu_id, true_k_id, lambda_id, true_p_id, true_pi_id = -1, -1, -1, -1, -1
        momenta = { "Muon (nu)": -1.0, "Kaon (nu)": -1.0, "Proton (Lambda)": -1.0, "Pion (Lambda)": -1.0 }
        
        for i, pdg in enumerate(tree.geant_PDG):
            mid = tree.geant_MotherID[i]
            tid = tree.geant_TrackID[i]
            p_mag = math.sqrt(tree.geant_Px[i]**2 + tree.geant_Py[i]**2 + tree.geant_Pz[i]**2)
            
            if pdg == 13 and mid in [0, 10000000]: 
                true_mu_id = tid
                momenta["Muon (nu)"] = p_mag
            if pdg == 321 and mid in [0, 10000000]: 
                true_k_id = tid
                momenta["Kaon (nu)"] = p_mag
            if pdg == 3122 and mid in [0, 10000000]: 
                lambda_id = tid

        if lambda_id != -1:
            for i, pdg in enumerate(tree.geant_PDG):
                if tree.geant_MotherID[i] == lambda_id:
                    p_mag = math.sqrt(tree.geant_Px[i]**2 + tree.geant_Py[i]**2 + tree.geant_Pz[i]**2)
                    if pdg == 2212: 
                        true_p_id = tree.geant_TrackID[i]
                        momenta["Proton (Lambda)"] = p_mag
                    if pdg == -211: 
                        true_pi_id = tree.geant_TrackID[i]
                        momenta["Pion (Lambda)"] = p_mag

        # Check is primary, secondary, or not reconstructed for each particle
        id_map = { "Muon (nu)": true_mu_id, "Kaon (nu)": true_k_id, "Proton (Lambda)": true_p_id, "Pion (Lambda)": true_pi_id }
        for p_name, t_id in id_map.items():
            if t_id != -1:
                status = get_hierarchy_status(t_id, tree)
                hierarchy_results[p_name][status] += 1
                hierarchy_results[p_name]["Total"] += 1
                
                # Fill the momentum histogram
                if momenta[p_name] >= 0:
                    h_hier_mom[p_name][status].Fill(momenta[p_name])

        # Get mah daughters
        n_p_tracks, n_p_showers = get_reco_counts(true_p_id, tree)
        n_pi_tracks, n_pi_showers = get_reco_counts(true_pi_id, tree)
        
        if n_p_showers > 0 or n_pi_showers > 0: daughter_modes["V-Vertex Showered"] += 1
        elif n_p_tracks > 0 and n_pi_tracks > 0: daughter_modes["Two Distinct Tracks"] += 1
        elif n_p_tracks > 0 and n_pi_tracks == 0: daughter_modes["Pion Lost/Merged"] += 1
        elif n_pi_tracks > 0 and n_p_tracks == 0: daughter_modes["Proton Lost/Merged"] += 1
        else: daughter_modes["Neither Reconstructed"] += 1

        # Are we fractured?
        if signal_slice_count == 1: 
            contained_slices += 1
            if completeness >= 0: h_comp["Perfectly Contained"].Fill(completeness)
        elif signal_slice_count > 1: 
            fractured_slices += 1
            slices_mu = get_slices(true_mu_id, tree); slices_k = get_slices(true_k_id, tree)
            slices_p = get_slices(true_p_id, tree); slices_pi = get_slices(true_pi_id, tree)

            is_mu_chopped = len(slices_mu) > 1
            is_p_chopped = len(slices_p) > 1
            is_pi_chopped = len(slices_pi) > 1
            
            p_sep = len(slices_p) > 0 and slices_p.isdisjoint(slices_mu)
            pi_sep = len(slices_pi) > 0 and slices_pi.isdisjoint(slices_mu)
            k_sep = len(slices_k) > 0 and slices_k.isdisjoint(slices_mu)

            assigned_cat = "Junk / Unreconstructed"
            if is_mu_chopped: assigned_cat = "Muon Track Chopped"
            elif is_p_chopped or is_pi_chopped: assigned_cat = "Proton/Pion Chopped"
            elif p_sep and pi_sep and not slices_p.isdisjoint(slices_pi): assigned_cat = "Intact V-Vertex Split"
            elif p_sep and pi_sep and slices_p.isdisjoint(slices_pi): assigned_cat = "Shattered V-Vertex"
            elif p_sep and not pi_sep: assigned_cat = "Only Proton Split"
            elif pi_sep and not p_sep: assigned_cat = "Only Pion Split"
            elif k_sep and not p_sep and not pi_sep: assigned_cat = "Kaon Split"

            fracture_modes[assigned_cat] += 1
            if completeness >= 0: h_comp[assigned_cat].Fill(completeness)


    c1 = ROOT.TCanvas("c1", "Diagnostics", 800, 600)
    
    c1.Clear(); h_slices.Draw("HIST"); c1.SaveAs(f"{output_tag}_signal_slices.pdf")
    
    c1.Clear(); ROOT.gPad.SetRightMargin(0.15); h2_slices_vs_dist.Draw("COLZ"); c1.SaveAs(f"{output_tag}_slices_vs_dist_2D.pdf"); ROOT.gPad.SetRightMargin(0.05)
    
    c1.Clear(); p_slices_vs_dist.SetMinimum(0); p_slices_vs_dist.Draw("PE"); c1.SaveAs(f"{output_tag}_slices_vs_dist_profile.pdf")

    c1.Clear(); ROOT.gPad.SetBottomMargin(0.2); ROOT.gPad.SetLeftMargin(0.15); h_fractures.SetFillColorAlpha(ROOT.kMagenta-4, 0.8)
    for i, (cat, count) in enumerate(fracture_modes.items()): h_fractures.GetXaxis().SetBinLabel(i+1, cat); h_fractures.SetBinContent(i+1, count)
    h_fractures.GetXaxis().SetLabelSize(0.045); h_fractures.SetBarWidth(0.8); h_fractures.SetBarOffset(0.1); h_fractures.Draw("BAR"); c1.SaveAs(f"{output_tag}_fracture_modes.pdf"); ROOT.gPad.SetBottomMargin(0.1)

    c1.Clear(); ROOT.gPad.SetLeftMargin(0.12); ROOT.gPad.SetLogy(1)
    hs_comp = ROOT.THStack("hs_comp", "Primary Slice Completeness by Fracture Mode;Completeness (Primary Slice NuHits / Total Event NuHits);Events")
    for cat in reversed(categories): hs_comp.Add(h_comp[cat])
    hs_comp.Draw("HIST")
    leg = ROOT.TLegend(0.15, 0.45, 0.45, 0.88); leg.SetBorderSize(0)
    for cat in categories: leg.AddEntry(h_comp[cat], cat, "f")
    leg.Draw(); c1.SaveAs(f"{output_tag}_slice_completeness.pdf"); ROOT.gPad.SetLogy(0)

    c1.Clear(); ROOT.gPad.SetBottomMargin(0.2); ROOT.gPad.SetLeftMargin(0.15)
    for i, (cat, count) in enumerate(daughter_modes.items()): h_daughters.GetXaxis().SetBinLabel(i+1, cat); h_daughters.SetBinContent(i+1, count)
    h_daughters.GetXaxis().SetLabelSize(0.045); h_daughters.SetBarWidth(0.8); h_daughters.SetBarOffset(0.1); h_daughters.Draw("BAR"); c1.SaveAs(f"{output_tag}_daughter_tracking.pdf")

    c1.Clear()
    ROOT.gPad.SetBottomMargin(0.15); ROOT.gPad.SetRightMargin(0.2)
    n_bins = len(hierarchy_particles)
    h_prim = ROOT.TH1F("h_prim", "Pandora PFP Hierarchy Mapping;;Percentage of Events", n_bins, 0, n_bins)
    h_sec  = ROOT.TH1F("h_sec", "Secondary", n_bins, 0, n_bins)
    h_miss = ROOT.TH1F("h_miss", "Not Reco", n_bins, 0, n_bins)
    h_prim.SetFillColor(ROOT.kAzure+2); h_sec.SetFillColor(ROOT.kOrange+1); h_miss.SetFillColor(ROOT.kGray+1)

    for i, p_name in enumerate(hierarchy_particles):
        tot = hierarchy_results[p_name]["Total"]
        if tot > 0:
            h_prim.SetBinContent(i+1, (hierarchy_results[p_name]["Primary"] / tot) * 100.0)
            h_sec.SetBinContent(i+1, (hierarchy_results[p_name]["Secondary"] / tot) * 100.0)
            h_miss.SetBinContent(i+1, (hierarchy_results[p_name]["Not Reco"] / tot) * 100.0)
        h_prim.GetXaxis().SetBinLabel(i+1, p_name)

    hs_hier = ROOT.THStack("hs_hier", "Pandora PFP Hierarchy Mapping;;Percentage (%)")
    hs_hier.Add(h_miss); hs_hier.Add(h_sec); hs_hier.Add(h_prim)
    hs_hier.Draw("BAR"); hs_hier.SetMaximum(100.0); hs_hier.GetXaxis().SetLabelSize(0.05)
    
    leg_hier = ROOT.TLegend(0.81, 0.4, 0.99, 0.8); leg_hier.SetBorderSize(0)
    leg_hier.AddEntry(h_prim, "Primary PFP", "f"); leg_hier.AddEntry(h_sec, "Secondary PFP", "f"); leg_hier.AddEntry(h_miss, "Not Reco'd", "f")
    leg_hier.Draw(); c1.SaveAs(f"{output_tag}_pfp_hierarchy.pdf")

    # Mom stacks
    c2 = ROOT.TCanvas("c2", "Hierarchy vs Momentum", 1200, 1000)
    c2.Divide(2, 2)
    
    leg_mom = ROOT.TLegend(0.45, 0.65, 0.88, 0.88)
    leg_mom.SetBorderSize(0)
    leg_mom.AddEntry(h_hier_mom["Muon (nu)"]["Primary"], "Primary PFP", "f")
    leg_mom.AddEntry(h_hier_mom["Muon (nu)"]["Secondary"], "Secondary PFP", "f")
    leg_mom.AddEntry(h_hier_mom["Muon (nu)"]["Not Reco"], "Not Reco'd", "f")

    momentum_stacks = []
    draw_objects = []

    # LArTPC Approximate Momentum Thresholds for a ~1.5cm track
    thresholds = {
        "Muon (nu)": 0.05,
        "Pion (Lambda)": 0.06,
        "Kaon (nu)": 0.13,
        "Proton (Lambda)": 0.20
    }

    for i, p_name in enumerate(hierarchy_particles):
        c2.cd(i+1)
        hs = ROOT.THStack(f"hs_mom_{i}", f"{p_name} Hierarchy vs True Momentum;True Momentum [GeV/c];Events")
        # Stack order: Not Reco (bottom), Secondary, Primary (top)
        hs.Add(h_hier_mom[p_name]["Not Reco"])
        hs.Add(h_hier_mom[p_name]["Secondary"])
        hs.Add(h_hier_mom[p_name]["Primary"])
        hs.Draw("HIST")
        momentum_stacks.append(hs)
        
        ROOT.gPad.Update()
        max_y = hs.GetMaximum()
        
        thresh = thresholds[p_name]
        line = ROOT.TLine(thresh, 0, thresh, max_y)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2) # Dashed line
        line.SetLineWidth(2)
        line.Draw()
        draw_objects.append(line)
        
        text = ROOT.TLatex(thresh + (hs.GetXaxis().GetXmax() * 0.02), max_y * 0.85, f"Threshold ~{thresh:.2f} GeV/c")
        text.SetTextSize(0.04)
        text.SetTextColor(ROOT.kRed)
        text.Draw()
        draw_objects.append(text)

        leg_mom.Draw()
        
    c2.SaveAs(f"{output_tag}_hierarchy_vs_momentum.pdf")


    c3 = ROOT.TCanvas("c3", "Tracking Efficiency vs Momentum", 1200, 1000)
    c3.Divide(2, 2)
    
    efficiencies = []
    
    for i, p_name in enumerate(hierarchy_particles):
        c3.cd(i+1)
        ROOT.gPad.SetGridy(1)
        ROOT.gPad.SetGridx(1)
        
        # 1. Total Particles = Primary + Secondary + Not Reco
        h_total = h_hier_mom[p_name]["Primary"].Clone(f"h_tot_{i}")
        h_total.Add(h_hier_mom[p_name]["Secondary"])
        h_total.Add(h_hier_mom[p_name]["Not Reco"])
        
        # 2. Reconstructed Particles = Primary + Secondary
        h_pass = h_hier_mom[p_name]["Primary"].Clone(f"h_pass_{i}")
        h_pass.Add(h_hier_mom[p_name]["Secondary"])
        
        if ROOT.TEfficiency.CheckConsistency(h_pass, h_total):
            pEff = ROOT.TEfficiency(h_pass, h_total)
            pEff.SetTitle(f"{p_name} Reco Efficiency vs True Momentum;True Momentum [GeV/c];Reconstruction Efficiency")
            pEff.SetMarkerStyle(20)
            pEff.SetMarkerColor(ROOT.kAzure+2)
            pEff.SetLineColor(ROOT.kAzure+2)
            pEff.Draw("AP")
            efficiencies.append(pEff)
            

            ROOT.gPad.Update()
            graph = pEff.GetPaintedGraph()
            if graph:
                graph.GetYaxis().SetRangeUser(0.0, 1.05)
                

                thresh = thresholds[p_name]
                line = ROOT.TLine(thresh, 0, thresh, 1.05)
                line.SetLineColor(ROOT.kRed)
                line.SetLineStyle(2)
                line.SetLineWidth(2)
                line.Draw()
                draw_objects.append(line)

    c3.SaveAs(f"{output_tag}_efficiency_vs_momentum.pdf")


    print("\n" + "="*85)
    print(" PANDORA RECONSTRUCTION MASTER DIAGNOSTICS (TRUE SIGNAL IN FV)")
    print("="*85)
    print(f"Total Signal Events Evaluated : {total_events}")
    
    print("\n" + "-"*85)
    print(" 1. SLICE FRACTURING & ANATOMY")
    print("-" * 85)
    if total_events > 0:
        print(f" Perfectly Contained (1 slice) : {contained_slices} ({(contained_slices/total_events)*100:.1f}%)")
        print(f" Fractured (>1 slices)         : {fractured_slices} ({(fractured_slices/total_events)*100:.1f}%)")
        if fractured_slices > 0:
            for cat, count in fracture_modes.items():
                print(f"    -> {cat:<24} : {count} ({(count/fractured_slices)*100:.1f}%)")
                
    print("\n" + "-"*85)
    print(" 2. PFP HIERARCHY MAPPING (Pandora's internal 'Family Tree')")
    print("-" * 85)
    print(f" {'Particle':<20} | {'Total Events':<12} | {'Primary %':<12} | {'Secondary %':<12} | {'Not Reco %'}")
    print("-" * 85)
    for p_name in hierarchy_particles:
        tot = hierarchy_results[p_name]["Total"]
        if tot == 0: continue
        pct_prim = (hierarchy_results[p_name]["Primary"] / tot) * 100.0
        pct_sec  = (hierarchy_results[p_name]["Secondary"] / tot) * 100.0
        pct_miss = (hierarchy_results[p_name]["Not Reco"] / tot) * 100.0
        print(f" {p_name:<20} | {tot:<12} | {pct_prim:>8.1f} %   | {pct_sec:>8.1f} %   | {pct_miss:>8.1f} %")

    print("\n" + "-"*85)
    print(" 3. DAUGHTER TRACKING BREAKDOWN (Why is Lambda -> p + pi- 'Not Reco'?)")
    print("-" * 85)
    if total_events > 0:
        for cat, count in daughter_modes.items():
            print(f" {cat:<27} : {count} ({(count/total_events)*100:.1f}%)")

    print("="*85 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output_tag", default="reco")
    args = parser.parse_args()
    run_diagnostics(args.input, args.output_tag)