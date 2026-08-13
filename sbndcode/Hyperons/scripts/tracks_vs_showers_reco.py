#!/usr/bin/env python3
import ROOT
import argparse
import math

def get_pot(filename):
    """Extracts the total simulated POT from the subrun tree."""
    f = ROOT.TFile.Open(filename, "READ")
    if not f:
        return 1e18
    subrun_tree = f.Get("ana/subrunTree") or f.Get("subrunTree")
    pot = 0.0
    if subrun_tree:
        for entry in range(subrun_tree.GetEntries()):
            subrun_tree.GetEntry(entry)
            pot += getattr(subrun_tree, "pot", 0.0)
    f.Close()
    return pot if pot > 0 else 1e18

def pass_reco_fv(tree, slice_idx):
    """Applies the reconstructed Fiducial Volume cut to accurately evaluate cosmic rays."""
    vtx_x = tree.slice_VtxX[slice_idx]
    vtx_y = tree.slice_VtxY[slice_idx]
    vtx_z = tree.slice_VtxZ[slice_idx]
    
    if vtx_x == -999.0:
        return False
        
    in_x = 5.0 < abs(vtx_x) < 180.0
    in_y = -180.0 < vtx_y < 180.0
    in_z = 20.0 < vtx_z < 470.0
    return in_x and in_y and in_z

def process_file(filename, is_signal_file, h_sig, h_bkg, h_comb, h_prim_sig, h_prim_bkg, h_prim_comb, scale_factor):
    """
    Scans a file for the primary slice track/shower phase space.
    Fills the histograms using the POT scale_factor as a weight.
    Routes events to Signal or Background histograms based on sample stitching logic.
    """
    f = ROOT.TFile.Open(filename, "READ")
    tree = f.Get("ana/tree") or f.Get("tree")
    if not tree:
        print(f"Error: Could not find tree in {filename}")
        return 0, 0, 0

    num_entries = tree.GetEntries()
    max_val = h_sig.GetNbinsX()
    
    sig_count = 0
    rare_bkg_count = 0
    bulk_bkg_count = 0

    for entry in range(num_entries):
        tree.GetEntry(entry)

        is_true_signal = (getattr(tree, "IsAssocLambdaKPlusWithProtonPiMinus", 0) == 1)

        has_any_hyperon = False
        for pdg in tree.geant_PDG:
            if abs(pdg) in [3122, 3222, 3212, 3112]:
                has_any_hyperon = True
                break

        is_valid_for_this_file = False
        is_background_for_this_file = False
        
        if is_signal_file:
            is_valid_for_this_file = True
            if not is_true_signal:
                is_background_for_this_file = True
        else:
            if has_any_hyperon:
                continue
            is_valid_for_this_file = True
            is_background_for_this_file = True

        if not is_valid_for_this_file:
            continue

        best_slice_idx = -1
        max_nuscore = -1
        
        for idx, score in enumerate(tree.slice_NuScore):
            if score > max_nuscore:
                max_nuscore = score
                best_slice_idx = idx

        #if max_nuscore <= 0.56: # Commented out for now, just using the highest scored slice
        #    continue
        in_fv = pass_reco_fv(tree,best_slice_idx)
        if not in_fv:
            continue

        # Only proceed if we found a valid primary slice
        if best_slice_idx != -1:
            primary_slice_id = tree.slice_ID[best_slice_idx]

            t_slice_ids = getattr(tree, "track_SliceID", [])
            t_scores = getattr(tree, "track_TrackScore", [])
            t_IsPrimary = getattr(tree, "track_IsPrimary", [])
            s_IsPrimary = getattr(tree, "shower_IsPrimary", [])

            n_tracks = 0
            n_showers = 0
            n_tracks_primary = 0
            n_showers_primary = 0
            
            # Loop over all PFParticles via the track arrays
            for j, t_sid in enumerate(t_slice_ids):
                if t_sid == primary_slice_id:
                    score = t_scores[j]
                    if score > 0.5:
                        n_tracks += 1
                        if t_IsPrimary[j] == 1:
                            n_tracks_primary += 1
                    elif 0.0 <= score <= 0.5:
                        n_showers += 1
                        if s_IsPrimary[j] == 1:
                            n_showers_primary += 1

            plot_tracks = min(n_tracks, max_val - 1)
            plot_showers = min(n_showers, max_val - 1)
            plot_tracks_prim = min(n_tracks_primary, max_val  - 1 )
            plot_showers_prim = min(n_showers_primary, max_val  - 1 )


            h_comb.Fill(plot_tracks, plot_showers, scale_factor)
            h_prim_comb.Fill(plot_tracks_prim, plot_showers_prim, scale_factor)
            
            if is_signal_file and not is_background_for_this_file:
                h_sig.Fill(plot_tracks, plot_showers, scale_factor)
                h_prim_sig.Fill(plot_tracks_prim, plot_showers_prim, scale_factor)
                sig_count += 1
            elif is_signal_file and is_background_for_this_file:
                h_bkg.Fill(plot_tracks, plot_showers, scale_factor)
                h_prim_bkg.Fill(plot_tracks_prim, plot_showers_prim, scale_factor)
                rare_bkg_count += 1
            elif not is_signal_file:
                h_bkg.Fill(plot_tracks, plot_showers, scale_factor)
                h_prim_bkg.Fill(plot_tracks_prim, plot_showers_prim, scale_factor)
                bulk_bkg_count += 1

    f.Close()
    return sig_count, rare_bkg_count, bulk_bkg_count

def run_phase_space(signal_filename, background_filename, target_pot=1e21, output_tag="topology"):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0) 
    
    ROOT.gStyle.SetPaintTextFormat(".0f")

    print("\n" + "="*80)
    print(" INITIALIZING POT-SCALED PHASE SPACE (WITH SAMPLE STITCHING)")
    print("="*80)

    sig_pot = get_pot(signal_filename)
    bkg_pot = get_pot(background_filename)
    
    sig_scale = target_pot / sig_pot
    bkg_scale = target_pot / bkg_pot
    
    print(f"Target POT      : {target_pot:.2e}")
    print(f"Signal POT      : {sig_pot:.2e} (Scale Factor: {sig_scale:.4f})")
    print(f"Background POT  : {bkg_pot:.2e} (Scale Factor: {bkg_scale:.4f})")
    print("-" * 80)

    n_bins = 12
    max_val = 12

    h2_comb = ROOT.TH2F("h2_comb", f"Primary Slice: Combined (Scaled to {target_pot:.0e} POT);Number of Tracks;Number of Showers", n_bins, 0, max_val, n_bins, 0, max_val)
    h2_sig = ROOT.TH2F("h2_sig", f"Primary Slice: True #Lambda^{{0}} Signal;Number of Tracks;Number of Showers", n_bins, 0, max_val, n_bins, 0, max_val)
    h2_bkg = ROOT.TH2F("h2_bkg", f"Primary Slice: Total Background;Number of Tracks;Number of Showers", n_bins, 0, max_val, n_bins, 0, max_val)

    h2_prim_comb = ROOT.TH2F("h_prim_comb", f"Primary-Only: Combined (Scaled to {target_pot:.0e} POT);Number of Primary Tracks;Number of Primary Showers", n_bins, 0, max_val, n_bins, 0, max_val)
    h2_prim_sig  = ROOT.TH2F("h_prim_sig", f"Primary-Only: True #Lambda^{{0}} Signal;Number of Primary Tracks;Number of Primary Showers", n_bins, 0, max_val, n_bins, 0, max_val)
    h2_prim_bkg  = ROOT.TH2F("h_prim_bkg", f"Primary-Only: Total Background;Number of Primary Tracks;Number of Primary Showers", n_bins, 0, max_val, n_bins, 0, max_val)
    
    
    h2_comb.SetMarkerSize(0.85)
    h2_sig.SetMarkerSize(0.85)
    h2_bkg.SetMarkerSize(0.85)
    h2_prim_comb.SetMarkerSize(0.85)
    h2_prim_sig.SetMarkerSize(0.85)
    h2_prim_bkg.SetMarkerSize(0.85)
    

    print("Processing Signal Sample (Extracting Signal + Background Hyperons)...")
    sig_found, rare_bkg_found, _ = process_file(signal_filename, True, h2_sig, h2_bkg, h2_comb, h2_prim_sig, h2_prim_bkg, h2_prim_comb, sig_scale)
    
    print("Processing Background Sample (Extracting Bulk Background)...")
    _, _, bulk_bkg_found = process_file(background_filename, False, h2_sig, h2_bkg, h2_comb, h2_prim_sig, h2_prim_bkg, h2_prim_comb, bkg_scale)

    c1 = ROOT.TCanvas("c1", "Topology Phase Space", 1800, 600)
    c1.Divide(3, 1)
    
    c1.cd(1)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    h2_comb.Draw("COLZ TEXT")
    
    c1.cd(2)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    h2_sig.Draw("COLZ TEXT")
    
    c1.cd(3)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    h2_bkg.Draw("COLZ TEXT")

    c1.SaveAs(f"{output_tag}_stitched_phasespace.pdf")

    c1.cd(1)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    h2_prim_comb.Draw("COLZ TEXT")
    
    c1.cd(2)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    h2_prim_sig.Draw("COLZ TEXT")
    
    c1.cd(3)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    h2_prim_bkg.Draw("COLZ TEXT")

    c1.SaveAs(f"{output_tag}_stitched_phasespace_primary_PFPs_only.pdf")

    # Optimisation, dumb
    print("\n" + "="*80)
    print(" TOPOLOGICAL CUT OPTIMIZER (Ranked by S / sqrt(S + B))")
    print("="*80)
    
    total_sig_pot = h2_sig.Integral()
    total_sig_pot_prim = h2_prim_sig.Integral()
    results = []
    results_prim = []

    for min_t in range(1, 7):     # Scan minimum tracks from 1 to 6
        for max_s in range(0, 4): # Scan maximum showers from 0 to 3
            
            x_start = min_t + 1   # ROOT bin for min_t
            x_end = h2_sig.GetNbinsX()
            
            y_start = 1           # ROOT bin for 0 showers
            y_end = max_s + 1     # ROOT bin for max_s
            
            s = h2_sig.Integral(x_start, x_end, y_start, y_end)
            b = h2_bkg.Integral(x_start, x_end, y_start, y_end)
            
            if s + b > 0:
                significance = s / math.sqrt(s + b)
                purity = (s / (s + b)) * 100.0
                efficiency = (s / total_sig_pot) * 100.0 if total_sig_pot > 0 else 0
                
                results.append({
                    "cut_str": f"Tracks >= {min_t} & Showers <= {max_s}",
                    "s": s,
                    "b": b,
                    "sig": significance,
                    "pur": purity,
                    "eff": efficiency
                })

    for min_t in range(1, 7):     # Scan minimum tracks from 1 to 6
        for max_s in range(0, 4): # Scan maximum showers from 0 to 3
            
            x_start = min_t + 1   # ROOT bin for min_t
            x_end = h2_prim_sig.GetNbinsX()
            
            y_start = 1           # ROOT bin for 0 showers
            y_end = max_s + 1     # ROOT bin for max_s
            
            s = h2_prim_sig.Integral(x_start, x_end, y_start, y_end)
            b = h2_prim_bkg.Integral(x_start, x_end, y_start, y_end)
            
            if s + b > 0:
                significance = s / math.sqrt(s + b)
                purity = (s / (s + b)) * 100.0
                efficiency = (s / total_sig_pot_prim) * 100.0 if total_sig_pot_prim > 0 else 0
                
                results_prim.append({
                    "cut_str": f"Tracks >= {min_t} & Showers <= {max_s}",
                    "s": s,
                    "b": b,
                    "sig": significance,
                    "pur": purity,
                    "eff": efficiency
                })

    results.sort(key=lambda x: x["sig"], reverse=True)
    results_prim.sort(key=lambda x: x["sig"], reverse=True)

    print(f" {'Proposed Cut':<25} | {'Signal':<8} | {'Background':<10} | {'Eff %':<6} | {'Pur %':<6} | {'Significance'}")
    print("-" * 80)
    for i, res in enumerate(results[:5]): # Print top 5 cuts
        star = " *" if i == 0 else ""     # Highlight the best one
        print(f" {res['cut_str']:<25} | {res['s']:<8.1f} | {res['b']:<10.1f} | {res['eff']:<6.1f} | {res['pur']:<6.1f} | {res['sig']:.2f}{star}")

    print("-"*80)
    print("FOR PRIMARY PFPS ONLY")
    print("-" * 80)
    for i, res in enumerate(results_prim[:5]): # Print top 5 cuts
        star = " *" if i == 0 else ""     # Highlight the best one
        print(f" {res['cut_str']:<25} | {res['s']:<8.1f} | {res['b']:<10.1f} | {res['eff']:<6.1f} | {res['pur']:<6.1f} | {res['sig']:.2f}{star}")


    print("\n" + "="*80)
    print(" TOPOLOGICAL PHASE SPACE SUMMARY")
    print("="*80)
    print(f"Raw Signal Primary Slices            : {sig_found}")
    print(f"Raw Background Hyperons (from Sig)   : {rare_bkg_found}")
    print(f"Raw Bulk Background (from Bkg)       : {bulk_bkg_found}")

    print("\n" + "="*80)

    print(f"-> Plots successfully stitched and scaled to {target_pot:.2e} POT.")
    print("="*80 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--signal", required=True, help="Filtered Hyperon ROOT file")
    parser.add_argument("-b", "--background", required=True, help="Inclusive Background ROOT file")
    parser.add_argument("-p", "--pot", type=float, default=1e21, help="Target POT to scale to")
    parser.add_argument("-o", "--output_tag", default="topology")
    args = parser.parse_args()
    run_phase_space(args.signal, args.background, args.pot, args.output_tag)