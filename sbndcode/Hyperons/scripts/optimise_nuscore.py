#!/usr/bin/env python3
import ROOT
import argparse
import math

def get_pot(filename):
    """Extracts the total simulated POT from the subrun tree."""
    f = ROOT.TFile.Open(filename, "READ")
    subrun_tree = f.Get("ana/subrunTree") or f.Get("subrunTree")
    pot = 0.0
    if subrun_tree:
        for entry in range(subrun_tree.GetEntries()):
            subrun_tree.GetEntry(entry)
            pot += getattr(subrun_tree, "pot", 0.0)
    f.Close()
    return pot if pot > 0 else 1e18

def pass_true_fv(tree):
    """Checks if the TRUE neutrino interaction is within the Fiducial Volume."""
    if len(tree.TrueVtxX) > 0:
        vx, vy, vz = tree.TrueVtxX[0], tree.TrueVtxY[0], tree.TrueVtxZ[0]
        in_x = (-190.0 < vx < -5.0) or (5.0 < vx < 190.0)
        in_y = (-180.0 < vy < 180.0)
        in_z = (20.0 < vz < 470.0)
        return in_x and in_y and in_z
    return False

def pass_reco_fv(tree, slice_idx):
    """Checks if the RECONSTRUCTED slice is within the FV to remove dirt/edge cosmics."""
    vtx_x = tree.slice_VtxX[slice_idx]
    vtx_y = tree.slice_VtxY[slice_idx]
    vtx_z = tree.slice_VtxZ[slice_idx]
    if vtx_x == -999.0: return False
    in_x = 5.0 < abs(vtx_x) < 180.0
    in_y = -180.0 < vtx_y < 180.0
    in_z = 20.0 < vtx_z < 470.0
    return in_x and in_y and in_z

def run_2d_optimization(signal_file, bkg_file, target_pot=1e21):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    
    print("\n" + "="*90)
    print(" 2D PRE-SELECTION OPTIMIZATION (NuScore & OpT0 Flash Match)")
    print("="*90)

    sig_pot = get_pot(signal_file)
    bkg_pot = get_pot(bkg_file)
    
    sig_scale = target_pot / sig_pot
    bkg_scale = target_pot / bkg_pot
    
    print(f"Target POT          : {target_pot:.2e}")
    print(f"Signal Scale Factor : {sig_scale:.4f} (from {sig_pot:.2e} POT)")
    print(f"Bkg Scale Factor    : {bkg_scale:.4f} (from {bkg_pot:.2e} POT)")
    print("-" * 90)

    nx, ny = 50, 50
    h2_sig = ROOT.TH2F("h2_sig", "True #Lambda^{0} Signal;Max NuScore in Event;Optical Flash Match Score (OpT0)", nx, 0, 1.0, ny, 0, 500)
    h2_bkg = ROOT.TH2F("h2_bkg", "Total Background;Max NuScore in Event;Optical Flash Match Score (OpT0)", nx, 0, 1.0, ny, 0, 500)
    h2_cosmic = ROOT.TH2F("h2_cosmic", "Pure Cosmics (True Origin = 2);Max NuScore in Event;Optical Flash Match Score (OpT0)", nx, 0, 1.0, ny, 0, 500)
    
    h1_nu_sig = ROOT.TH1F("h1_nu_sig", "Signal", nx, 0, 1.0)
    h1_nu_beam = ROOT.TH1F("h1_nu_beam", "Beam Background & Dirt", nx, 0, 1.0)
    h1_nu_cosmic = ROOT.TH1F("h1_nu_cosmic", "Pure Cosmics", nx, 0, 1.0)
    
    h1_opt0_sig = ROOT.TH1F("h1_opt0_sig", "Signal", ny, 0, 500)
    h1_opt0_beam = ROOT.TH1F("h1_opt0_beam", "Beam Background & Dirt", ny, 0, 500)
    h1_opt0_cosmic = ROOT.TH1F("h1_opt0_cosmic", "Pure Cosmics", ny, 0, 500)
    
    for h in [h1_nu_sig, h1_opt0_sig]:
        h.SetFillColor(ROOT.kAzure+1)
        h.SetLineColor(ROOT.kAzure+1) 
    for h in [h1_nu_beam, h1_opt0_beam]:
        h.SetFillColor(ROOT.kOrange+1)
        h.SetLineColor(ROOT.kOrange+1)
    for h in [h1_nu_cosmic, h1_opt0_cosmic]:
        h.SetFillColor(ROOT.kGray+1)
        h.SetLineColor(ROOT.kGray+2) 
    
    base_signal = 0
    base_bkg = 0

    print("Processing Signal Sample (Extracting Signal in FV)...")
    f_sig = ROOT.TFile.Open(signal_file, "READ")
    t_sig = f_sig.Get("ana/tree") or f_sig.Get("tree")
    
    for entry in range(t_sig.GetEntries()):
        t_sig.GetEntry(entry)
        
        is_signal = (getattr(t_sig, "IsAssocLambdaKPlusWithProtonPiMinus", 0) == 1)
        if not is_signal or not pass_true_fv(t_sig):
            continue
            
        base_signal += 1
            
        # Find the slice with the maximum NuScore
        best_idx = -1
        max_nu = -1.0
        for i, score in enumerate(t_sig.slice_NuScore):
            if score > max_nu:
                max_nu = score
                best_idx = i
                
        # ONLY plot if it passes the Fiducial Volume cut!
        if best_idx >= 0 and pass_reco_fv(t_sig, best_idx):
            opt0 = t_sig.slice_Opt0Score[best_idx]
            h2_sig.Fill(max_nu, opt0, sig_scale)
            
            h1_nu_sig.Fill(max_nu, sig_scale)
            h1_opt0_sig.Fill(opt0, sig_scale)
                
    f_sig.Close()

    print("Processing Background Sample (Extracting Bulk Beam & Cosmics)...")
    f_bkg = ROOT.TFile.Open(bkg_file, "READ")
    t_bkg = f_bkg.Get("ana/tree") or f_bkg.Get("tree")
    
    for entry in range(t_bkg.GetEntries()):
        t_bkg.GetEntry(entry)
        
        # Do not double count true signal from the inclusive file
        is_signal = (getattr(t_bkg, "IsAssocLambdaKPlusWithProtonPiMinus", 0) == 1)
        if is_signal and pass_true_fv(t_bkg): 
            continue 
            
        base_bkg += 1
            
        best_idx = -1
        max_nu = -1.0
        for i, score in enumerate(t_bkg.slice_NuScore):
            if score > max_nu:
                max_nu = score
                best_idx = i
                
        if best_idx >= 0 and pass_reco_fv(t_bkg, best_idx):
            opt0 = t_bkg.slice_Opt0Score[best_idx]
            
            is_cosmic = False
            origin_array = getattr(t_bkg, "slice_TrueOrigin", [])
            if len(origin_array) > best_idx and origin_array[best_idx] == 2:
                is_cosmic = True
            
            h2_bkg.Fill(max_nu, opt0, bkg_scale)
            
            if is_cosmic:
                h2_cosmic.Fill(max_nu, opt0, bkg_scale)
                h1_nu_cosmic.Fill(max_nu, bkg_scale)
                h1_opt0_cosmic.Fill(opt0, bkg_scale)
            else:
                h1_nu_beam.Fill(max_nu, bkg_scale)
                h1_opt0_beam.Fill(opt0, bkg_scale)
            
    f_bkg.Close()

    base_sig_scaled = base_signal * sig_scale
    base_bkg_scaled = base_bkg * bkg_scale
    
    results = []
    
    h2_significance = ROOT.TH2F("h2_sig_map", "Significance Heatmap (S / #sqrt{S+B});NuScore Threshold (Cut > X);OpT0 Threshold (Cut > Y)", nx, 0, 1.0, ny, 0, 500)

    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            nu_cut = h2_sig.GetXaxis().GetBinLowEdge(i)
            opt0_cut = h2_sig.GetYaxis().GetBinLowEdge(j)
            
            S = h2_sig.Integral(i, nx + 1, j, ny + 1)
            B = h2_bkg.Integral(i, nx + 1, j, ny + 1)
            C = h2_cosmic.Integral(i, nx + 1, j, ny + 1) 
            
            eff = S / base_sig_scaled if base_sig_scaled > 0 else 0
            pur = S / (S + B) if (S + B) > 0 else 0
            
            sig_cosmic = S / math.sqrt(S + C) if (S + C) > 0 else 0
            
            h2_significance.SetBinContent(i, j, sig_cosmic)
            
            results.append({
                "nu": nu_cut,
                "opt0": opt0_cut,
                "S": S,
                "B": B,
                "C_surv": C,
                "eff": eff * 100,
                "pur": pur * 100,
                "sig_cosmic": sig_cosmic
            })

    # Sort by the new Cosmic-only Significance
    results.sort(key=lambda x: x["sig_cosmic"], reverse=True)

    print("\n" + "="*115)
    print(" TOP 10 PRE-SELECTION COMBINATIONS (Ranked by Cosmic-Only Significance)")
    print("="*115)
    print(f" {'NuScore >':<10} | {'OpT0 >':<8} | {'Signal':<8} | {'Total Bkg':<10} | {'Sig Eff %':<10} | {'Cosmics Surv':<12} | {'Cosmic Sig':<12}")
    print("-" * 115)
    
    for res in results[:10]:
        print(f" {res['nu']:<10.2f} | {res['opt0']:<8.1f} | {res['S']:<8.1f} | {res['B']:<10.1f} | {res['eff']:<10.2f} | {res['C_surv']:<12.1f} | {res['sig_cosmic']:<12.2f}")
        
    print("\n" + "="*115)
    best = results[0]
    print(f" Optimal 2D Cut: NuScore > {best['nu']:.2f}, OpT0 > {best['opt0']:.0f}")
    print("="*115 + "\n")

    c1 = ROOT.TCanvas("c1", "2D Pre-Selection Optimization", 1800, 1200)
    c1.Divide(2, 2)
    
    ROOT.gStyle.SetPalette(ROOT.kRainBow) # Good contrast palette
    
    # 1. Signal Distribution
    c1.cd(1)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    if h2_sig.GetMaximum() > 0: h2_sig.SetMinimum(0.1)
    h2_sig.Draw("COLZ")
    
    # 2. Background Distribution
    c1.cd(2)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    if h2_bkg.GetMaximum() > 0: h2_bkg.SetMinimum(0.1)
    h2_bkg.Draw("COLZ")
    
    # 3. Cosmic Distribution
    c1.cd(3)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetLogz(1)
    if h2_cosmic.GetMaximum() > 0: h2_cosmic.SetMinimum(0.1)
    h2_cosmic.Draw("COLZ")
    
    # 4. Significance Heatmap
    c1.cd(4)
    ROOT.gPad.SetRightMargin(0.15)
    h2_significance.Draw("COLZ")
    
    marker = ROOT.TMarker(best["nu"], best["opt0"], 29)
    marker.SetMarkerColor(ROOT.kRed)
    marker.SetMarkerSize(3)
    marker.Draw()
    
    label = ROOT.TLatex()
    label.SetTextSize(0.04)
    label.SetTextColor(ROOT.kRed)
    label.DrawLatex(best["nu"] + 0.05, best["opt0"] + 20, f"Optimal: Nu>{best['nu']:.2f}, OpT0>{best['opt0']:.0f}")

    c1.SaveAs("preselection_2d_optimization.pdf")
    print("[*] Saved 2D diagnostic plots to preselection_2d_optimization.pdf")

    c2 = ROOT.TCanvas("c2", "1D Pre-Selection Projections", 1400, 600)
    c2.Divide(2, 1)

    # From Bears CC inclusive
    chosen_nu = 0.50
    chosen_opt0 = 320.0

    for h in [h1_nu_cosmic, h1_nu_beam, h1_nu_sig, h1_opt0_cosmic, h1_opt0_beam, h1_opt0_sig]:
        if h.Integral() > 0:
            h.Scale(1.0 / h.Integral())
        h.SetLineWidth(3)
        h.SetFillStyle(0) # Remove solid fill so we can see overlapping lines

    # Plot 1: NuScore
    c2.cd(1)
    ROOT.gPad.SetLogy(0) 
    
    max_y_nu = max(h1_nu_cosmic.GetMaximum(), h1_nu_beam.GetMaximum(), h1_nu_sig.GetMaximum())
    h1_nu_cosmic.SetMaximum(max_y_nu * 1.3)
    h1_nu_cosmic.SetTitle("Max NuScore (Area Normalized Shapes);Pandora NuScore;Probability Density")
    
    h1_nu_cosmic.Draw("HIST")
    h1_nu_beam.Draw("HIST SAME")
    h1_nu_sig.Draw("HIST SAME")
    

    l_nu = ROOT.TLine(chosen_nu, 0, chosen_nu, max_y_nu * 1.1)
    l_nu.SetLineColor(ROOT.kRed)
    l_nu.SetLineWidth(3)
    l_nu.SetLineStyle(2)
    l_nu.Draw()

    leg_nu = ROOT.TLegend(0.40, 0.70, 0.88, 0.88)
    leg_nu.SetBorderSize(0)
    leg_nu.AddEntry(h1_nu_sig, "True #Lambda^{0} Signal", "l")
    leg_nu.AddEntry(h1_nu_beam, "Beam Background & Dirt", "l")
    leg_nu.AddEntry(h1_nu_cosmic, "Pure Cosmics (Origin=2)", "l")
    leg_nu.AddEntry(l_nu, f"Chosen Cut (> {chosen_nu})", "l")
    leg_nu.Draw()

    # Plot 2: OpT0
    c2.cd(2)
    ROOT.gPad.SetLogy(0)
    
    max_y_opt0 = max(h1_opt0_cosmic.GetMaximum(), h1_opt0_beam.GetMaximum(), h1_opt0_sig.GetMaximum())
    h1_opt0_cosmic.SetMaximum(max_y_opt0 * 1.3)
    h1_opt0_cosmic.SetTitle("Optical Flash Match (Area Normalized Shapes);OpT0 Score;Probability Density")
    
    h1_opt0_cosmic.Draw("HIST")
    h1_opt0_beam.Draw("HIST SAME")
    h1_opt0_sig.Draw("HIST SAME")

    l_opt0 = ROOT.TLine(chosen_opt0, 0, chosen_opt0, max_y_opt0 * 1.1)
    l_opt0.SetLineColor(ROOT.kRed)
    l_opt0.SetLineWidth(3)
    l_opt0.SetLineStyle(2)
    l_opt0.Draw()
    
    leg_opt0 = ROOT.TLegend(0.40, 0.70, 0.88, 0.88)
    leg_opt0.SetBorderSize(0)
    leg_opt0.AddEntry(h1_opt0_sig, "True #Lambda^{0} Signal", "l")
    leg_opt0.AddEntry(h1_opt0_beam, "Beam Background & Dirt", "l")
    leg_opt0.AddEntry(h1_opt0_cosmic, "Pure Cosmics (Origin=2)", "l")
    leg_opt0.AddEntry(l_opt0, f"Chosen Cut (> {chosen_opt0:.0f})", "l")
    leg_opt0.Draw()

    c2.SaveAs("preselection_1d_stacks.pdf")
    print("[*] Saved 1D shape plots to preselection_1d_stacks.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--signal", required=True, help="Filtered Hyperon ROOT file")
    parser.add_argument("-b", "--background", required=True, help="Inclusive Background ROOT file")
    parser.add_argument("-p", "--pot", type=float, default=1e21)
    args = parser.parse_args()
    
    run_2d_optimization(args.signal, args.background, args.pot)