import ROOT
import argparse
import sys
import math
import array

def main(sig_file, bkg_file):
    print(f"Loading Signal:     {sig_file}")
    print(f"Loading Background: {bkg_file}")
    
    df_sig = ROOT.RDataFrame("ana/tree", sig_file)
    df_bkg = ROOT.RDataFrame("ana/tree", bkg_file)


    weight_sig = 4.18
    weight_bkg = 50.39


    fv_cut = """
        std::abs(SliceVertexX.back()) > 5.0 && std::abs(SliceVertexX.back()) < 180.0 &&
        SliceVertexY.back() > -180.0 && SliceVertexY.back() < 180.0 &&
        SliceVertexZ.back() > 20.0 && SliceVertexZ.back() < 470.0
    """
    
    # Apply baseline cuts to both datasets
    df_sig_base = df_sig.Filter("NTracks > 0").Filter(fv_cut).Filter("NTracks >= 4")
    df_bkg_base = df_bkg.Filter("NTracks > 0").Filter(fv_cut).Filter("NTracks >= 4")

    # SIGNAL: Must be the actual Lambda-producing neutrino slice, not a cosmic overlay
    df_signal = df_sig_base.Filter("SliceIsTrueNeutrino.back() == 1", "Genuine Lambda Neutrino")
    
    # BACKGROUND: The entire inclusive file (cosmics + boring neutrinos)
    df_background = df_bkg_base

    # Get the TOTAL signal before any NuScore cut is applied (for Efficiency math)
    total_raw_S = df_signal.Count().GetValue()
    total_scaled_S = total_raw_S * weight_sig

    print("\n" + "="*85)
    print(f"{'NuScore (>x)':<12} | {'Raw S':<8} | {'Raw B':<8} | {'Scaled S':<10} | {'Scaled B':<10} | {'Significance':<10} | {'efficiency':<10} | {'Purity':<10} | {'Raw Cosmic B':<10} | {'Raw Beam B':<10}")
    print("-" * 115)

    best_sig = 0.0
    best_cut = 0.0

    cuts = array.array('d')
    sigs = array.array('d')

    # Scan NuScore thresholds from 0.00 to 1.00
    for i in range(51):
        threshold = i * 0.02
        
        # Count passing events
        raw_S = df_signal.Filter(f"SliceNuScore.back() > {threshold}").Count().GetValue()
        raw_B = df_background.Filter(f"SliceNuScore.back() > {threshold}").Count().GetValue()

        raw_bkg_cosmics = df_background.Filter(f"SliceNuScore.back() > {threshold} && SliceIsTrueNeutrino.back() == 0").Count().GetValue()
        raw_bkg_beam    = df_background.Filter(f"SliceNuScore.back() > {threshold} && SliceIsTrueNeutrino.back() == 1").Count().GetValue()


        # Apply POT scaling
        scaled_S = raw_S * weight_sig
        scaled_B = raw_B * weight_bkg
        
        # Calculate Math Metrics
        if (scaled_S + scaled_B) > 0:
            significance = scaled_S / math.sqrt(scaled_S + scaled_B)
            purity = (scaled_S / (scaled_S + scaled_B)) * 100.0
        else:
            significance = 0.0
            purity = 0.0
            
        if total_scaled_S > 0:
            efficiency = (scaled_S / total_scaled_S) * 100.0
        else:
            efficiency = 0.0
            
        print(f"> {threshold:<10.2f} | {raw_S:<8} | {raw_B:<8} | {scaled_S:<10.2f} | {scaled_B:<10.2f} | {significance:<10.3f} | {efficiency:<10.3f}% | {purity:<10.3f}% | {raw_bkg_cosmics:<10.3f} | {raw_bkg_beam:.3f}")

        # Save points for plotting
        cuts.append(threshold)
        sigs.append(significance)

        # Track the absolute maximum significance
        if significance > best_sig:
            best_sig = significance
            best_cut = threshold

    print("="*85)
    print(f"OPTIMAL NUSCORE CUT: > {best_cut:.2f} (Max Significance: {best_sig:.3f})")
    print("="*85 + "\n")

    c1 = ROOT.TCanvas("c1", "NuScore Optimization", 800, 600)
    c1.SetGrid()
    
    graph = ROOT.TGraph(len(cuts), cuts, sigs)
    graph.SetTitle("NuScore Optimization for Lambda Selection;Pandora NuScore Threshold (> x);Significance (S / #sqrt{S+B})")
    graph.SetLineColor(ROOT.kBlue)
    graph.SetLineWidth(2)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue)
    
    graph.Draw("APL")
    c1.SaveAs("NuScore_Significance.pdf")
    print("Saved optimization curve to NuScore_Significance.pdf!")

    # Are the cosmics dead?
    print("\n" + "="*85)
    print("BACKGROUND BREAKDOWN AT OPTIMAL CUT")
    print("="*85)
    
    # Get all background events that survived the best cut
    df_surviving_bkg = df_background.Filter(f"SliceNuScore.back() > {best_cut}")
    
    # Split them using the Truth Flag
    raw_bkg_cosmics = df_surviving_bkg.Filter("SliceIsTrueNeutrino.back() == 0").Count().GetValue()
    raw_bkg_beam    = df_surviving_bkg.Filter("SliceIsTrueNeutrino.back() == 1").Count().GetValue()
    
    total_raw_survivors = raw_bkg_cosmics + raw_bkg_beam
    
    if total_raw_survivors > 0:
        pct_cosmic = (raw_bkg_cosmics / total_raw_survivors) * 100.0
        pct_beam   = (raw_bkg_beam / total_raw_survivors) * 100.0
        
        print(f"Total Surviving Background (Raw): {total_raw_survivors}")
        print(f" -> True Cosmics:         {raw_bkg_cosmics} ({pct_cosmic:.1f}%)")
        print(f" -> True Generic Beam:    {raw_bkg_beam} ({pct_beam:.1f}%)")
    print("="*85 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Optimize Pandora NuScore for Lambdas')
    parser.add_argument('-s', '--signal', required=True, help="Signal files wildcard (Must be in quotes!)")
    parser.add_argument('-b', '--background', required=True, help="Background files wildcard (Must be in quotes!)")
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    main(args.signal, args.background)