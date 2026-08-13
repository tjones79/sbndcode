#!/usr/bin/env python3
import ROOT
import argparse

def analyze_slicing_inclusive(input_filename, output_filename="slicing_metrics_inclusive.root", hit_threshold=5):
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("ana/tree") or input_file.Get("tree")
    
    if not tree:
        print("Error: Could not find tree in file.")
        return

    output_file = ROOT.TFile.Open(output_filename, "RECREATE")
    
    # Histograms to assess the "Inclusive" Lambda K+ Slicing
    h_nslices = ROOT.TH1F("h_nslices", "N Slices with Nu Hits (Inclusive Lambda K+);Slices;Events", 6, -0.5, 5.5)
    h_purity = ROOT.TH1F("h_purity", "Primary Slice Purity;Purity;Slices", 40, 0, 1.0)
    h_completeness = ROOT.TH1F("h_completeness", "Primary Slice Completeness;Completeness;Slices", 40, 0, 1.0)

    print(f"Assessing slicing for inclusive IsAssocLambdaKPlus channel...")
    
    n_analyzed = 0
    for event in tree:
        # Key: Filter only for the inclusive Lambda K+ signal
        if event.IsAssocLambdaKPlus != 1 or event.IsInsideFV < 1:
            continue
            
        n_analyzed += 1
        total_true_hits = event.event_TotalTrueNuHits
        print(f"Total true hits are: {total_true_hits}")
        if total_true_hits <= 0: continue

        # Identify the slice with the most signal (Primary Slice)
        max_true_hits = -1
        primary_slice = None
        slices_with_signal_count = 0
        
        for i in range(len(event.slice_ID)):
            true_hits = event.slice_TrueNuHits[i]
            if true_hits >= hit_threshold:
                slices_with_signal_count += 1
                if true_hits > max_true_hits:
                    max_true_hits = true_hits
                    primary_slice = i

        h_nslices.Fill(slices_with_signal_count)
        print(slices_with_signal_count)
        
        if primary_slice is not None:
            total_hits = event.slice_TotalHits[primary_slice]
            purity = event.slice_TrueNuHits[primary_slice] / total_hits if total_hits > 0 else 0
            completeness = event.slice_TrueNuHits[primary_slice] / total_true_hits
            
            h_purity.Fill(purity)
            h_completeness.Fill(completeness)

    print(f"Finished. Analyzed {n_analyzed} inclusive signal events.")
    output_file.Write()
    output_file.Close()
    input_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    args = parser.parse_args()
    analyze_slicing_inclusive(args.input)