import ROOT
import glob
import argparse

ROOT.gInterpreter.Declare('#include "larcoreobj/SummaryData/POTSummary.h"')

def main(file_glob):
    ROOT.ROOT.EnableImplicitMT()

    files = glob.glob(file_glob, recursive=True)
    if not files:
        raise FileNotFoundError(f"No files matched the glob: {file_glob}")

    chain_subruns = ROOT.TChain("SubRuns")
    chain_events = ROOT.TChain("Events")

    valid_files_count = 0
    excluded_files = []

    for f in files:
        try:
            tf = ROOT.TFile.Open(f, "READ")
        except OSError as e:
            # Incase file is truncated/not finished copying
            excluded_files.append(f"{f} (Currently Copying / Truncated)")
            continue

        # Don't count genuinely corrupted files
        if not tf or tf.IsZombie():
            excluded_files.append(f"{f} (Zombie / Corrupted)")
            if tf: 
                tf.Close()
            continue
        
        if not tf.GetListOfKeys().Contains("SubRuns") or not tf.GetListOfKeys().Contains("Events"):
            excluded_files.append(f"{f} (Missing 'SubRuns' or 'Events' tree)")
            tf.Close()
            continue

        tf.Close()
        chain_subruns.Add(f)
        chain_events.Add(f)
        valid_files_count += 1

    if valid_files_count == 0:
        print("\n--- ERROR ---")
        print("No valid files could be opened.")
        if excluded_files:
            print("Excluded files:")
            for ef in excluded_files:
                print(f"  - {ef}")
        raise RuntimeError("Zero valid files to process.")

    df_subruns = ROOT.RDataFrame(chain_subruns)
    df_events = ROOT.RDataFrame(chain_events)

    num_events_ptr = df_events.Count()

    df_subruns = df_subruns.Alias("POTObj", "sumdata::POTSummary_generator__GenieGen.obj")
    df_subruns = df_subruns.Define("POTValue", "POTObj.totpot")
    totalPOT_ptr = df_subruns.Sum("POTValue")

    num_events = num_events_ptr.GetValue()
    totalPOT = totalPOT_ptr.GetValue()

    print("\nFinished getting values")
    print("========================")
    print(f"For files in: {file_glob}")
    print(f"Processed files: {valid_files_count}/{len(files)}")
    print(f"POT: {totalPOT}")
    print(f"Number of events is: {num_events}")
    print("========================")

    if excluded_files:
        print(f"\n[WARNING] Excluded {len(excluded_files)} file(s) because they could not be opened:")
        for ef in excluded_files:
            print(f"  - {ef}")
    else:
        print("\nAll matched files were processed successfully!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_glob', required=True, help="Input glob",
                        default="/pnfs/sbnd/scratch/users/tjones/reco2_test/out/*/reco1-detsim-g4-gen-*.root")
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    main(args.input_glob)
