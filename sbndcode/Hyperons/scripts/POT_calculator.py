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

    for f in files:
        chain_subruns.Add(f)
        chain_events.Add(f)

    if chain_subruns.GetNtrees() == 0:
        raise RuntimeError("TChain could not open the 'SubRuns' tree in the provided files.")

    df_subruns = ROOT.RDataFrame(chain_subruns)
    df_events = ROOT.RDataFrame(chain_events)

    num_events_ptr = df_events.Count()

    df_subruns = df_subruns.Alias("POTObj", "sumdata::POTSummary_generator__GenieGen.obj")
    df_subruns = df_subruns.Define("POTValue", "POTObj.totpot")
    totalPOT_ptr = df_subruns.Sum("POTValue")


    num_events = num_events_ptr.GetValue()
    totalPOT = totalPOT_ptr.GetValue()

    print("Finished getting values")

    print(f"For files in {file_glob}")
    print(f"POT: {totalPOT}")
    print(f"Number of events is: {num_events}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_glob', required=True, help="Input glob",
                        default="/pnfs/sbnd/scratch/users/tjones/reco2_test/out/*/reco1-detsim-g4-gen-*.root")
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    main(args.input_glob)
