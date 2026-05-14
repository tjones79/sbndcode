import ROOT
import argparse
import sys

# 1. The Event-Level Truth Macro (Using TruePart arrays)
ROOT.gInterpreter.Declare("""
#include <ROOT/RVec.hxx>
using ROOT::VecOps::RVec;

bool is_true_event(const RVec<int>& pdg, const RVec<int>& mother) {
    bool has_muon = false;
    bool has_kaon = false;
    bool has_lambda = false;

    // We loop over ALL particles Geant4 created
    for (size_t i = 0; i < pdg.size(); i++) {
        // Did it come from the neutrino?
        if (mother[i] == 10000000) {
            if (pdg[i] == 13) has_muon = true;
            else if (pdg[i] == 321) has_kaon = true;
            else if (pdg[i] == 3122) has_lambda = true;
        }
    }
    return (has_muon && has_kaon && has_lambda);
}
""")

def main(file_name):
    f = ROOT.TFile.Open(file_name, "READ")
    if not f or f.IsZombie():
        print(f"Error: Could not open {file_name}!")
        sys.exit(1)

    tree = f.Get("ana/tree")
    if not tree:
        print("Error: Could not find ana/tree in the file!")
        print("Check this is the right file you silly bugger")
        sys.exit(1)

    df = ROOT.RDataFrame("ana/tree", file_name)

    df_base = df.Filter("NTracks > 0", "Events with > 0 Nu Tracks")
    
    fv_cut = """
        std::abs(SliceVertexX.back()) > 5.0 && std::abs(SliceVertexX.back()) < 180.0 &&
        SliceVertexY.back() > -180.0 && SliceVertexY.back() < 180.0 &&
        SliceVertexZ.back() > 20.0 && SliceVertexZ.back() < 470.0
    """

    df_base = df_base.Filter(fv_cut, "Fiducial Volume")


    c1 = ROOT.TCanvas()
    nuscore_histo = df_base.Histo1D("SliceNuScore")
    nuscore_histo.Draw()
    c1.SaveAs("nuscore.pdf")

    df_selected = df_base.Filter("NTracks >= 4", "NTracks >= 4") # All neutrino slices with 4 or more tracks

    report = df_selected.Report()
    report.Print()
  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='make_vertex_plots',
                    description='Makes a plot of the primary vertices in x,y and z',
                    epilog='Pure silliness')
    parser.add_argument('-i', '--input_filename', required=True, help="Input file init m8")
    args = parser.parse_args()
    ROOT.gROOT.SetBatch(True)
    main(args.input_filename)