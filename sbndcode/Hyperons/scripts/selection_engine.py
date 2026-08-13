#!/usr/bin/env python3
import ROOT
import argparse
import math

def pass_true_fv(tree):
    """Manually checks if the true neutrino interaction is within the FV."""
    if len(tree.TrueVtxX) > 0:
        vx, vy, vz = tree.TrueVtxX[0], tree.TrueVtxY[0], tree.TrueVtxZ[0]
        in_x = (-190.0 < vx < -5.0) or (5.0 < vx < 190.0)
        in_y = (-180.0 < vy < 180.0)
        in_z = (20.0 < vz < 470.0)
        return in_x and in_y and in_z
    return False

def get_category(tree):
    """Categorises the event into Signal, Cosmic/Dirt, or specific beam backgrounds."""
    in_fv = pass_true_fv(tree)
    is_signal = (getattr(tree, "IsAssocLambdaKPlusWithProtonPiMinus", 0) == 1)
    
    if is_signal and in_fv:
        return "Signal"
        
    if not in_fv:
        return "Cosmic/Dirt"
        
    # If it is in the FV but isn't signal, it's a Beam Background. Identify the mode!
    mode = -1
    if hasattr(tree, "nu_numode") and len(tree.nu_numode) > 0:
        mode = tree.nu_numode[0]
    elif hasattr(tree, "InteractionModeNuMode"):
        mode = getattr(tree, "InteractionModeNuMode")
        
    if mode == 0: return "Beam_QE"
    elif mode == 1: return "Beam_RES"
    elif mode == 2: return "Beam_DIS"
    elif mode == 10: return "Beam_MEC"
    elif mode == 3: return "Beam_COH"
    else: return "Beam_Other"

class SelectionEngine:
    def __init__(self, target_pot=1e21):
        self.target_pot = target_pot
        self.stages = []
        self.base_stage_name = "0. Generated (True Signal in FV)"
        
        self.categories = [
            "Signal", "Cosmic/Dirt", "Beam_QE", "Beam_RES", 
            "Beam_DIS", "Beam_MEC", "Beam_COH", "Beam_Other"
        ]

    def add_stage(self, name, cut_function):
        self.stages.append({
            "name": name,
            "func": cut_function
        })
        print(f"[*] Registered Selection Stage: {name}")

    def get_pot(self, filename):
        f = ROOT.TFile.Open(filename, "READ")
        if not f: return 1e18
        subrun_tree = f.Get("ana/subrunTree") or f.Get("subrunTree")
        pot = sum(getattr(sr, "pot", 0.0) for sr in subrun_tree) if subrun_tree else 0.0
        f.Close()
        return pot if pot > 0 else 1e18

    def run(self, signal_file, bkg_file):
        ROOT.gROOT.SetBatch(True)
        
        print("\n" + "="*80)
        print(" INITIALIZING SELECTION ENGINE")
        print("="*80)
        
        sig_pot = self.get_pot(signal_file)
        bkg_pot = self.get_pot(bkg_file)
        sig_scale = self.target_pot / sig_pot
        bkg_scale = self.target_pot / bkg_pot
        
        print(f" Target POT     : {self.target_pot:.2e}")
        print(f" Signal POT     : {sig_pot:.2e} (Scale: {sig_scale:.4f})")
        print(f" Background POT : {bkg_pot:.2e} (Scale: {bkg_scale:.4f})")
        print("-" * 80)

        # Initialise tracking dictionaries
        self.results = {}
        stage_names = [self.base_stage_name] + [s["name"] for s in self.stages]
        for name in stage_names:
            self.results[name] = {"raw_sig": 0, "raw_bkg": 0}
            for cat in self.categories:
                self.results[name][f"scaled_{cat}"] = 0.0

        # Signal file
        print("Processing Signal Sample (Extracting Signal + Rare Backgrounds)...")
        f_sig = ROOT.TFile.Open(signal_file, "READ")
        t_sig = f_sig.Get("ana/tree") or f_sig.Get("tree")
        
        for entry in range(t_sig.GetEntries()):
            t_sig.GetEntry(entry)
            
            cat = get_category(t_sig)
            if cat == "Signal": self.results[self.base_stage_name]["raw_sig"] += 1
            else: self.results[self.base_stage_name]["raw_bkg"] += 1
            self.results[self.base_stage_name][f"scaled_{cat}"] += sig_scale

            context = {} 
            for stage in self.stages:
                if stage["func"](t_sig, context):
                    if cat == "Signal": self.results[stage["name"]]["raw_sig"] += 1
                    else: self.results[stage["name"]]["raw_bkg"] += 1
                    self.results[stage["name"]][f"scaled_{cat}"] += sig_scale
                else: break
        f_sig.Close()

        # Background file
        print("Processing Background Sample (Extracting Bulk Backgrounds)...")
        f_bkg = ROOT.TFile.Open(bkg_file, "READ")
        t_bkg = f_bkg.Get("ana/tree") or f_bkg.Get("tree")
        
        for entry in range(t_bkg.GetEntries()):
            t_bkg.GetEntry(entry)
            
            cat = get_category(t_bkg)
            if cat == "Signal": continue # Vetoed! Avoid double counting
            
            self.results[self.base_stage_name]["raw_bkg"] += 1
            self.results[self.base_stage_name][f"scaled_{cat}"] += bkg_scale
                
            context = {}
            for stage in self.stages:
                if stage["func"](t_bkg, context):
                    self.results[stage["name"]]["raw_bkg"] += 1
                    self.results[stage["name"]][f"scaled_{cat}"] += bkg_scale
                else: break
        f_bkg.Close()
        
        self.print_table(sig_scale, bkg_scale)
        self.print_breakdown_table()
        self.plot_detailed_cutflow()

    def print_table(self, sig_scale, bkg_scale):
        print("\n" + "="*128)
        print(f" {'Selection Stage':<35} | {'Raw Sig':<8} | {'Raw Bkg':<8} | {'POT Sig':<10} | {'POT Bkg':<10} | {'Sig Eff %':<9} | {'Bkg Rej %':<9} | {'Pur %':<7}")
        print("-" * 128)
        
        base_sig_scaled = self.results[self.base_stage_name]["scaled_Signal"]
        base_bkg_scaled = sum(self.results[self.base_stage_name][f"scaled_{c}"] for c in self.categories if c != "Signal")
        
        stage_names = [self.base_stage_name] + [s["name"] for s in self.stages]
        
        for name in stage_names:
            s_raw = self.results[name]["raw_sig"]
            b_raw = self.results[name]["raw_bkg"]
            
            s_scaled = self.results[name]["scaled_Signal"]
            b_scaled = sum(self.results[name][f"scaled_{c}"] for c in self.categories if c != "Signal")
            
            sig_eff = (s_scaled / base_sig_scaled) * 100 if base_sig_scaled > 0 else 0
            bkg_rej = (1.0 - (b_scaled / base_bkg_scaled)) * 100 if base_bkg_scaled > 0 else 0
            pur = (s_scaled / (s_scaled + b_scaled)) * 100 if (s_scaled + b_scaled) > 0 else 0
            
            print(f" {name:<35} | {s_raw:<8} | {b_raw:<8} | {s_scaled:<10.1f} | {b_scaled:<10.1f} | {sig_eff:>7.1f} % | {bkg_rej:>7.1f} % | {pur:>5.1f} %")
            
        print("="*128 + "\n")

    def print_breakdown_table(self):
        print("="*128)
        print(" DETAILED BACKGROUND BREAKDOWN (Expected Events Scaled to POT)")
        print("="*128)
        print(f" {'Selection Stage':<35} | {'Cosmic/Dirt':<12} | {'Beam QE':<10} | {'Beam RES':<10} | {'Beam DIS':<10} | {'Beam MEC':<10} | {'Beam COH':<10} | {'Other':<10}")
        print("-" * 128)
        
        stage_names = [self.base_stage_name] + [s["name"] for s in self.stages]
        for name in stage_names:
            r = self.results[name]
            print(f" {name:<35} | {r['scaled_Cosmic/Dirt']:<12.1f} | {r['scaled_Beam_QE']:<10.1f} | {r['scaled_Beam_RES']:<10.1f} | {r['scaled_Beam_DIS']:<10.1f} | {r['scaled_Beam_MEC']:<10.1f} | {r['scaled_Beam_COH']:<10.1f} | {r['scaled_Beam_Other']:<10.1f}")
        print("="*128 + "\n")

    def plot_detailed_cutflow(self, output_tag="cutflow_detailed"):
        ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "Detailed Cut Flow", 1400, 800)
        ROOT.gPad.SetBottomMargin(0.3)
        ROOT.gPad.SetRightMargin(0.2)
        ROOT.gPad.SetLogy(1) 
        
        stage_names = [self.base_stage_name] + [s["name"] for s in self.stages]
        n_bins = len(stage_names)
        
        hs = ROOT.THStack("hs", f"Associated Hyperon Selection Cut-Flow (Scaled to {self.target_pot:.0e} POT);;Expected Events")
        
        colors = {
            "Cosmic/Dirt": ROOT.kGray+1,
            "Beam_Other": ROOT.kGray+2,
            "Beam_COH": ROOT.kRed+2,
            "Beam_QE": ROOT.kAzure+2,
            "Beam_MEC": ROOT.kMagenta+2,
            "Beam_RES": ROOT.kOrange+7,
            "Beam_DIS": ROOT.kGreen+2,
            "Signal": ROOT.kBlue
        }
        
        h_dict = {}
        for cat in colors.keys():
            h_dict[cat] = ROOT.TH1F(f"h_{cat}", f"{cat}", n_bins, 0, n_bins)
            h_dict[cat].SetFillColor(colors[cat])
            h_dict[cat].SetLineColor(ROOT.kBlack)
            
        for i, name in enumerate(stage_names):
            label = name.split(". ")[-1] if ". " in name else name
            for cat in colors.keys():
                h_dict[cat].GetXaxis().SetBinLabel(i+1, label)
                h_dict[cat].SetBinContent(i+1, self.results[name]["scaled_" + cat])
                
        # Stack from bottom to top
        stack_order = ["Cosmic/Dirt", "Beam_Other", "Beam_COH", "Beam_QE", "Beam_MEC", "Beam_RES", "Beam_DIS", "Signal"]
        for cat in stack_order:
            hs.Add(h_dict[cat])
            
        hs.Draw("BAR")
        hs.GetXaxis().SetLabelSize(0.045)
        
        leg = ROOT.TLegend(0.81, 0.3, 0.99, 0.9)
        leg.SetBorderSize(0)
        # Legend from top to bottom
        for cat in reversed(stack_order):
            leg.AddEntry(h_dict[cat], cat.replace("_", " "), "f")
        leg.Draw()
        
        c1.SaveAs(f"{output_tag}.pdf")
        print(f"[*] Saved detailed cut-flow visual to {output_tag}.pdf")


# ==============================================================================
# SELECTION CUT DEFINITIONS
# ==============================================================================

def cut_nuscore_slice(tree, context):
    best_idx, max_score = -1, -1.0
    for i, score in enumerate(tree.slice_NuScore):
        if score > max_score:
            max_score, best_idx = score, i
    if max_score > 0.5:
        context["primary_slice_idx"] = best_idx 
        return True
    return False

def cut_reco_fv(tree, context):
    idx = context.get("primary_slice_idx")
    if idx is None: return False
    
    vtx_x, vtx_y, vtx_z = tree.slice_VtxX[idx], tree.slice_VtxY[idx], tree.slice_VtxZ[idx]
    if vtx_x == -999.0: return False
        
    in_x = 5.0 < abs(vtx_x) < 180.0
    in_y = -180.0 < vtx_y < 180.0
    in_z = 20.0 < vtx_z < 470.0
    
    if in_x and in_y and in_z:
        context["primary_slice_id"] = tree.slice_ID[idx]
        return True
    return False

def cut_topology(tree, context):
    slice_id = context.get("primary_slice_id")
    if slice_id is None: return False
    
    n_tracks, n_showers = 0, 0
    t_slice_ids, t_scores = getattr(tree, "track_SliceID", []), getattr(tree, "track_TrackScore", [])
    
    for i, sid in enumerate(t_slice_ids):
        if sid == slice_id:
            if t_scores[i] > 0.5: n_tracks += 1
            elif 0.0 <= t_scores[i] <= 0.5: n_showers += 1
                
    return (n_tracks >= 4 and n_showers <= 2)

def cut_flash_match(tree, context):
    slice_idx = context.get("primary_slice_idx")
    if slice_idx is None: return False
        
    opt0_scores = getattr(tree, "slice_Opt0Score", [])
    if len(opt0_scores) > slice_idx and opt0_scores[slice_idx] > 320.0:
        return True
    return False

def cut_muon_candidate(tree, context):
    slice_idx, slice_id = context.get("primary_slice_idx"), context.get("primary_slice_id")
    if slice_idx is None or slice_id is None: return False
        
    vtx_x, vtx_y, vtx_z = tree.slice_VtxX[slice_idx], tree.slice_VtxY[slice_idx], tree.slice_VtxZ[slice_idx]
        
    best_muon_idx, max_length, max_vtx_dist = -1, -1.0, 5.0
    
    t_slice_ids, t_scores = getattr(tree, "track_SliceID", []), getattr(tree, "track_TrackScore", [])
    t_lengths = getattr(tree, "track_Length", [])
    t_start_x, t_start_y, t_start_z = getattr(tree, "track_StartX", []), getattr(tree, "track_StartY", []), getattr(tree, "track_StartZ", [])
    t_chi2_p, t_chi2_mu = getattr(tree, "track_Chi2Proton", []), getattr(tree, "track_Chi2Muon", [])
    
    for i, sid in enumerate(t_slice_ids):
        if sid == slice_id and t_scores[i] > 0.6:
            dist = math.sqrt((t_start_x[i] - vtx_x)**2 + (t_start_y[i] - vtx_y)**2 + (t_start_z[i] - vtx_z)**2)
            if dist < max_vtx_dist:
                if t_lengths[i] > 32.0 and t_chi2_mu[i] < 18.0 and t_chi2_p[i] > 87.0:
                    if t_lengths[i] > max_length:
                        max_length, best_muon_idx = t_lengths[i], i
                    
    if best_muon_idx != -1:
        context["muon_track_idx"] = best_muon_idx
        return True
    return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--signal", required=True, help="Filtered Hyperon ROOT file")
    parser.add_argument("-b", "--background", required=True, help="Inclusive Background ROOT file")
    parser.add_argument("-p", "--pot", type=float, default=1e21)
    args = parser.parse_args()
    
    engine = SelectionEngine(target_pot=args.pot)
    
    engine.add_stage("1. NuScore Cut", cut_nuscore_slice)
    engine.add_stage("2. Optical Flash Match", cut_flash_match)
    engine.add_stage("3. Reco Vertex in FV", cut_reco_fv)
    engine.add_stage("4. Tracks and showers topology cut", cut_topology)
    engine.add_stage("5. Muon Candidate Identification", cut_muon_candidate)
    
    engine.run(args.signal, args.background)