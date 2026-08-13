#!/usr/bin/env python3
import ROOT
import argparse
from collections import Counter, defaultdict

def get_particle_name(pdg):
    mapping = {
        13: "#mu^{-}", -13: "#mu^{+}",
        14: "#nu_{#mu}", -14: "#bar{#nu}_{#mu}",
        11: "e^{-}", -11: "e^{+}",
        12: "#nu_{e}", -12: "#bar{#nu}_{e}",
        3122: "#Lambda^{0}",
        321: "K^{+}", -321: "K^{-}",
        311: "K^{0}", 130: "K^{0}_{L}", 310: "K^{0}_{S}",
        3222: "#Sigma^{+}", 3212: "#Sigma^{0}", 3112: "#Sigma^{-}",
        2212: "p", 2112: "n",
        211: "#pi^{+}", -211: "#pi^{-}", 111: "#pi^{0}",
        22: "#gamma"
    }
    return mapping.get(pdg, str(pdg))

def sort_particles(pdg_list):
    def sort_key(p):
        if abs(p) < 20: return (0, abs(p))       
        if abs(p) > 3000: return (1, abs(p))     
        return (2, abs(p))                       
    return sorted(pdg_list, key=sort_key)

def run_mc_study(input_filename, output_tag="mc_study"):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("ana/tree") or input_file.Get("tree")
    subrun_tree = input_file.Get("ana/subrunTree") or input_file.Get("subrunTree")
    
    if not tree:
        print("Error: Could not find event tree in file.")
        return

    total_pot = 0.0
    if subrun_tree:
        for entry in range(subrun_tree.GetEntries()):
            subrun_tree.GetEntry(entry)
            total_pot += getattr(subrun_tree, "pot", 0.0)
    
    if total_pot == 0.0:
        print("Warning: Total POT is 0. Using a placeholder for scaling (1e18).")
        total_pot = 1e18

    pot_scale_factor = 1e21 / total_pot

    pdgs_to_process = [14, -14]
    labels = {14: "nu", -14: "antinu"}
    titles = {14: "#nu_{#mu}", -14: "#bar{#nu}_{#mu}"}
    
    mode_map = {
        0:  ("QE", ROOT.kAzure+2),
        1:  ("RES", ROOT.kOrange+7),
        2:  ("DIS", ROOT.kGreen+2),
        10: ("MEC", ROOT.kMagenta+2),
        3:  ("COH", ROOT.kRed+2),
        -1: ("Unknown", ROOT.kGray+2)
    }
    stack_order = [3, 10, 2, 1, 0, -1] 

    h_inc = {p: {} for p in pdgs_to_process}
    h_sig = {p: {} for p in pdgs_to_process}
    h_flight_dist = {p: {} for p in pdgs_to_process}
    h_k_mom = {p: {} for p in pdgs_to_process}
    h_open_angle = {p: {} for p in pdgs_to_process}
    h_p_mom = {p: {} for p in pdgs_to_process}
    h_pi_mom = {p: {} for p in pdgs_to_process}

    for p in pdgs_to_process:
        for m, (name, color) in mode_map.items():
            # Energy Histograms
            h_inc[p][m] = ROOT.TH1F(f"h_inc_{p}_{m}", f"{titles[p]} {name} Inclusive;True Neutrino Energy [GeV];Events", 30, 0, 6)
            h_inc[p][m].SetFillColorAlpha(color, 0.9)
            
            h_sig[p][m] = ROOT.TH1F(f"h_sig_{p}_{m}", f"{titles[p]} {name} Signal;True Neutrino Energy [GeV];Events", 30, 0, 6)
            h_sig[p][m].SetFillColorAlpha(color, 0.9)

            # Flight Distance Histograms
            h_flight_dist[p][m] = ROOT.TH1F(f"h_flight_dist_{p}_{m}", f"{titles[p]} {name} True #Lambda^{{0}} Flight Distance;Flight Distance [cm];Events", 50, 0, 30)
            h_flight_dist[p][m].SetFillColorAlpha(color, 0.9)

            # Kinemati
            h_k_mom[p][m] = ROOT.TH1F(f"h_k_mom_{p}_{m}", f"{titles[p]} {name} True Primary K^{{+}} Momentum;Momentum [GeV/c];Events", 40, 0, 2)
            h_k_mom[p][m].SetFillColorAlpha(color, 0.9)
            
            h_open_angle[p][m] = ROOT.TH1F(f"h_open_angle_{p}_{m}", f"{titles[p]} {name} #Lambda^{{0}} Daughter Opening Angle (p, #pi^{{-}});3D Opening Angle [Degrees];Events", 36, 0, 180)
            h_open_angle[p][m].SetFillColorAlpha(color, 0.9)
            
            h_p_mom[p][m] = ROOT.TH1F(f"h_p_mom_{p}_{m}", f"{titles[p]} {name} True Daughter Proton Momentum;Momentum [GeV/c];Events", 40, 0, 2)
            h_p_mom[p][m].SetFillColorAlpha(color, 0.9)
            
            h_pi_mom[p][m] = ROOT.TH1F(f"h_pi_mom_{p}_{m}", f"{titles[p]} {name} True Daughter Pion (#pi^{{-}}) Momentum;Momentum [GeV/c];Events", 40, 0, 1)
            h_pi_mom[p][m].SetFillColorAlpha(color, 0.9)

    topology_counter = {p: defaultdict(Counter) for p in pdgs_to_process}
    decay_counter = {p: {m: Counter() for m in mode_map} for p in pdgs_to_process}

    print("Running MC study...")

    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        
        # Grab vectors
        nu_pdgs = list(tree.nu_pdg)
        nu_modes = list(tree.nu_numode)
        nu_energies = list(tree.nu_energy)
        
        pdgs = list(tree.particlePDG)
        status = list(tree.status_code)
        interaction_indices = list(tree.particle_interaction_index) 
        
        for i in range(len(nu_pdgs)):
            p_pdg = nu_pdgs[i]
            if p_pdg in pdgs_to_process:
                m = nu_modes[i] if nu_modes[i] in mode_map else -1
                e = nu_energies[i]
                if e >= 0:
                    h_inc[p_pdg][m].Fill(e)

        hyperon_interaction_indices = set()
        for p, s, idx in zip(pdgs, status, interaction_indices):
            if s == 1 and abs(p) in [3122, 3222, 3212, 3112]:
                hyperon_interaction_indices.add(idx)

        for target_idx in hyperon_interaction_indices:
            if target_idx >= len(nu_pdgs): continue
            
            nu_pdg = nu_pdgs[target_idx]
            if nu_pdg not in pdgs_to_process: continue
            
            mode = nu_modes[target_idx]
            if mode not in mode_map: mode = -1
            energy = nu_energies[target_idx]
            
            if energy >= 0:
                h_sig[nu_pdg][mode].Fill(energy)

            # Isolate topology
            interaction_particles = [p for p, s, idx in zip(pdgs, status, interaction_indices) if s == 1 and idx == target_idx]
            key_particles = [p for p in interaction_particles if abs(p) in [11, 12, 13, 14, 321, 311, 130, 310, 3122, 3222, 3212, 3112]]
            unique_pdgs = list(set(key_particles))
            sorted_keys = sort_particles(unique_pdgs)
            topology_str = "".join([get_particle_name(p) for p in sorted_keys])
            if not topology_str: topology_str = "Other"
            topology_counter[nu_pdg][topology_str][mode] += 1

            # tree.particleP is the momentum of the GENIE primary particle
            if hasattr(tree, "particleP"):
                genie_momenta = list(tree.particleP)
                for p_id, s, p_mom, idx in zip(pdgs, status, genie_momenta, interaction_indices):
                    if s == 1 and p_id == 321 and idx == target_idx:
                        h_k_mom[nu_pdg][mode].Fill(p_mom)

            # Lambda analysis: flight distance, decay daughters, kinematics
            if hasattr(tree, "geant_PDG"):
                geant_pdgs = list(tree.geant_PDG)
                geant_track_ids = list(tree.geant_TrackID)
                geant_mother_ids = list(tree.geant_MotherID)
                geant_sx = list(tree.geant_startX)
                geant_sy = list(tree.geant_startY)
                geant_sz = list(tree.geant_startZ)
                geant_ex = list(tree.geant_endX)
                geant_ey = list(tree.geant_endY)
                geant_ez = list(tree.geant_endZ)
                geant_px = list(tree.geant_Px)
                geant_py = list(tree.geant_Py)
                geant_pz = list(tree.geant_Pz)
                
                for i, gpdg in enumerate(geant_pdgs):
                    # Is it a Lambda0?
                    if abs(gpdg) == 3122:
                        tid = geant_track_ids[i]
                        
                        # Flight
                        dx = geant_ex[i] - geant_sx[i]
                        dy = geant_ey[i] - geant_sy[i]
                        dz = geant_ez[i] - geant_sz[i]
                        flight_dist = (dx**2 + dy**2 + dz**2)**0.5
                        h_flight_dist[nu_pdg][mode].Fill(flight_dist)
                        
                        # Find da duaghter
                        daughters = []
                        p_idx = -1
                        pi_idx = -1
                        
                        for j, mid in enumerate(geant_mother_ids):
                            if mid == tid:
                                daughters.append(geant_pdgs[j])
                                if geant_pdgs[j] == 2212: p_idx = j
                                if geant_pdgs[j] == -211: pi_idx = j
                        
                        # If this Lambda decayed into p + pi-
                        if p_idx != -1 and pi_idx != -1:
                            v_p = ROOT.TVector3(geant_px[p_idx], geant_py[p_idx], geant_pz[p_idx])
                            v_pi = ROOT.TVector3(geant_px[pi_idx], geant_py[pi_idx], geant_pz[pi_idx])
                            
                            angle = v_p.Angle(v_pi) * 180.0 / ROOT.TMath.Pi()
                            h_open_angle[nu_pdg][mode].Fill(angle)
                            
                            h_p_mom[nu_pdg][mode].Fill(v_p.Mag())
                            h_pi_mom[nu_pdg][mode].Fill(v_pi.Mag())

                        if daughters:
                            # Ignore generic argon nuclei fragments if Geant4 interacted weirdly
                            daughters = [d for d in daughters if abs(d) < 10000]
                            daughters.sort(key=lambda x: abs(x), reverse=True) 
                            decay_str = " + ".join([get_particle_name(d) for d in daughters])
                            decay_counter[nu_pdg][mode][decay_str] += 1
                        else:
                            decay_counter[nu_pdg][mode]["No Daughters (Escaped/Absorbed)"] += 1


    c1 = ROOT.TCanvas("c1", "MC Study", 1400, 600)

    for p in pdgs_to_process:
        tag = labels[p]
        title_tag = titles[p]
        
        for is_scaled in [False, True]:
            if is_scaled:
                for m in mode_map.keys():
                    h_inc[p][m].Scale(pot_scale_factor)
                    h_sig[p][m].Scale(pot_scale_factor)
                    h_flight_dist[p][m].Scale(pot_scale_factor)
                suffix = "_POT_scaled"
                y_title = "Expected Events / 1e21 POT"
            else:
                suffix = ""
                y_title = "Events"
            
            c1.Clear()
            ROOT.gPad.SetLeftMargin(0.12)
            hs_inclusive = ROOT.THStack(f"hs_inc_{p}{suffix}", f"{title_tag} True Energy by Interaction Mode;True Neutrino Energy [GeV];{y_title}")
            for m in stack_order:
                hs_inclusive.Add(h_inc[p][m])
            hs_inclusive.Draw("HIST")
            
            if hs_inclusive.GetStack():
                h_inclusive_total = hs_inclusive.GetStack().Last().Clone(f"h_inc_total_{p}{suffix}")
            else:
                h_inclusive_total = ROOT.TH1F(f"h_inc_total_{p}{suffix}", "", 30, 0, 6)

            leg1 = ROOT.TLegend(0.7, 0.55, 0.88, 0.88)
            leg1.SetBorderSize(0)
            for m in reversed(stack_order):
                leg1.AddEntry(h_inc[p][m], mode_map[m][0], "f")
            leg1.Draw()
            c1.SaveAs(f"{output_tag}_{tag}_stacked_Energy_Inclusive{suffix}.pdf")

            c1.Clear()
            ROOT.gPad.SetLeftMargin(0.12)
            hs_hyperons = ROOT.THStack(f"hs_sig_{p}{suffix}", f"{title_tag} Hyperon Signal by Interaction Mode;True Neutrino Energy [GeV];{y_title}")
            for m in stack_order:
                hs_hyperons.Add(h_sig[p][m])
            hs_hyperons.Draw("HIST")
            
            leg2 = ROOT.TLegend(0.7, 0.55, 0.88, 0.88)
            leg2.SetBorderSize(0)
            for m in reversed(stack_order):
                leg2.AddEntry(h_sig[p][m], mode_map[m][0], "f")
            leg2.Draw()
            c1.SaveAs(f"{output_tag}_{tag}_stacked_Energy_Hyperons{suffix}.pdf")

            c1.Clear()
            c1.SetLogy()
            h_inclusive_total.SetTitle(f"{title_tag} Signal vs Inclusive Overlay;True Neutrino Energy [GeV];{y_title}")
            h_inclusive_total.SetFillStyle(0)
            h_inclusive_total.SetLineColor(ROOT.kBlack)
            
            if h_inclusive_total.GetMaximum() > 0:
                h_inclusive_total.SetMaximum(h_inclusive_total.GetMaximum() * 5)
                h_inclusive_total.SetMinimum(0.5)

            h_inclusive_total.Draw("HIST")
            hs_hyperons.Draw("HIST SAME")

            leg3 = ROOT.TLegend(0.7, 0.5, 0.88, 0.88)
            leg3.SetBorderSize(0)
            for m in reversed(stack_order):
                leg3.AddEntry(h_sig[p][m], mode_map[m][0], "f")
            leg3.AddEntry(h_inclusive_total, "Inclusive Total", "l")
            leg3.Draw()
            c1.SaveAs(f"{output_tag}_{tag}_stacked_Energy_Hyperons_overlayed_inclusive{suffix}.pdf")
            c1.SetLogy(0)

            top_topologies = []
            for topo_name, mode_counts in topology_counter[p].items():
                total_count = sum(mode_counts.values())
                top_topologies.append((topo_name, total_count, mode_counts))
            top_topologies.sort(key=lambda x: x[1], reverse=True)
            top_topologies = top_topologies[:10]

            c1.Clear()
            ROOT.gPad.SetBottomMargin(0.15)
            ROOT.gPad.SetLeftMargin(0.12)
            
            n_bins = max(len(top_topologies), 1)
            h_top = ROOT.TH1F(f"h_top_{p}{suffix}", f"{title_tag} Final State Topology Distribution;;{y_title}", n_bins, 0, n_bins)
            h_top.SetFillColorAlpha(ROOT.kAzure+1, 0.9)
            
            for i, (topo_name, total_count, _) in enumerate(top_topologies):
                h_top.GetXaxis().SetBinLabel(i+1, topo_name)
                h_top.SetBinContent(i+1, total_count * pot_scale_factor if is_scaled else total_count)
                
            h_top.GetXaxis().SetLabelSize(0.06)
            h_top.SetBarWidth(0.8); h_top.SetBarOffset(0.1)
            h_top.Draw("BAR")
            c1.SaveAs(f"{output_tag}_{tag}_topology_count{suffix}.pdf")
            
            def draw_stacked(h_dict, title_str, filename):
                c1.Clear()
                ROOT.gPad.SetLeftMargin(0.12)
                ROOT.gPad.SetBottomMargin(0.12)
                hs = ROOT.THStack("hs", title_str)
                for m in stack_order:
                    if is_scaled: h_dict[m].Scale(pot_scale_factor)
                    hs.Add(h_dict[m])
                hs.Draw("HIST")
                
                leg = ROOT.TLegend(0.7, 0.55, 0.88, 0.88)
                leg.SetBorderSize(0)
                for m in reversed(stack_order):
                    leg.AddEntry(h_dict[m], mode_map[m][0], "f")
                leg.Draw()
                c1.SaveAs(f"{output_tag}_{tag}_{filename}{suffix}.pdf")
            
            draw_stacked(h_k_mom[p], f"{title_tag} Primary K^{{+}} Momentum;Momentum [GeV/c];{y_title}", "k_momentum")
            draw_stacked(h_p_mom[p], f"{title_tag} #Lambda^{{0}} Daughter Proton Momentum;Momentum [GeV/c];{y_title}", "proton_momentum")
            draw_stacked(h_pi_mom[p], f"{title_tag} #Lambda^{{0}} Daughter Pion Momentum;Momentum [GeV/c];{y_title}", "pion_momentum")
            draw_stacked(h_open_angle[p], f"{title_tag} #Lambda^{{0}} Daughter Opening Angle (p, #pi^{{-}});Opening Angle [Degrees];{y_title}", "daughter_opening_angle")

        c1.Clear()
        c1.SetLogy(0)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetBottomMargin(0.12)
        
        hs_fd = ROOT.THStack(f"hs_fd_{p}", f"{title_tag} True #Lambda^{{0}} Flight Distance by Mode;Flight Distance [cm];Events")
        for m in stack_order:
            hs_fd.Add(h_flight_dist[p][m])
        hs_fd.Draw("HIST")
        
        line = ROOT.TLine(7.89, 0, 7.89, hs_fd.GetMaximum() if hs_fd.GetMaximum() > 0 else 1)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()
        
        leg_fd = ROOT.TLegend(0.6, 0.55, 0.88, 0.88)
        leg_fd.SetBorderSize(0)
        for m in reversed(stack_order):
            leg_fd.AddEntry(h_flight_dist[p][m], mode_map[m][0], "f")
        leg_fd.AddEntry(line, "c#tau (7.89 cm)", "l")
        leg_fd.Draw()
        c1.SaveAs(f"{output_tag}_{tag}_lambda_flight_distance.pdf")

        c1.Clear()
        ROOT.gPad.SetBottomMargin(0.2)
        
        # Aggregate decays for this neutrino flavor
        flavor_decays = Counter()
        for m in mode_map.keys():
            flavor_decays.update(decay_counter[p][m])
            
        top_decays = flavor_decays.most_common(5)
        n_bins = max(len(top_decays), 1)
        h_decay = ROOT.TH1F(f"h_decay_{p}", f"{title_tag} True #Lambda^{{0}} Decay Channels;;Events", n_bins, 0, n_bins)
        h_decay.SetFillColorAlpha(ROOT.kMagenta-4, 0.8)
        
        for i, (name, count) in enumerate(top_decays):
            h_decay.GetXaxis().SetBinLabel(i+1, name)
            h_decay.SetBinContent(i+1, count)
            
        h_decay.GetXaxis().SetLabelSize(0.06)
        h_decay.SetBarWidth(0.7)
        h_decay.SetBarOffset(0.15)
        h_decay.Draw("BAR")
        c1.SaveAs(f"{output_tag}_{tag}_lambda_decay_modes.pdf")

    print("\n" + "="*85)
    print(" SBND ASSOCIATED HYPERON FEASIBILITY STUDY")
    print("="*85)
    print(f" Total POT Processed : {total_pot:.2e}")
    print(f" Scale to SBND Goal  : x {pot_scale_factor:.2f} (Assuming 1e21 POT)")
    
    mode_print_order = [0, 1, 2, 3, 10, -1] # QE, RES, DIS, COH, MEC, nonsense

    for p in pdgs_to_process:
        header_title = "NEUTRINO (nu_mu) INTERACTIONS" if p == 14 else "ANTINEUTRINO (anti-nu_mu) INTERACTIONS"
        print("\n" + "="*85)
        print(f" {header_title}")
        print("="*85)
        
        print(f" {'Primary Topology':<22} | {'Expected/1e21 POT':<17} | {'Mode Breakdown (Raw MC)'}")
        print("-" * 85)
        
        top_topos = []
        for topo_name, mode_counts in topology_counter[p].items():
            top_topos.append((topo_name, sum(mode_counts.values()), mode_counts))
        top_topos.sort(key=lambda x: x[1], reverse=True)
        
        for topo_name, total_count, mode_counts in top_topos[:10]:
            clean_name = topo_name.replace("^{0}","").replace("^{-}","-").replace("^{+}","+").replace("#","")
            expected_events = int(total_count * pot_scale_factor)
            
            breakdown_list = []
            for m in mode_print_order:
                if mode_counts[m] > 0:
                    breakdown_list.append(f"{mode_map[m][0]}: {mode_counts[m]}")
            
            breakdown_str = ", ".join(breakdown_list)
            print(f" {clean_name:<22} | {expected_events:<17} | {breakdown_str}")
    
    print("\n" + "-" * 85)
    print(" LAMBDA0 DECAY BRANCHING RATIOS BY MODE (GEANT4 TRUTH)")
    print("-" * 85)
    
    for p in pdgs_to_process:
        header = "NEUTRINO (nu_mu)" if p == 14 else "ANTINEUTRINO (anti-nu_mu)"
        print(f"\n >>> {header} <<<")
        
        flavor_decays = Counter()
        for m in mode_map.keys():
            flavor_decays.update(decay_counter[p][m])
            
        tot_lambdas = sum(flavor_decays.values())
        if tot_lambdas == 0:
            print("    No Lambdas found.")
            continue
            
        print(f"  TOTAL LAMBDAS: {tot_lambdas}")
        for decay_str, count in flavor_decays.most_common(3):
            clean = decay_str.replace("^{0}","").replace("^{-}","-").replace("^{+}","+").replace("#","")
            print(f"    {clean:<25} : {(count/tot_lambdas)*100:>5.1f}%  ({count} events)")
            
        print("\n  Breakdown by Top Producing Modes:")
        modes_by_volume = sorted(mode_map.keys(), key=lambda m: sum(decay_counter[p][m].values()), reverse=True)
        
        for m in modes_by_volume[:3]: # Print top 3 modes
            m_lambdas = sum(decay_counter[p][m].values())
            if m_lambdas > 0:
                print(f"    -> {mode_map[m][0]} (Produced {m_lambdas} Lambdas)")
                for decay_str, count in decay_counter[p][m].most_common(2):
                    clean = decay_str.replace("^{0}","").replace("^{-}","-").replace("^{+}","+").replace("#","")
                    print(f"         {clean:<23} : {(count/m_lambdas)*100:>5.1f}%")

    print("\n" + "="*85 + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output_tag", default="mc_study")
    args = parser.parse_args()
    run_mc_study(args.input, args.output_tag)