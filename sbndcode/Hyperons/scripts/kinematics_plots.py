#!/usr/bin/env python3
import ROOT
import math
import argparse

def analyze_truth_physics(input_filename, output_tag="truth"):
    # Suppress pop-ups and stats boxes
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0) 
    
    input_file = ROOT.TFile.Open(input_filename, "READ")
    tree = input_file.Get("ana/tree") or input_file.Get("tree")
    if not tree:
        print("Error: Could not find tree in file.")
        return

    # Original Topological Limits
    h_flight_dist = ROOT.TH1F("h_flight_dist", "#Lambda Flight Distance;Distance [cm];Events", 50, 0, 30)
    h_flight_dist.SetLineColor(ROOT.kBlack)
    h_flight_dist.SetLineWidth(2)
    
    h_fsi_gen = ROOT.TH1F("h_fsi_gen", "#Lambda Decay Generation;Generation (1 = Primary);Events", 5, 0.5, 5.5)
    h_fsi_gen.SetFillColorAlpha(ROOT.kAzure+1, 0.5)
    h_fsi_gen.SetLineColor(ROOT.kAzure+3)
    h_fsi_gen.SetLineWidth(2)

    h2_hits_len_p = ROOT.TH2F("h2_hits_len_p", "Proton Visibility;True Track Length [cm];Expected Hit Count", 40, 0, 40, 40, 0, 200)
    h2_hits_len_pi = ROOT.TH2F("h2_hits_len_pi", "Pion Visibility;True Track Length [cm];Expected Hit Count", 40, 0, 40, 40, 0, 200)

    h_opening_angle = ROOT.TH1F("h_opening_angle", "Decay Opening Angle (p, #pi^{-});Angle [degrees];Events", 36, 0, 180)
    h_opening_angle.SetLineColor(ROOT.kGreen+2)
    h_opening_angle.SetLineWidth(2)

    # Kinematic phasespace
    h_nu_e = ROOT.TH1F("h_nu_e", "True Neutrino Energy;E_{#nu} [GeV];Events", 40, 0, 8.0)
    h_nu_e.SetLineColor(ROOT.kBlack)
    h_nu_e.SetLineWidth(2)

    h_k_ke = ROOT.TH1F("h_k_ke", "Primary K^{+} Kinetic Energy;Kinetic Energy [GeV];Events", 40, 0, 6)
    h_lambda_p = ROOT.TH1F("h_lambda_p", "Primary #Lambda Momentum;Momentum [GeV/c];Events", 40, 0, 6)
    h_k_ke.SetLineColor(ROOT.kBlue)
    h_lambda_p.SetLineColor(ROOT.kRed)
    h_k_ke.SetLineWidth(2)
    h_lambda_p.SetLineWidth(2)

    h_p_mom = ROOT.TH1F("h_p_mom", "Decay Proton Momentum;p [GeV/c];Events", 40, 0, 5.0)
    h_pi_mom = ROOT.TH1F("h_pi_mom", "Decay Pion Momentum;p [GeV/c];Events", 40, 0, 5.0)
    h_p_mom.SetLineColor(ROOT.kRed)
    h_pi_mom.SetLineColor(ROOT.kBlue)
    h_p_mom.SetLineWidth(2)
    h_pi_mom.SetLineWidth(2)

    h_pt = ROOT.TH1F("h_pt", "Transverse Momentum (p_{T, #Lambda} + p_{T, K^{+}});p_{T} [GeV/c];Events", 40, 0, 4.0)
    h_pt.SetLineColor(ROOT.kMagenta+2)
    h_pt.SetLineWidth(2)

    print(f"Processing {tree.GetEntries()} events...")


    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        
        # Truth tagging
        if getattr(tree, "IsAssocLambdaKPlusWithProtonPiMinus", 0) != 1 or getattr(tree, "IsInsideFV", 0) != 1:
            continue
            
        h_nu_e.Fill(tree.TrueNuE)
        h_fsi_gen.Fill(tree.LambdaDecayGen)

        pt_k, pt_lam = 0.0, 0.0
        lambda_track_id_gen = -1

        # Find the primary K+ and Lambda
        for i, pdg in enumerate(tree.particlePDG):

            if pdg == 321:
                h_k_ke.Fill(tree.particleKE[i])
                pt_k = math.sqrt(tree.particlePx[i]**2 + tree.particlePy[i]**2)
            elif pdg == 3122:
                p_mag = math.sqrt(tree.particlePx[i]**2 + tree.particlePy[i]**2 + tree.particlePz[i]**2)
                h_lambda_p.Fill(p_mag)
                pt_lam = math.sqrt(tree.particlePx[i]**2 + tree.particlePy[i]**2)
                lambda_track_id_gen = tree.particleTrackID[i]

        h_pt.Fill(pt_k + pt_lam)

        # G4 plots
        lambda_track_ids = [tree.geant_TrackID[i] for i, pdg in enumerate(tree.geant_PDG) if pdg == 3122]
        decaying_lambda_idx = -1
        proton_idx, pion_idx = -1, -1
        
        for l_id in lambda_track_ids:
            p_tmp, pi_tmp = -1, -1
            for i, m_id in enumerate(tree.geant_MotherID):
                if m_id == l_id:
                    if tree.geant_PDG[i] == 2212: p_tmp = i
                    if tree.geant_PDG[i] == -211: pi_tmp = i
            
            if p_tmp != -1 and pi_tmp != -1:
                decaying_lambda_idx = list(tree.geant_TrackID).index(l_id)
                proton_idx, pion_idx = p_tmp, pi_tmp
                break
        
        if decaying_lambda_idx != -1:
            # Flight Distance
            dx = tree.geant_endX[decaying_lambda_idx] - tree.geant_startX[decaying_lambda_idx]
            dy = tree.geant_endY[decaying_lambda_idx] - tree.geant_startY[decaying_lambda_idx]
            dz = tree.geant_endZ[decaying_lambda_idx] - tree.geant_startZ[decaying_lambda_idx]
            h_flight_dist.Fill(math.sqrt(dx**2 + dy**2 + dz**2))
            
            # Daughter Track Lengths
            p_dx = tree.geant_endX[proton_idx] - tree.geant_startX[proton_idx]
            p_dy = tree.geant_endY[proton_idx] - tree.geant_startY[proton_idx]
            p_dz = tree.geant_endZ[proton_idx] - tree.geant_startZ[proton_idx]
            p_len = math.sqrt(p_dx**2 + p_dy**2 + p_dz**2)
            
            pi_dx = tree.geant_endX[pion_idx] - tree.geant_startX[pion_idx]
            pi_dy = tree.geant_endY[pion_idx] - tree.geant_startY[pion_idx]
            pi_dz = tree.geant_endZ[pion_idx] - tree.geant_startZ[pion_idx]
            pi_len = math.sqrt(pi_dx**2 + pi_dy**2 + pi_dz**2)
            
            # Since true hits are broken currently, I will approx
            expected_hits_p = p_len / 0.3 if p_len > 0 else 0
            expected_hits_pi = pi_len / 0.3 if pi_len > 0 else 0

            h_p_mom.Fill(tree.geant_P[proton_idx])
            h_pi_mom.Fill(tree.geant_P[pion_idx])

            h2_hits_len_p.Fill(p_len, expected_hits_p)
            h2_hits_len_pi.Fill(pi_len, expected_hits_pi)
            
            # Opening Angle
            if p_len > 0 and pi_len > 0:
                dot_product = (p_dx * pi_dx + p_dy * pi_dy + p_dz * pi_dz)
                cos_theta = max(-1.0, min(1.0, dot_product / (p_len * pi_len))) 
                h_opening_angle.Fill(math.acos(cos_theta) * (180.0 / math.pi))

    # Save shizzle
    c1 = ROOT.TCanvas("c1", "Kinematics Canvas", 800, 600)
    
    # True Neutrino Energy
    h_nu_e.Draw("HIST")
    c1.SaveAs(f"{output_tag}_neutrino_energy.pdf")

    # Transverse Momentum Sum
    h_pt.Draw("HIST")
    c1.SaveAs(f"{output_tag}_pt_sum.pdf")

    # Daughter Momenta (Overlay)
    h_p_mom.Draw("HIST")
    h_pi_mom.Draw("HIST SAME")
    leg_mom = ROOT.TLegend(0.6, 0.7, 0.85, 0.85)
    leg_mom.AddEntry(h_p_mom, "Proton", "l")
    leg_mom.AddEntry(h_pi_mom, "Pion", "l")
    leg_mom.Draw()
    c1.SaveAs(f"{output_tag}_daughter_momenta.pdf")

    # Flight Distance
    h_flight_dist.Draw("HIST")
    c1.SaveAs(f"{output_tag}_flight_distance.pdf")
    
    # FSI Generation
    h_fsi_gen.Draw("HIST")
    c1.SaveAs(f"{output_tag}_fsi_generation.pdf")
    
    # Opening Angle
    h_opening_angle.Draw("HIST")
    c1.SaveAs(f"{output_tag}_opening_angle.pdf")
    
    # Phase Space (Overlay)
    h_k_ke.Draw("HIST")
    h_lambda_p.Draw("HIST SAME")
    leg_phase = ROOT.TLegend(0.5, 0.7, 0.85, 0.85)
    leg_phase.AddEntry(h_lambda_p, "#Lambda Momentum", "l")
    leg_phase.AddEntry(h_k_ke, "K^{+} Kinetic Energy", "l")
    leg_phase.Draw()
    c1.SaveAs(f"{output_tag}_phase_space.pdf")
    
    # Visibility 2D Plots
    h2_hits_len_p.Draw("COLZ")
    c1.SaveAs(f"{output_tag}_proton_visibility.pdf")
    
    h2_hits_len_pi.Draw("COLZ")
    c1.SaveAs(f"{output_tag}_pion_visibility.pdf")
    
    input_file.Close()
    print("All kinematics plots generated successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output_tag", default="truth")
    args = parser.parse_args()
    analyze_truth_physics(args.input, args.output_tag)