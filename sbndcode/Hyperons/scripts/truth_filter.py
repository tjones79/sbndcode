#!/usr/bin/env python3
import ROOT
import sys
import argparse
import array

FV_Y_MIN, FV_Y_MAX = -180.0, 180.0
FV_Z_MIN, FV_Z_MAX =  20.0, 470.0

PDG_NAMES = {
    22: "gamma", 111: "pi0", 211: "pi+", -211: "pi-",
    13: "mu-", -13: "mu+", 2112: "n", 2212: "p", 
    3122: "Lambda", 321: "K+", -321: "K-",
    1000010020: "Deuteron", 1000020040: "Alpha"
}

def get_particle_name(pdg):
    return PDG_NAMES.get(pdg, f"PDG({pdg})")

def is_inside_true_fv(vx, vy, vz):
    in_x = (-190.0 < vx < -5.0) or (5.0 < vx < 190.0)
    in_y = (FV_Y_MIN < vy < FV_Y_MAX)
    in_z = (FV_Z_MIN < vz < FV_Z_MAX)
    return in_x and in_y and in_z

def process_and_clone(input_filename, output_filename, target_pot, verbose=False):
    input_file = ROOT.TFile.Open(input_filename, "READ")
    if not input_file or input_file.IsZombie():
        print(f"Error: Could not open input file {input_filename}")
        return

    input_tree = input_file.Get("ana/tree")
    input_subrun_tree = input_file.Get("ana/subrunTree")
    
    if not input_tree:
        print("Error: Could not find main TTree 'ana/tree' in input file.")
        return

    print(f"Opened input file: {input_filename}")
    
    total_sim_pot = 0.0
    if input_subrun_tree:
        for sr in input_subrun_tree:
            total_sim_pot += sr.pot

    if total_sim_pot == 0:
        print("Warning: Total simulated POT is 0. Scaling factor will default to 1.0.")
        scale_factor = 1.0
    else:
        scale_factor = target_pot / total_sim_pot

    print(f"Total Simulated POT : {total_sim_pot:.3e}")
    print(f"Target Exposure POT : {target_pot:.3e}")
    print(f"POT Scale Factor    : {scale_factor:.5f}")
    print(f"Cloning structure for {input_tree.GetEntries()} events...\n")

    output_file = ROOT.TFile.Open(output_filename, "RECREATE")
    output_dir = output_file.mkdir("ana")
    output_dir.cd()
    output_tree = input_tree.CloneTree(0)

    # Output array branches
    is_assoc_lambda_k_plus = array.array('i', [0])
    is_assoc_lambda_k_plus_with_proton_pi_minus = array.array('i', [0])
    is_inside_fv_arr = array.array('i', [0])
    n_mu_minus_arr = array.array('i', [0])
    n_k_plus_arr = array.array('i', [0]) 
    n_lambdas_arr = array.array('i', [0]) 
    lambda_decay_gen_arr = array.array('i', [0])  

    output_tree.Branch("IsAssocLambdaKPlus", is_assoc_lambda_k_plus, "IsAssocLambdaKPlus/I")
    output_tree.Branch("IsAssocLambdaKPlusWithProtonPiMinus", is_assoc_lambda_k_plus_with_proton_pi_minus, "IsAssocLambdaKPlusWithProtonPiMinus/I")
    output_tree.Branch("IsInsideFV", is_inside_fv_arr, "IsInsideFV/I")
    output_tree.Branch("n_mu_minus", n_mu_minus_arr, "n_mu_minus/I")
    output_tree.Branch("n_k_plus", n_k_plus_arr, "n_k_plus/I")
    output_tree.Branch("n_lambdas", n_lambdas_arr, "n_lambdas/I")
    output_tree.Branch("LambdaDecayGen", lambda_decay_gen_arr, "LambdaDecayGen/I")

    for entry in range(input_tree.GetEntries()):
        if entry > 0 and entry % 50000 == 0:
            print(f"{entry}/{input_tree.GetEntries()} events analysed")

        input_tree.GetEntry(entry)
        
        is_assoc_lambda_k_plus[0] = 0
        is_assoc_lambda_k_plus_with_proton_pi_minus[0] = 0
        is_inside_fv_arr[0] = 0
        n_mu_minus_arr[0] = 0
        n_k_plus_arr[0] = 0
        n_lambdas_arr[0] = 0
        lambda_decay_gen_arr[0] = 0
        event = input_tree

        # 1. Group stable GENIE particles by their interaction index to prevent pile-up mixing
        interaction_dict = {}
        for pdg, status, idx in zip(event.particlePDG, event.status_code, event.particle_interaction_index):
            if status == 1: # Only look at final-state particles
                if idx not in interaction_dict:
                    interaction_dict[idx] = []
                interaction_dict[idx].append(pdg)

        in_fv = False

        # 2. Check if ANY single interaction has the exact topology
        for idx, pdgs in interaction_dict.items():
            # Safety check: ensure index exists and parent was a muon neutrino
            if idx < len(event.nu_pdg) and event.nu_pdg[idx] == 14:
                
                c_mu = pdgs.count(13)
                c_k = pdgs.count(321)
                c_lambda = pdgs.count(3122)
                c_sigma0 = pdgs.count(3212)
                
                if c_mu >= 1 and c_k >= 1 and c_lambda >= 1 and c_sigma0 == 0:
                    is_assoc_lambda_k_plus[0] = 1
                    
                    # Verify if THIS specific signal interaction occurred inside the FV
                    if idx < len(event.TrueVtxX):
                        in_fv = is_inside_true_fv(event.TrueVtxX[idx], event.TrueVtxY[idx], event.TrueVtxZ[idx])
                        if in_fv:
                            is_inside_fv_arr[0] = 1

                    # Record the multiplicities of this specific signal interaction
                    n_mu_minus_arr[0] = c_mu
                    n_k_plus_arr[0] = c_k
                    n_lambdas_arr[0] = c_lambda
                    
                    # Once we found the signal interaction, we can stop searching this event
                    break

        # 3. If we have exactly our signal topology inside the FV, trace the Lambda decay
        if (is_assoc_lambda_k_plus[0] == 1 and in_fv):

            # Track lambda generations
            lambda_gens = {}
            for i, pdg in enumerate(event.geant_PDG):
                if pdg == 3122 and (event.geant_MotherID[i] == 10000000 or event.geant_MotherID[i] == 0):
                    lambda_gens[event.geant_TrackID[i]] = 1

            # find re-scattered child Lambdas 
            made_progress = True
            while made_progress:
                made_progress = False
                for i, pdg in enumerate(event.geant_PDG):
                    t_id = event.geant_TrackID[i]
                    m_id = event.geant_MotherID[i]
                    if pdg == 3122 and (m_id in lambda_gens) and (t_id not in lambda_gens):
                        lambda_gens[t_id] = lambda_gens[m_id] + 1
                        made_progress = True

            signal_decay_found = False
            decayed_track_id = None
            decayed_gen = 0
            
            for lambda_id, gen in sorted(lambda_gens.items(), key=lambda x: x[1]):
                daughter_pdgs = [event.geant_PDG[i] for i, m_id in enumerate(event.geant_MotherID) if m_id == lambda_id]
                
                # Check for charged decay node: Lambda -> p + pi-
                # Eliminate background nuclear splits by enforcing a clean 2-body decay topology
                if 2212 in daughter_pdgs and -211 in daughter_pdgs and len(daughter_pdgs) <= 3:
                    signal_decay_found = True
                    decayed_track_id = lambda_id
                    decayed_gen = gen
                    break 

            if signal_decay_found:
                is_assoc_lambda_k_plus_with_proton_pi_minus[0] = 1
                lambda_decay_gen_arr[0] = decayed_gen

                if verbose:
                    print(f"\n[SIGNAL MATCH] Run {event.run}, Event {event.event}")
                    print(f"  --> Signal decay captured at Lambda Generation {decayed_gen}")
                    print("  --> Complete Geant4 History Chain:")
                    for l_id, gen in sorted(lambda_gens.items(), key=lambda x: x[1]):
                        daughters = [event.geant_PDG[i] for i, m_id in enumerate(event.geant_MotherID) if m_id == l_id]
                        daughter_names = [get_particle_name(pdg) for pdg in daughters]
                        
                        marker = "==>" if l_id == decayed_track_id else "-->"
                        print(f"    {marker} [Gen {gen}] Lambda Track ID {l_id}")
                        print(f"         Daughters: {', '.join(daughter_names) if daughter_names else 'None (Escaped TPC)'}")
        output_tree.Fill()

    if input_subrun_tree:
        output_subrun_tree = input_subrun_tree.CloneTree()
        output_subrun_tree.Write()

    # Final metrics
    topo_count = output_tree.GetEntries("IsAssocLambdaKPlus == 1")
    topo_fv_count = output_tree.GetEntries("IsAssocLambdaKPlus == 1 && IsInsideFV >= 1")
    topo_with_decay_count = output_tree.GetEntries("IsAssocLambdaKPlusWithProtonPiMinus == 1 && IsInsideFV >= 1")

    output_tree.Write() 
    output_file.Close()
    input_file.Close()

    # yield summary
    print("\n" + "="*70)
    print(f"   INCLUSIVE SIGNAL TRUTH (Scaled to {target_pot:.2e} POT)")
    print("="*70)
    print("Target Topology: nu_mu + Ar -> mu- + K+ + Lambda (-> p + pi-)")
    print("-" * 70)
    print(f"Events matching base topology  : {topo_count} | (Scaled: {topo_count * scale_factor:.2f})")
    print(f"Events with topology in FV     : {topo_fv_count} | (Scaled: {topo_fv_count * scale_factor:.2f})")
    print(f"Events with topology + decay    : {topo_with_decay_count} | (Scaled: {topo_with_decay_count * scale_factor:.2f})")
    print("="*70)
    print(f"Output written to: {output_filename}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, required=True)
    parser.add_argument("-o", "--output_file", type=str, required=True)
    parser.add_argument("-p", "--target_pot", type=float, default=1e21)
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()
    process_and_clone(args.input_file, args.output_file, args.target_pot, verbose=args.verbose)