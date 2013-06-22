#!/usr/bin/python2.5
"""
ab-initio folding of zinc finger protein 
"""

import os, rosetta

simulation_iter = 1000
frag_iter = 500000

def low_res_folding(pose, move_9, move_3, monte_carlo):
    for _ in range(frag_iter):
        move_9.apply(pose)
        monte_carlo.boltzmann(pose)
    
    for _ in range(frag_iter):
        move_3.apply(pose)
        monte_carlo.boltzmann(pose)
        
    monte_carlo.show_state()

def main():
    
    parent_path = "/Users/yanxia/Documents/Workspace/PyRosetta_Practice/"
    resource_path = parent_path + "resources"
    os.chdir(resource_path)
    
    rosetta.init()
    # initiate pose and two score functions
    pose = rosetta.Pose()
    rosetta.make_pose_from_sequence(pose, 
                                    "GSSGSSGTGVKPYGCSQCAKTFSLKSQLIVHQRSHTGVKPSGPSSG", 
                                    "centroid")
    
    fa_scorefxn = rosetta.create_score_function("standard")
    ct_scorefxn = rosetta.create_score_function("score3")
    kt_value = 1
    
    # initiate fragment set
    fragmentSet9 = rosetta.ConstantLengthFragSet(9)
    fragmentSet3 = rosetta.ConstantLengthFragSet(3)
    fragmentSet9.read_fragment_file("zf_9mer.txt")
    fragmentSet3.read_fragment_file("zf_3mer.txt")
    
    # set up movemap and Fragment Mover
    movemap = rosetta.MoveMap()
    movemap.set_bb(True)
    move_9mer = rosetta.ClassicFragmentMover(fragmentSet9, movemap)
    move_3mer = rosetta.ClassicFragmentMover(fragmentSet3, movemap)
    
    # Monte Carlo
    mc_low = rosetta.MonteCarlo(pose, ct_scorefxn, kt_value)
    
    #set up small and shear movers    
    n_moves = 5
    small_mover = rosetta.SmallMover(movemap, kt_value, n_moves)
    shear_mover = rosetta.ShearMover(movemap, kt_value, n_moves)
    
    #set up minimize mover
    min_mover = rosetta.MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(fa_scorefxn)
    min_mover.min_type("linmin")
    min_mover.tolerance(0.5)
    
    #set up sequence mover and repeat mover
    seq_mover = rosetta.SequenceMover()
    seq_mover.add_mover(small_mover) 
    seq_mover.add_mover(min_mover)
    seq_mover.add_mover(shear_mover)
    seq_mover.add_mover(min_mover)
    
    # folding
    # first low resolution
    
    #ct_switch = rosetta.SwitchResidueTypeSetMover("centroid")
    #ct_switch.apply(pose)
    low_res_folding(pose, move_9mer, move_3mer, mc_low)
    
    # high resolution
    fa_switch = rosetta.SwitchResidueTypeSetMover("fa_standard")
    fa_switch.apply(pose)
    mc_high = rosetta.MonteCarlo(pose, fa_scorefxn, kt_value)
    
    for i in range(5):
        print "before: ", fa_scorefxn(pose)
        max_angle = 25 - 5 * i
        small_mover.angle_max("H", max_angle)
        small_mover.angle_max("E", max_angle)
        small_mover.angle_max("S", max_angle)
        shear_mover.angle_max("H", max_angle)
        shear_mover.angle_max("E", max_angle)
        shear_mover.angle_max("S", max_angle)
        
        for _ in range(simulation_iter):
            seq_mover.apply(pose)
            mc_high.boltzmann(pose)
        
        print "after: ", fa_scorefxn(pose)
    
    result_path = parent_path + "results/"
    os.chdir(result_path)
    pose.dump_pdb("ara.pdb")
    print "Done!"

if __name__ == "__main__":
    main()
