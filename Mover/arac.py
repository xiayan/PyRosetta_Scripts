#!/usr/bin/python2.5
"""
This program does ab initio folding
"""

from rosetta import *
import os.path

def main():
    """
    This program does ab inito folding
    """
    rosetta.init()
    pose = Pose("model.pdb")
    original_pose = Pose("model.pdb")
    
    #set up small and shear movers    
    kt_value = 1.0
    n_moves = 5
    movemap = MoveMap()
    # movemap.set_bb(False)
    movemap.set_bb_true_range(1, 15)
    movemap.set_chi(True)

    small_mover = SmallMover(movemap, kt_value, n_moves)
    shear_mover = ShearMover(movemap, kt_value, n_moves)
    
    #set up minimize mover
    scorefxn = create_score_function("standard")
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)
    min_mover.min_type("linmin")
    min_mover.tolerance(0.5)

    #set up monte-carlo mover
    mc_mover = MonteCarlo(pose, scorefxn, kt_value)

    #set up sequence mover and repeat mover
    seq_mover = SequenceMover()
    seq_mover.add_mover(small_mover) 
    seq_mover.add_mover(min_mover)
    seq_mover.add_mover(shear_mover)
    seq_mover.add_mover(min_mover)

    simulation_iter = 100

    for i in range(5):
        print "First pose:", scorefxn(pose)
        max_angle = 25 - 5 * i
        print "max angle: ", str(max_angle)
        small_mover.angle_max("H", max_angle)
        small_mover.angle_max("E", max_angle)
        small_mover.angle_max("S", max_angle)
        shear_mover.angle_max("H", max_angle)
        shear_mover.angle_max("E", max_angle)
        shear_mover.angle_max("S", max_angle)
        
        for i in range(simulation_iter):
            seq_mover.apply(pose)
            mc_mover.boltzmann(pose)
        
        mc_mover.boltzmann(pose)
        mc_mover.show_scores()
        mc_mover.show_counters()
        mc_mover.show_state()
        print "Last pose:", scorefxn(pose)
        
    print "ara original score: ", scorefxn(original_pose)
    print "ara final score: ", scorefxn(pose)
    mc_mover.recover_low(pose)
    print "ara low score: ", scorefxn(pose)

    original_pose.dump_pdb("ara_start.pdb")
    pose.dump_pdb("ara_finish.pdb")

if __name__ == "__main__": 
    main()
