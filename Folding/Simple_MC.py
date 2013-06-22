#!/usr/bin/python2.5
import argparse
import math, time
from rosetta import *
from random import gauss, randint, random, seed

simIteration = 10001

def perturb_pose(aPose, tempPose, residueTotal):
    tempPose.assign(aPose)
    testingResidue = randint(1, residueTotal - 1)
    tempPose.set_phi(testingResidue, gauss(tempPose.phi(testingResidue), 25.0))
    tempPose.set_psi(testingResidue, gauss(tempPose.psi(testingResidue), 25.0))
    tempPose.set_omega(testingResidue, gauss(tempPose.omega(testingResidue), 25.0))
    
def fragmentMove(aPose, tempPose, mover):
    tempPose.assign(aPose)
    mover.apply(tempPose)
    
def compareScore(currentScore, newScore):
    deltaE = newScore - currentScore
    if (deltaE < 0) :
        return True
    else:
        prob = random()
        newProb = math.e**(-deltaE)
        if (newProb > prob):
#            print "random, math", str(prob), str(newProb)
            return True
        else:
#            print "random, math", str(prob), str(newProb)
            return False

def setUpScoreFunction(scorefxn, runningMode):
    if runningMode:
        scorefxn = create_score_function("standard")
    else:
        scorefxn = create_score_function("score3")
    return scorefxn

def folding(pose, tempPose, sequence, runningMode, fragMode, index):
    residueTotal = pose.total_residue() + 1
    if not runningMode:
        switch1 = SwitchResidueTypeSetMover("centroid")
        switch1.apply(pose)
        switch2 = SwitchResidueTypeSetMover("fa_standard")
    
    #outfile = open("data.txt", "w")
    
    scorefxn = ScoreFunction()
    scorefxn = setUpScoreFunction(scorefxn, runningMode)
    
    lowScore = scorefxn(pose)
    acceptPose = False
    
    if fragMode == True:
        fragset = ConstantLengthFragSet(3)
        fragset.read_fragment_file("aat000_03_05.200_v1_3.txt")
        movemap = MoveMap()
        movemap.set_bb(True)
        mover_3mer = ClassicFragmentMover(fragset, movemap)
    
    for i in range(simIteration+1):
        if fragMode == False:
            perturb_pose(pose, tempPose, residueTotal)
        else: 
            fragmentMove(pose, tempPose, mover_3mer)
        
        newScore = scorefxn(tempPose)
        
        acceptPose = compareScore(lowScore, newScore)
        if acceptPose:
            pose.assign(tempPose)
            lowScore = newScore
            
        #outfile.write(str(i + 1) + "\t" + str(lowScore) + "\n")
        
        #if i == 0:
         #   if not runningMode:
          #      switch2.apply(tempPose)
           # tempPose.dump_pdb("first.pdb")
        #elif i == simIteration:
         #   if not runningMode:
          #      switch2.apply(tempPose)
           # tempPose.dump_pdb("last.pdb")
            
    
    if not runningMode:
        switch2.apply(tempPose)
    print index, " index"
    outputName = "decoy" + str(index) + ".pdb"
    pose.dump_pdb(outputName)
#    outfile.close()
    
def main():
    parser = argparse.ArgumentParser(description="Toggle between full-atom and centroid")
    parser.add_argument('-mode', dest='atomMode', default='fa_input', choices=['fa_input', 'centroid'])
    parser.add_argument('-s', dest='sequence', default='AAAAAAAAAA', type=str)
    parser.add_argument('-frag', dest='frag', default=False, action='store_true')
    parser.add_argument('-i', dest='iter', default=101, type=int)
    runningMode = True

    inputMode = parser.parse_args()
    if inputMode.atomMode == 'centroid':
        print "centroid mode"
        runningMode = False
    else:
        print "full atom mode"
    
    print "Fragment Mode: ", inputMode.frag
    
    sequence = inputMode.sequence
    
    sequence = sequence.upper()
    print sequence
    
    rosetta.init()
    seed()
    
    for index in range(inputMode.iter):
        startTime = time.time()
        print "Simulating: ", str(index+1)
        
        pose = Pose()
        tempPose = Pose()
        make_pose_from_sequence(pose, sequence, "fa_standard")
        for i in range(1, pose.total_residue() + 1):
            pose.set_omega(i, 180.0)
        
        folding(pose, tempPose, sequence, runningMode, inputMode.frag, index)
        print "total time: ", time.time() - startTime
    
if __name__ == "__main__" : main()
