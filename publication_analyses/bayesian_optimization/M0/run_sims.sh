#!/bin/sh
../../../build/main -treeOut ./trees.trees -mcmcOut analysis.log -branchOut branches.log -simulateM0 -M0 -simulationOutput simulation.log -sequentialTuningSim -bayesOpt 40000 -bayesOptFreq 1000 -numSimulations 20 -burnInIter 10000 -numIter 10000 -tuneFreq 250 -sampleFreq 25
