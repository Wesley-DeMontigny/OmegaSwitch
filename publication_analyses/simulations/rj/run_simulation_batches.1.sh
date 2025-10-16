#!/bin/sh
../../../build/main -treeOut ./trees.1.trees -mcmcOut analysis.1.log -tipsOut tips.1.log -simulateRJ -RJ -simulationOutput simulation.1.log -numSimulations 25 -burnInIter 5000 -numIter 30000 -tuneFreq 100 -sampleFreq 50 &
sleep 1
../../../build/main -treeOut ./trees.2.trees -mcmcOut analysis.2.log -tipsOut tips.2.log -simulateRJ -RJ -simulationOutput simulation.2.log -numSimulations 25 -burnInIter 5000 -numIter 30000 -tuneFreq 100 -sampleFreq 50 &
