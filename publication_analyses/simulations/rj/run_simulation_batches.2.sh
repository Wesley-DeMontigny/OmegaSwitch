#!/bin/sh
../../../build/main -treeOut ./trees.3.trees -mcmcOut analysis.3.log -tipsOut tips.3.log -simulateRJ -RJ -simulationOutput simulation.3.log -numSimulations 25 -burnInIter 5000 -numIter 30000 -tuneFreq 100 -sampleFreq 50 &
sleep 1
../../../build/main -treeOut ./trees.4.trees -mcmcOut analysis.4.log -tipsOut tips.4.log -simulateRJ -RJ -simulationOutput simulation.4.log -numSimulations 25 -burnInIter 5000 -numIter 30000 -tuneFreq 100 -sampleFreq 50 &
