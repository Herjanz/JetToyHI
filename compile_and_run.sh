make
./runCreatePythiaEvents -nev 1000 -pthat 120 -tune 14
./runSimpleJetAnalysis -hard PythiaEventsTune14PtHat120.pu14  -nev 1000
root JetToyHIResultSimpleJetAnalysis.root -l
