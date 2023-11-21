make
./runCreatePythiaEvents -nev 10000 -pthat 120 -tune 14
./runSimpleJetAnalysis -hard PythiaEventsTune14PtHat120.pu14  -nev 10000
root JetToyHIResultSimpleJetAnalysis.root -l

