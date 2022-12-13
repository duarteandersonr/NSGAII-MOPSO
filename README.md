# NSGAII-MOPSO
Codes and input files used in the article "Post-processing Improvements in Multi-objective Optimization of General Single-server Finite Queuing Networks"

The file gapsoMG1k.for is the source code in Fortran language (Fortran77 compiler) for the execution of NSGA-II (elitist non-dominated sorting genetic algorithm) and MOPSO (multi-objective particle swarm optimization algorithm).

The input files have the following encoding, Me for Merge topology, Se for Serie topology, Sp for Split topology, and Mx for Mixed topology. In addition, NqXX is for the number of queues in the queueing network, for queue number, LaX, arrival rate value, and SXX for the coefficient of variation squared value. Example: SeNq5La5S05.txt, queueing network with series topology, five queues in the queueing network, arrival rate equal to 5, squared coefficient of variation equal to 0.5.

The execution procedure is carried out through the command prompt in the directory with the gapsoMG1k.exe file and the input files through the following command line: gapsoMG1k.exe < input_file >> output_file

The gapsoMG1k.exe file is the executable file generated after Fortran compilation.
