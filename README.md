# generalised_smoluchowski

Numerical solution of the Smoluchowski polymer condensation equations with an additional sink term to accound for ring formation. The explicit Euler scheme is used for the integration of the population equations. 

==================== input and output files =======================================

input files: average_length.txt. This file is generated after the topology reconstruction (see MariaPanoukidou/topology_reconstruction repository).

output files: none

output variables: topological number kappa, a dimensionless number related to the ratio of ring rate formation/ linear rate formation. The ring and linear rates account for the linear and ring polymer chains forming during a condesation process. 

===================== technical details ============================================

To run: In MATLAB version 2020a (or 2019) run the smol_fitting.m file. The input parameters are Nframes = number of lines in the average_length.txt files vol = volume used during the simulation (with LAMMPS or any other software). 

==================== information on the codes ======================================

smol_fitting.m: The main script where the least squares fitting (LSQ) is performed on the average_length calculated after a simulation where polymer condensation was happening. The topological parameter kappa is exported.

Obj_smoluchowski.m: The objective function used for the least squares fitting. This function is custom made and calls the exEuler_smoluchowski.m in every LSQ iteration. 

exEuler_smoluchowski.m: Explicit Euler scheme for solving the ordinary differential equations of the generalised Smoluchowski equation. 
