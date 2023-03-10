# Spherical_TA
Routine used in Rafaela Filippozzi's thesis

# Folders and Files

- Algorithms_and_Routines: In this folder are the algorithms and Routines used in the experiments.
    * `SPHERICALTAPLUSHEU.m`: Triangle Algorithm for spherical problem with option to have heuristic
    * `RoutineEmpiricalInvestigation.m`: Routine used to generate the results for Empirical investigation of iteration complexity for the Spherical CHMP with p in convS
    * `Routine_Artificial_Instances_pinconvS.m`: Routine used to generate the results for the Spherical CHMP with p in convS
    * `Routine_Artificial_Instances_pnotinconvS.m`: Routine used to generate the results for the Spherical CHMP with p notin convS
    * `solverheuristic.m`:  Heuristic for solve min (Stilx)^Ts_j; e^Tx=1 x\geq 0 \norm{Stilx}^2 \leq \norm{p_k}^2
    Auxiliary functions used by Routines:
    * `Random_cvx.m` : Randomly generate data points in convex hull using vertices
    * `Random_pts.m`: Randomly generate data points 
    * `generateArandom.m`: Generates A randomly according to a uniform distribution on the unit ball of Rm

- Experiment_data_Routine_Empirical_Investigation:Data included in section 4.4.1
    
- Experiment_data_Routine_Artificial_Instances: Data included in section 4.4.2

