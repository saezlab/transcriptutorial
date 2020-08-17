
## Questions related to CARNIVAL input 
- What is the number of Dorothea scores to feed CARNIVAL (through _measObj_)?
It feels that the smaller number you use, 
the better the results will be. However, this is not the general experience.
The solver needs to check a lot more possible solutions when having few inputs for the prior knowledge network(PKN).
This situation could lead to long processing times and big gap value in the ILP solution.

- How the maping of progeny scores (_weightObj_) is done to the "representative" genes?
The list of representative genes contain 33 entries covering the 14 pathways calculated with progeny.
This list has been curated (by Panuwat and Bence) to select which genes are the most upstream activators
of each progeny pathways and/or the actual targets of the drug perturbation that were used to build the original progeny matrix.
The idea is that those are the genes whose activities are actually responsible for the activity of the pathway. 
So in a way, the progeny scores that are estimated from the RESPONSIVE genes can represent a proxy of the activity of the top upstream genes. 
So the progeny pathway activities are mapped on those upstream genes and help constrain the carnival solution space.

## Questions related to CARNIVAL solutions
- How large are the "family of models" that you usually get in a regular run with CARNIVAL?
The ILP solver returns up to 100 models with replacement.
That is, the algorithm will replace similar solutions with different ones to give a wider range of posibilities.
The number of potential optimal solutions can be reduced when more information is provided by user - perturbations and pathway activities.

- Can I get very different networks after multiple CARNIVAL runs on the same data and parameters?
<<<<<<< HEAD
No, CARNIVAL is deterministic. It should lead to the same solutions everytime you run it.

- If I don't know the state ( active / inactive ) of the perturbed node, what can I do?
You can run inverse CARNIVAL to test both states with all the possible perturbed nodes,
or you can give the name of all the nodes you are interested on without value.
See the transcriptutorial for an example of this.

- can I use CARNIVAL with another species that is not human?
You can! As long as all the inputs you give are consistent.
=======
ILP algorithm is deterministic by nature but when the problem is big (many constraints and large input), solvers can use various optimization algorithms which can lead to non-deterministic results. For IBM cplex in parallel mode, several option exist: auto, deterministic and opportunistic (see official documentation [here](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/ParallelMode.html)

- What are NodeType values "S" and "T" in nodesAttributes in the results of CARNIVAL? 
"S" means species, or perturbed nodes. "T" means targets, for which the measured values were provided as input. Empty value "" means intermediate nodes for which the activity will be solved by ILP. 
>>>>>>> ea6f3ccd87ad9397a6f6035e6bf3b02c5e91cc8b
