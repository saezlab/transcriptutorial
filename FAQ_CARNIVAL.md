
- What is the number of Dorothea scores to feed CARNIVAL (through _measObj_)?
It feels that the smaller number you use, 
the better the results will be. However, this is not the general experience.
The solver needs to check a lot more possible solutions when having few inputs for the PKN.
This situation could lead to long processing times and big gaps.

- How large are the "family of models" that you usually get in a regular run with CARNIVAL?
The solver returns up to 100 models with replacement.
That is, the algorithm will replace similar solutions with different ones to give a wider range of posibilities.
The number of potential optimal solutions can be reduced when using perturbations and pathway activities,
and you can end up with a couple of solutions.

- How the maping of progeny scores (_weightObj_) is done to the "representative" genes?
The list of representative genes contain 33 entries covering the 14 pathways calculated with progeny.
This list has been curated (by Panuwat and Bence) to select which genes are the most upstream activators
of each progeny pathways and/or the actual targets of the drug perturbation that were used to build the original progeny matrix.
The idea is that those are the genes whose activities are actually responsible for the activity of the pathway. 
So in a way, the progeny scores that are estimated from the RESPONSIVE genes can represent a proxy of the activity of the top upstream genes. 
So the progeny pathway activities are mapped on those upstream genes and help constrain the carnival solution space.

- Can I get very different networks after multiple CARNIVAL runs on the same data and parameters?
No, CARNIVAL is deterministic. It should lead to the same solutions everytime you run it.

- If I don't know the state ( active / inactive ) of the perturbed node, what can I do?
You can run inverse CARNIVAL to test both states with all the possible perturbed nodes,
or you can give the name of all the nodes you are interested on without value.
See the transcriptutorial for an example of this.

- can I use CARNIVAL with another species that is not human?
You can! As long as all the inputs you give are consistent.
