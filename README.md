# Fitness cost associated with cell phenotypic switching drives population diversification dynamics and controllability

-----

Creator & Maintainer: _Lucas Henrion_ -- Contact at lucas.henrion@uliege.be if there are issues with the code.

Contributors: see CONTRIBUTING.md

Creation date: 2023-07-12

-----

To run [simulations](./simulation.py), [FlowStoCKS](./FlowStoCKS_main.py) takes the [process](./ProcessParam.py) and the [biological](./BioModel.py) parameters as input and returns a time scatter plot representing the induction level of the cells over time as well as a matrix containing the raw results of cells sizes and fluorescences.
Additionaly [KITopologyAnalysis](./KITopologyAnalysis.py) was constructed to analyse the influence of the inhibition strength associated to switching (KI value) for different environmental regimes (Chemostat and Segregostat). To run these simulations, the process conditions (dilution rate, glucose feed concentration ...) were set based on the experimental ones. For the Segregostat simulations, the pulsing time was extracted from the summary file of the *S. cerevisiae Glc3* Segregostat experiment.

FlowStoCKS requires the following packages to run:

1. numpy
2. pandas
3. random
4. matplotlib
5. itertools

## Publication history

This code was used to generate the results of the following publication:

"Fitness cost associated with cell phenotypic switching drives population diversification dynamics and controllability", Henrion et al., 2023, [pre-print](https://doi.org/10.1101/2023.04.06.535654) (version used for the pre-print: [v1.0.1](../../tree/v1.0.1), version used for the final text: [v1.0.2](../../tree/v1.0.2)) -- _LH and VV  are supported by the FRS-FNRS (Fond National pour la Recherche Scientifique, Belgium) through a FRIA PhD grant. MD is supported by a FNRS PhD grant in the context of an Era-Net Aquatic Pollutant project (ARENA). JAM is supported by a postdoctoral grant in the context of an Era-Cobiotech project (Contibio)._