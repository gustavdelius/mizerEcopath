This package provides functions that help to create a mizer model based on the
parameters of an Ecopath model.

An Ecopath model requires the user to specify for each group $i$ included in the
model the following parameters:

* Biomass $B_i$
* Consumption $Q_i$
* Production $P_i$
* Fishery catch $C_i$
* The diet composition $D_{ij}$ that gives the proportion of the diet of group
  $i$ that is made up of group $j$

One of the first three parameters above can be replaced by the "Ecotrophic
Efficiency" $EE$, because Ecopath can then deduce the missing parameter, see
below.

We want to use this information to create a mizer model that in its steady state
has the same values for the above parameters. The challenge is that mizer,
which is a size-structured model, is parametrized differently. Instead of
specifying properties at the group level, we specify parameters at the level of
the individual and these parameters are size-dependent. In this document we
will discuss how we deduce mizer parameters from the Ecopath parameters.

An Ecopath model usually includes a large number of ecosystem components
(called groups), only some of which it would be appropriate to model explicitly
in mizer. We will therefore need to select a subset of the Ecopath groups to
include in the mizer model. Everything else we will treat as external forcing
on the mizer model. So the groups that we do not include in the mizer model
will still be sources of consumption and of mortality for the groups that we do
include.

In cases where a species in the Ecopath model is split into several stanzas,
the rates for the stanzas can be added together to get the total rate for the
species. Mizer then decides itself how the rate is distributed over the mizer
size classes.

The Ecopath parameters give us no direct information on the size structure of
the populations. Therefore we supplement the Ecopath parameters with information
about the size structure of the catches. This information is usually readily
available from fishery-based observations.

We furthermore make some allometric assumptions, i.e., we assume that in the
steady state the following parameters for an individual of weight $w$ scale as a
power of that weight:

* The consumption rate
* The metabolic respiration rate
* The mortality rate
* The rate at which a mature individual invests in reproduction

Of course, when we turn on the dynamics of the mizer model (which corresponds to
going from Ecopath to Ecosim) and leave the steady state, the rates will change
in response to changes in the abundances of prey and predators.

## Definition of the Ecopath parameters

To start, we need to understand the definition of the Ecopath parameters and
how to calculate their values in a mizer model.

### Biomass
The biomass parameter $B$ in Ecopath is the total biomass of all the individuals
in the group in the model. In mizer, the total biomass of a species is given by
integrating over all sizes:
$$B_i = \int_0^\infty N_i(w) w dw$$
where $N_i(w)$ is the number density of the species $i$ at size $w$.

### Consumption
The consumption parameter $Q$ in Ecopath is the total rate at which all the
individuals taken together consume food. In mizer, This is calculated as 
$$Q_i = \int E_{e.i}(w) (1 - f_i(w)) N_i(w) dw$$ 
where $E_{e.i}(w)$ is the encounter rate of an individual of species $i$ and
weight $w$ (calculated with getEncounter), $f_i(w)$ is the feeding level of
an individual of species $i$ and weight $w$ (calculated with
`getFeedingLevel()`) and $N_i(w)$ is the number density of species $i$ at
weight $w$.
