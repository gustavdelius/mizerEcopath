### **Appendix: Mathematical Formulation of the Simulation Method**

This appendix details the mathematical procedure used to generate a simulated sample from a size-structured population model that is directly comparable to a length-stratified empirical dataset collected across multiple survey dates. The method convolves a single-cohort simulation with an annual spawning distribution and an age-to-ring mapping function to predict the expected age-length structure on any given day.

#### **1. Model Components**

Let the state of the fish population be described by the number of individuals in discrete size classes $s \in \{s_1, ..., s_{n_s}\}$ at a continuous age $a$.

**1.1. Single-Cohort Dynamics (Impulse Response)**

The user's linear, time-invariant, size-structured model provides the temporal evolution of a single cohort of fish born at the same instant. We define the output of this simulation as the **impulse response matrix**, $\mathbf{G}$, where the element $G(s, a)$ is the number of fish of size class $s$ that have survived to age $a$ from an initial pulse of $N_0$ recruits.

**1.2. Annual Spawning Distribution**

The relative intensity of spawning throughout the year is modeled by a von Mises probability density function, $S(\tau)$, where $\tau \in [0, 1]$ represents the fractional time of year. The function is defined as:

$$
S(\tau) = \frac{e^{\kappa \cos(2\pi(\tau - \mu))}}{2\pi I_0(\kappa)}
$$

where $\mu$ is the mean time of year for spawning, $\kappa$ is the concentration parameter, and $I_0(\kappa)$ is the modified Bessel function of the first kind of order zero, which serves as a normalization constant.

**1.3. Age-to-Ring Mapping**

The number of observed otolith rings, $K$, is a deterministic function of a fish's true age $a$ and the date of the survey, $t_{survey}$. This mapping, $\mathcal{K}(a, t_{survey})$, is defined as:

$$
\mathcal{K}(a, t_{survey}) = \sum_{y=Y_{birth}+1}^{Y_{survey}-1} \mathbb{I}(t_{ring,y} - t_{birth} \ge a_{min})
$$

where:

* $t_{birth} = t_{survey} - a$ is the birth date.

* $Y_{birth}$ and $Y_{survey}$ are the years of birth and survey, respectively.

* $t_{ring,y}$ is the date of ring formation in year $y$.

* $a_{min}$ is the minimum age required for the first ring to form.

* $\mathbb{I}(\cdot)$ is the indicator function, which is 1 if the condition is true and 0 otherwise.

#### **2. Population and Observation Synthesis**

**2.1. Population Structure Convolution**

The expected number of fish of size $s$ and true age $a$ in the population on a given survey date, $t_{survey}$, is found by weighting the impulse response by the spawning intensity at the corresponding birth date. Let $\tau(t)$ be the function that returns the fractional time of year for any date $t$. The population structure, $\mathbf{N}_{pop}$, is given by:

$$
N_{pop}(s, a | t_{survey}) = G(s, a) \cdot S(\tau(t_{survey} - a))
$$

This convolution integrates the fate of all cohorts born throughout the year into a single snapshot of the population.

**2.2. Observation Convolution**

The expected number of fish of size $s$ that would be observed with $K$ rings, $\mathbf{N}_{model}$, is obtained by aggregating the true population structure according to the age-to-ring mapping function:

$$
N_{model}(s, K | t_{survey}) = \int_0^{a_{max}} N_{pop}(s, a | t_{survey}) \cdot \mathbb{I}(\mathcal{K}(a, t_{survey}) = K) \, da
$$

In the discrete-time implementation, this integral becomes a sum over all age steps $a_i$ where $\mathcal{K}(a_i, t_{survey}) = K$.

#### **3. Simulated Sampling for Comparison**

To compare the model output with the length-stratified data, we first convert the model's expected counts into a conditional age-length key.

**3.1. Conditional Probabilities**

For a given survey date $t_{survey}$, the probability of a fish in size class $s$ having $K$ rings is:

$$
P(K | s, t_{survey}) = \frac{N_{model}(s, K | t_{survey})}{\sum_{K'} N_{model}(s, K' | t_{survey})}
$$

**3.2. Multinomial Sampling**

Let the empirical data for a survey on date $t_{survey,j}$ consist of a set of observations, where $n_{s,j}$ is the number of fish sampled from size class $s$. To generate a comparable simulated sample, we draw from a multinomial distribution for each size class:

$$
(K_{s,1}, K_{s,2}, ..., K_{s,n_K}) \sim \text{Multinomial}(n_{s,j}, \mathbf{p}_{s,j})
$$

where $\mathbf{p}_{s,j}$ is the vector of probabilities $[P(K=0|s, t_{survey,j}), P(K=1|s, t_{survey,j}), ...]$. This process is repeated for each size class $s$ and for each unique survey date $j$. The resulting simulated counts are then combined to form a final simulated data frame that has the exact same structure and sampling effort as the empirical data.
            
