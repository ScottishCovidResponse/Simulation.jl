# Model structure

*Description of the model as it currently stands.*

## Epidemiology

* We currently have an SEI3HRD compartmental model: Susceptible-Exposed-Asymptomatic-Presymptomatic-Symptomatic-Hospitalised-Recovered-Dead, as seen in Figure 1 below. Each compartment may be broken up into several age categories.

### Virus update loop
* Virus grows as the sum of several Binomial draws from each infection category, ``NV_{k, c} \sim \mathrm{Poisson}(f_c)`` where ``NV`` is new virus, ``k`` is the current grid square, ``c`` is the disease class, and ``f_c`` is the force of infection generation per disease class.
* The force of infection is dependent on the population size of that disease class (``P_{c,k}``), the growth parameter of that disease class (``g_c``): ``f_c = g_c * P_{c,k}``
* The virus decays similarly, according to a set probability ``d`` and the inverse match with environment at that location (``T_k^{-1}``): Decayed virus ``\sim \mathrm{Binomial}(V_k, d * T_k^{-1})``
* The newly generated virus is distributed 'locally' in space (i.e. grid square ``k \rightarrow j``) via a Gaussian kernel: ``D_c(k,j)*NV_{k,c}`` where ``D_c`` is a Gaussian kernel per infection class.
* An optional 'commuting' distribution of the virus is also possible, in addition to the local gaussian spread. Each can be weighted per compartment, and turned off/on through this parameter choice. In lockdown, all commuting movements are turned off.
* Currently, only the symptomatic and asymptomatic infectious classes can generate and spread virus, and ``f_c = 0`` for all others.
* Newly generated virus infects other susceptible groups either as an instantaneous force of infection (e.g. aerosol transmission) or through environmental transmission (e.g. through contact with a surface), both at separate rates. See below for more details.

### Disease class update loop
* Birth/death per class: Susceptibles are born into the population at a set probability per class. There is also a background probability of death from each disease class.
  * ``\mathrm{Births}_{k,c} \sim \mathrm{Binomial}(P_{k,c}, b_c)``
  * ``\mathrm{Deaths}_{k,c} \sim \mathrm{Binomial}(P_{k,c}, d_c)``
* Transitions per class: transitions between disease classes happen according to a transition matrix ``M``, which are calculated as moves into the disease class from other categories: ``in_{k, c} = \sum_{j}^{J} \mathrm{Binomial}(P_{k, j}, M_{c, j})`` where ``j \in{1,..,J}`` number of classes.
* ``M`` is constructed in advance and is altered for movement from Susceptible to Exposed categories by the amount of virus available in that location, ``V_k``, or the instantaneous force of infection, ``F_k``, which disappears at the end of each step. Each has its own associated parameter to weight the different types of infection; ``\beta_e`` and ``\beta_f``, respectively.
* You can optionally decide whether or not to weight the transmission by a separate parameter (set 0 - 1), ``w``, to decide whether aerosol or contact-based transmission is density or frequency dependent:
  * Environmental reservoir = ``\beta_f * V_k / N^w``
  * Force of infection = ``\beta_e * F_k / N^w``



![](SEI3HRD.svg)

*Figure 1: Current model structure.*

## Code Structure

* The virus update loop happens first, parallelised over disease category. The virus must move between locations, so each process must have full access to the entire space. This will parallelise better when we have more disease categories, like age.
* The disease class update loop happens second, parallelised over space. There is no movement between locations, but instead transitions between different categories.
* An overall matrix housing the abundances per grid square and per disease category is housed in the `EpiSystem`, along with information on environment and a lookup table of moves between different grid squares for each kernel.
* At every iteration of the update step, this matrix is updated in place to avoid additional memory allocation.
* Movement of the virus is now deterministic in the first stage of the virus update loop, and then drawn from a Poisson distribution based on this result in the second stage.
