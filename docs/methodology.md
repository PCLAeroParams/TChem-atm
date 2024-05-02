TChem solves a system of ordinary differential equations (ODEs) in order to time advance the volumetric mixing ratio (vmr,$\eta_k$) of gas species, $k$.
$$
\frac{\dif{} \eta_k}{\dif{} t}=\dot{\omega}_k
$$

The net production rate of species $k$, $\dot{\omega}_k$, or the right hand side of the vmr equation is computed using:

$$
  \dot{\omega}_k=\sum_{i=1}^{N_{\text{react}}}\nu_{ki}q_i,\quad \nu_{ki}=\nu''_{ki}-\nu'_{ki},
$$

where $q_i$ is the rate of progress of reaction $i$, $N_{\text{react}}$ is the number of reactions, $\nu''_{ki}$ and $\nu'_{ki}$ are the stoichiometric coefficients of species $k$ in reaction $i$ for the reactant and product sides of the reaction, respectively.


Finally, the rate of progress of reaction $i$ is computed as

\begin{equation}\label{eq:rate_of_progress}
  q_i={k_f}_i\prod_{j=1}^{N_{\text{spec}}}\eta_j^{\nu'_{ji}},
\end{equation}

where $N_{\text{spec}}$ is the number of species, ${k_f}_i$ is the reaction constant of reaction $i$.
Note that in \eee{}`s CAMPP solver only forward reaction calculations are employed, and the single reaction constant values depend on the type of reaction.
In TChem, we have implemented separate reaction types for Troe, Arrhenius, and JPL-Troe types.
