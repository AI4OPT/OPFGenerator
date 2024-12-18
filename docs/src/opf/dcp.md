# DC-OPF

## Mathematical Formulation

The DCOPF model considered in OPFGenerator is presented below.

```math
\begin{align}
    \min_{\PG, \PF, \VA} &\quad
        \sum_{i \in \NODES} \sum_{j \in \GENERATORS_{i}} c_j \PG_j + c_0 \label{model:dcopf:obj} \\
    \text{s.t.} \quad
    & \sum_{j\in\GENERATORS_i}\PG_j - \sum_{e \in \mathcal{E}^{+}_{i}}  \PF_{e} + \sum_{e \in \mathcal{E}^{-}_{i}} \PF_{e}
    = \sum_{j\in\LOADS_i}\PD_j + \GS_i 
        & \forall i &\in \NODES
    \label{model:dcopf:kirchhoff} \\
    & {-}b_{e}(\VA_{i} - \VA_{j}) - \PF_{e} = 0
        & \forall e = (i, j) &\in \EDGES
    \label{model:dcopf:ohm} \\
& \dvamin_{e} \leq \VA_i - \VA_j \leq \dvamax_{e}
        & \forall e = (i, j) &\in \EDGES
    \label{model:dcopf:angledifference} \\
    & \VA_\text{ref} = 0 \label{model:dcopf:slackbus} \\
    & \pgmin_{i} \leq \PG_i \leq \pgmax_{i}
        & \forall i &\in \GENERATORS
    \label{model:dcopf:pgbound} \\
    & {-}\overline{S}_{e} \leq  \PF_{e} \leq \overline{S}_{e}
        & \forall e &\in \EDGES
    \label{model:dcopf:thrmbound:from}
\end{align}
```


### Variables

* ``\PG \in \mathbb{R}^{G}``: active power dispatch
* ``\PF \in \mathbb{R}^{E}``: active power flow "from"
* ``\VA \in \mathbb{R}^{N}``: nodal voltage angle

### Objective

The objective function ``\eqref{model:dcopf:obj}`` minimizes the cost of active power generation.

!!! todo
    OPFGenerator currently supports only linear cost functions.
    Support for quadratic functions is planned for a later stage; please open an issue if 
    you would like to request this feature.

### Constraints

* ``\eqref{model:dcopf:kirchhoff}``: Kirchhoff's current law.
* ``\eqref{model:dcopf:ohm}``: Ohm's law expressing power flows as a function of nodal voltages.
* ``\eqref{model:dcopf:thrmbound:from}``: Thermal limits.
* ``\eqref{model:dcopf:angledifference}``: Voltage angle deviation constraints.
* ``\eqref{model:dcopf:slackbus}``: this constraint fixes the voltage angle of the reference (slack) bus to zero.
* ``\eqref{model:dcopf:pgbound}``: Active generation limits.



## Data Format

### Primal variables

| Symbol | Data | Size | Description 
|:-------|:-----|:-----|:------------|
| ``\PG`` | `pg` | ``G`` | Active power dispatch
| ``\VA`` | `va` | ``N`` | Nodal voltage angle
| ``\PF`` | `pf` | ``E`` | Active power flow

### Dual variables

| Associated constraint                             | Data         | Size  |
|:--------------------------------------------------|:-------------|:------|
| ``\eqref{model:dcopf:slackbus}``                  | `slack_bus`  | ``1`` |
| ``\eqref{model:dcopf:kirchhoff}``                 | `kcl_p`      | ``N`` |
| ``\eqref{model:dcopf:ohm}``                       | `ohm`        | ``E`` |
| ``\eqref{model:dcopf:angledifference}``           | `va_diff`    | ``E`` |
| ``\eqref{model:dcopf:pgbound}`` (lower)           | `pg_min`     | ``G`` |
| ``\eqref{model:dcopf:pgbound}`` (upper)           | `pg_max`     | ``G`` |
| ``\eqref{model:dcopf:thrmbound:from}`` (lower)    | `pf_min`     | ``E`` |
| ``\eqref{model:dcopf:thrmbound:from}`` (upper)    | `pf_max`     | ``E`` |
