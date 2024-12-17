# DC-OPF

## Mathematical formulation

The primal problem reads

```math
\begin{align}
    \min_{\PG, \PF, \VA} &\quad
        \sum_{i \in \NODES} \sum_{j \in \GENERATORS_{i}} c_j \PG_j \label{model:dcopf:obj} \\
    \text{s.t.} \quad
    & \sum_{j\in\GENERATORS_i}\PG_j - \sum_{e \in \mathcal{E}_{i}}  \PF_{e} + \sum_{e \in \mathcal{E}^R_{i}} \PF_{e}
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
    & \underline{\PG_i} \leq \PG_i \leq \overline{\PG_i}
        & \forall i &\in \GENERATORS
    \label{model:dcopf:pgbound} \\
    & {-}\overline{S_{e}} \leq  \PF_{e} \leq \overline{S_{e}}
        & \forall e &\in \EDGES
    \label{model:dcopf:thrmbound:from}
\end{align}
```

## Nomenclature

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
