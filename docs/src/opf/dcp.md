# DC-OPF

See [`DCPPowerModel`](https://lanl-ansi.github.io/PowerModels.jl/stable/formulation-details/#PowerModels.DCPPowerModel) formulation in [`PowerModels.jl`](https://lanl-ansi.github.io/PowerModels.jl/stable/).

## Mathematical formulation

The primal problem reads

```math
\begin{align}
    \min_{\mathbf{p}^{\text{g}}, \mathbf{pf}, \boldsymbol{\theta}} \quad 
    & \label{eq:DCP:objective}
        \sum_{g \in \mathcal{G}} c_{g} \mathbf{p}^{\text{g}}_{g} + c^{0}_{g}
        \\
    \text{s.t.} \quad
    & \label{eq:DCP:slack_bus}
        \boldsymbol{\theta}_{\text{slack}} = 0
        &
        &&& [\lambda^{\text{slack}}]\\
    & \label{eq:DCP:kirchhoff}
        \sum_{g \in \mathcal{G}_{i}} \mathbf{p}^{\text{g}}_{g} 
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{pf}_{e}
        + \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{pf}_{e}
        = 
        \sum_{l \in \mathcal{L}_{i}} p^{d}_{l}
        + \sum_{s \in \mathcal{S}_{i}} g^{s}_{s}
        & \forall i \in \mathcal{N}
        &&& [\lambda^{\text{kcl}}]\\
    & \label{eq:DCP:ohm}
        \mathbf{pf}_{e}
        =
        b_{e} (\boldsymbol{\theta}_{t(e)} - \boldsymbol{\theta}_{s(e)})
        & \forall e \in \mathcal{E}
        &&& [\lambda^{\text{ohm}}]\\
    & \label{eq:DCP:va_diff}
        \underline{\Delta \boldsymbol{\theta}}_{e}
        \leq
        \boldsymbol{\theta}_{t(e)} - \boldsymbol{\theta}_{s(e)}
        \leq 
        \overline{\Delta \boldsymbol{\theta}}_{e}
        & \forall e \in \mathcal{E}
        &&& [\mu^{\Delta \boldsymbol{\theta}}]\\
    & \label{eq:DCP:thermal}
        -\overline{f}_{e} \leq \mathbf{pf}_{e} \leq \overline{f}_{e}
        & \forall e \in \mathcal{E}
        &&& [\mu^{f}]\\
    & \label{eq:DCP:pg_bounds}
        \underline{p}^{g}_{g} \leq \mathbf{p}^{\text{g}}_{g} \leq \overline{p}^{g}_{g}
        & \forall g \in \mathcal{G}
        &&& [\mu^{pg}]
\end{align}
```

## Nomenclature

### Primal variables

| Symbol | Data | Size | Description 
|:-------|:-----|:-----|:------------|
| ``\mathbf{p}^{\text{g}}`` | `pg` | ``G`` | Active power dispatch
| ``\boldsymbol{\theta}`` | `va` | ``N`` | Nodal voltage angle
| ``\mathbf{pf}`` | `pf` | ``E`` | Active power flow

### Dual variables

| Symbol | Data | Size | Associated constraint 
|:-------|:-----|:-----|:------------|
| ``\lambda^{\text{kcl}}`` | `lam_kirchoff` | ``N`` | Nodal power balance ``\eqref{eq:DCP:kirchhoff}``
| ``\lambda^{\text{ohm}}`` | `lam_ohm` | ``E`` | Ohm's law ``\eqref{eq:DCP:ohm}``
| ``\mu^{\Delta \boldsymbol{\theta}}`` | `mu_va_diff` | ``E`` | Angle difference limit ``\eqref{eq:DCP:va_diff}``
| ``\mu^{pf}`` | `mu_sm` | ``E`` | Thermal limit ``\eqref{eq:DCP:thermal}``
| ``\mu^{pg}`` | `mu_pg` | ``G`` | Generation min/max limits ``\eqref{eq:DCP:pg_bounds}``
| ``\lambda^{\text{slack}}`` | -- | ``1`` | Slack bus voltage angle ``\eqref{eq:DCP:slack_bus}``

Dual variable ``\lambda^{\text{slack}}`` is always zero at the optimum, hence it is not exported.