# AC-OPF


## Mathematical Formulation

The ACOPF model considered in OPFGenerator is presented below.

```math
\begin{align}
    \min \quad 
    & \label{eq:ACOPF:objective}
        \sum_{g \in \mathcal{G}} c_{g} \mathbf{p}^{\text{g}}_{g} + c^{0}_{g}\\
    \text{s.t.} \quad
    & \label{eq:ACOPF:slack}
        \boldsymbol{\theta}_{i_{s}} = 0\\
    & \label{eq:ACOPF:kcl_p}
        \sum_{g \in \mathcal{G}_{i}} \mathbf{p}^{\text{g}}_{g}
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{p}^{\text{f}}_{e}
        - \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{p}^{\text{t}}_{e}
        - g^{s}_{i} \mathbf{v}_{i}^{2}
        =
        \sum_{l \in \mathcal{L}_{i}} p^{d}_{l}
        && \forall i \in \mathcal{N}
        && [\lambda^{p}]\\
    & \label{eq:ACOPF:kcl_q}
        \sum_{g \in \mathcal{G}_{i}} \mathbf{q}^{\text{g}}_{g}
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{q}^{\text{f}}_{e}
        - \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{q}^{\text{t}}_{e}
        + b^{s}_{i} \mathbf{v}_{i}^{2}
        =
        \sum_{l \in \mathcal{L}_{i}} q^{d}_{l}
        && \forall i \in \mathcal{N}
        && [\lambda^{q}]\\
    % Ohm's law
    & \label{eq:ACOPF:ohm_pf}
        g^{ff}_{e} \mathbf{v}_{i}^{2}
        + g^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        + b^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{p}^{\text{f}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{pf}]\\
    & \label{eq:ACOPF:ohm_qf}
        -b^{ff}_{e} \mathbf{v}_{i}^{2}
        - b^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        + g^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{q}^{\text{f}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{qf}]\\
    & \label{eq:ACOPF:ohm_pt}
        g^{tt}_{e} \mathbf{v}_{j}^{2}
        + g^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        - b^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{p}^{\text{t}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{pt}]\\
    & \label{eq:ACOPF:ohm_qt}
        -b^{tt}_{e} \mathbf{v}_{j}^{2}
        - b^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        - g^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{q}^{\text{t}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{qt}]\\
    % Thermal limits
    & \label{eq:ACOPF:sm_f}
        (\mathbf{p}^{\text{f}}_{e})^{2} + (\mathbf{q}^{\text{f}}_{e})^{2} \leq \bar{s}_{e}^{2}
        && \forall e \in \mathcal{E}
        && [\nu^{f}]\\
    & \label{eq:ACOPF:sm_t}
        (\mathbf{p}^{\text{t}}_{e})^{2} + (\mathbf{q}^{\text{t}}_{e})^{2} \leq \bar{s}_{e}^{2}
        && \forall e \in \mathcal{E}
        && [\nu^{t}]\\
    % Voltage angle deviation
    & \label{eq:ACOPF:va_diff}
        \underline{\Delta} \theta_{e} \leq \Delta \theta_{e}
        \leq \bar{\Delta} \theta_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{\Delta \theta}]\\
    % Variable bounds
    & \label{eq:ACOPF:vm_bounds}
        \underline{v}_{i} \leq \mathbf{v}_{i} \leq \bar{v}_{i}, 
        && \forall i \in \mathcal{N}
        && [\mu^{v}]\\ 
    & \label{eq:ACOPF:pg_bounds}
        \underline{p}^{g}_{i} \leq \mathbf{p}^{\text{g}}_{i} \leq \bar{p}^{g}_{i}, 
        && \forall i \in \mathcal{G}
        && [\mu^{pg}]\\
    & \label{eq:ACOPF:qg_bounds}
        \underline{q}^{g}_{i} \leq p^{q}_{i} \leq \bar{q}^{g}_{i},
        && \forall i \in \mathcal{G}
        && [\mu^{qg}]\\
    & \label{eq:ACOPF:pf_bounds}
        -\bar{s}_{e} \leq \mathbf{p}^{\text{f}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{pf}]\\
    & \label{eq:ACOPF:qf_bounds}
        -\bar{s}_{e} \leq \mathbf{q}^{\text{f}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{qf}]\\
    & \label{eq:ACOPF:pt_bounds}
        -\bar{s}_{e} \leq \mathbf{p}^{\text{t}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{pt}]\\
    & \label{eq:ACOPF:qt_bounds}
        -\bar{s}_{e} \leq \mathbf{q}^{\text{t}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{qt}]
\end{align}
```
where ``\Delta \theta_{e} = \boldsymbol{\theta}_{i} - \boldsymbol{\theta}_{j}`` for branch ``e = (i, j) \in \mathcal{E}``.


### Variables

* ``\mathbf{v} \in \mathbb{R}^{N}``: nodal voltage magnitude
* ``\boldsymbol{\theta} \in \mathbb{R}^{N}``: nodal voltage angle
* ``\mathbf{p}^{\text{g}} \in \mathbb{R}^{G}``: active power dispatch
* ``\mathbf{q}^{\text{g}} \in \mathbb{R}^{G}``: reactive power dispatch
* ``\mathbf{p}^{\text{f}} \in \mathbb{R}^{E}``: active power flow "from"
* ``\mathbf{q}^{\text{f}} \in \mathbb{R}^{E}``: reactive power flow "from"
* ``\mathbf{p}^{\text{t}} \in \mathbb{R}^{E}``: active power flow "to"
* ``\mathbf{q}^{\text{t}} \in \mathbb{R}^{E}``: reactive power flow "to"

### Objective

The objective function ``\eqref{eq:ACOPF:objective}`` minimizes the cost of active power generation.

!!! todo
    OPFGenerator currently supports only linear cost functions.
    Support for quadratic functions is planned for a later stage; please open an issue if 
    you would like to request this feature.

### Constraints

* ``\eqref{eq:ACOPF:slack}``: this constraint fixes the voltage angle of the 
    reference (slack) bus to zero.
* ``\eqref{eq:ACOPF:kcl_p}-\eqref{eq:ACOPF:kcl_q}``:
    active and reactive Kirchhoff's current law.
* ``\eqref{eq:ACOPF:ohm_pf}-\eqref{eq:ACOPF:ohm_qt}``:
    Ohm's law expressing power flows as a function of nodal voltages.
* ``\eqref{eq:ACOPF:sm_f}-\eqref{eq:ACOPF:sm_t}``: Thermal limits
* ``\eqref{eq:ACOPF:va_diff}``: Voltage angle deviation constraints.
* ``\eqref{eq:ACOPF:vm_bounds}``: Nodal voltage angle limits
* ``\eqref{eq:ACOPF:pg_bounds}-\eqref{eq:ACOPF:qg_bounds}``: Active/reactive generation limits
* ``\eqref{eq:ACOPF:pf_bounds}-\eqref{eq:ACOPF:qt_bounds}``: Power flow bounds, derived from thermal limits

!!! info
    Although power flow variable bounds``\eqref{eq:ACOPF:pf_bounds}-\eqref{eq:ACOPF:qt_bounds}``
    are redundant with thermal limits ``\eqref{eq:ACOPF:sm_f}-\eqref{eq:ACOPF:sm_t}``, 
    their inclusion improves the performance of interior-point solvers like Ipopt.


## Data format

### Primal solution

| Variable                  | Data | Size  | Description 
|:--------------------------|:-----|:------|:----------------------------------|
| ``\mathbf{v}``            | `vm` | ``N`` | Nodal voltage magnitude
| ``\boldsymbol{\theta}``   | `va` | ``N`` | Nodal voltage angle
| ``\mathbf{p}^{\text{g}}`` | `pg` | ``G`` | Active power generation
| ``\mathbf{q}^{\text{g}}`` | `pg` | ``G`` | Reactive power generation
| ``\mathbf{p}^{\text{f}}`` | `pf` | ``E`` | Active power flow (from)
| ``\mathbf{p}^{\text{t}}`` | `pt` | ``E`` | Active power flow (to)
| ``\mathbf{q}^{\text{f}}`` | `qf` | ``E`` | Reactive power flow (from)
| ``\mathbf{q}^{\text{t}}`` | `qt` | ``E`` | Reactive power flow (to)

### Dual solution

| Constraint                     | Data | Size  | Domain 
|:-------------------------------|:-----|:------|:----------------------------------|
| ``\eqref{eq:ACOPF:slack}``     | `slack_bus`  | ``N`` | Nodal voltage magnitude
| ``\eqref{eq:ACOPF:kcl_p}``     | `kcl_p`      | ``N`` | Nodal voltage angle
| ``\eqref{eq:ACOPF:kcl_q}``     | `kcl_q`      | ``N`` | Active power generation
| ``\eqref{eq:ACOPF:ohm_pf}``    | `ohm_pf`     | ``E`` | Reactive power generation
| ``\eqref{eq:ACOPF:ohm_qf}``    | `ohm_qf` | ``E`` | Active power flow (to)
| ``\eqref{eq:ACOPF:ohm_pt}``    | `ohm_pt` | ``E`` | Active power flow (from)
| ``\eqref{eq:ACOPF:ohm_qt}``    | `ohm_qt` | ``E`` | Reactive power flow (from)
| ``\eqref{eq:ACOPF:sm_f}``      | `sm_fr` | ``E`` | Reactive power flow (to)
| ``\eqref{eq:ACOPF:sm_t}``      | `sm_to` | ``E`` | Reactive power flow (to)
| ``\eqref{eq:ACOPF:va_diff}``   | `va_diff` | ``E`` | Reactive power flow (to)
| ``\eqref{eq:ACOPF:vm_bounds}`` (lower) | `vm_lb` | ``N`` | Voltage magnitude lower bounds
| ``\eqref{eq:ACOPF:vm_bounds}`` (upper) | `vm_ub` | ``N`` | Voltage magnitude upper bounds
| ``\eqref{eq:ACOPF:pg_bounds}`` (lower) | `pg_lb` | ``G`` | Active power generation lower bound
| ``\eqref{eq:ACOPF:pg_bounds}`` (upper) | `pg_ub` | ``G`` | Active power generation upper bound
| ``\eqref{eq:ACOPF:qg_bounds}`` (lower) | `qg_lb` | ``G`` | Reactive power generation lower bound
| ``\eqref{eq:ACOPF:qg_bounds}`` (upper) | `qg_ub` | ``G`` | Reactive power generation upper bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (lower) | `pf_lb` | ``E`` | Active power flow (from) lower bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (upper) | `pf_ub` | ``E`` | Active power flow (from) upper bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (lower) | `qf_lb` | ``E`` | Reactive power flow (from) lower bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (upper) | `qf_ub` | ``E`` | Reactive power flow (from) upper bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (lower) | `pt_lb` | ``E`` | Active power flow (to) lower bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (upper) | `pt_ub` | ``E`` | Active power flow (to) upper bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (lower) | `qt_lb` | ``E`` | Reactive power flow (to) lower bound
| ``\eqref{eq:ACOPF:pf_bounds}`` (upper) | `qt_ub` | ``E`` | Reactive power flow (to) upper bound