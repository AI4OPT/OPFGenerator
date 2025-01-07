# SOC-OPF

The SOC-OPF model considered in OPFGenerator is presented below.

### Definitions

The _second order cone_ and _rotated second-order cone_ or order ``n`` are defined as
```math
\begin{align*}
\mathcal{Q}^{n} &= \left\{
    x \in \mathbb{R}^{n}   
\ \middle| \
    x_{1} \geq \sqrt{x_{2}^{2} + ... + x_{n}^{2}}
\right\}\\
\mathcal{Q}_{r}^{n} &= \left\{
    x \in \mathbb{R}^{n}   
\ \middle| \
    2 x_{1} x_{2} \geq x_{3}^{2} + ... + x_{n}^{2},
    x_{1}, x_{2} \geq 0
\right\}.
\end{align*}
```

## Mathematical formulation

The SOC-OPF formulation in OPFGenerator is the Jabr relaxation of AC-OPF.
The formulation is obtained through the change of variable
```math
\begin{align*}
\mathbf{w}_{i} &= \mathbf{v}_{i}^{2} 
    && \forall i \in \mathcal{N}\\
\mathbf{w}^{\text{r}}_{e} &= \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e}) 
    && \forall e = (i, j) \in \mathcal{E}\\
\mathbf{w}^{\text{i}}_{e} &= \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e}) 
    && \forall e = (i, j) \in \mathcal{E}
\end{align*}
```
where ``\Delta \theta_{e} = \theta_{i} - \theta_{j}``.
Note that ``\mathbf{w}^{\text{r}}_{e}`` and ``\mathbf{w}^{\text{i}}_{e}`` correspond to 
the real and imaginary parts of the complex voltage product ``V_{i}V_{j}^{*}``, respectively.
This transformation implies the non-convex equality
```math
(\mathbf{w}^{\text{r}})^{2} + (\mathbf{w}^{\text{i}})^{2} = \mathbf{v}_{i}^{2} \times \mathbf{v}_{j}^{2} = \mathbf{w}_{i} \times \mathbf{w}_{j}
```
The Jabr relaxation is obtained by relaxing this into the convex constraint
```math
(\mathbf{w}^{\text{r}})^{2} + (\mathbf{w}^{\text{i}})^{2} \leq \mathbf{w}_{i} \times \mathbf{w}_{j}
```

The resulting SOC-OPF problem has the form
```math
\begin{align}
    \min \quad 
    & \label{eq:SOCOPF:objective}
        \sum_{g \in \mathcal{G}} c_{g} \mathbf{p}^{\text{g}}_{g} + c^{0}_{g}\\
    \text{s.t.} \quad
    & \label{eq:SOCOPF:kcl_p}
        \sum_{g \in \mathcal{G}_{i}} \mathbf{p}^{\text{g}}_{g}
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{p}^{\text{f}}_{e}
        - \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{p}^{\text{t}}_{e}
        - g^{s}_{i} \mathbf{w}_{i}
        =
        \sum_{l \in \mathcal{L}_{i}} p^{d}_{l}
        && \forall i \in \mathcal{N}
        && [\lambda^{p}]\\
    & \label{eq:SOCOPF:kcl_q}
        \sum_{g \in \mathcal{G}_{i}} \mathbf{q}^{\text{g}}_{g}
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{q}^{\text{f}}_{e}
        - \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{q}^{\text{t}}_{e}
        + b^{s}_{i} \mathbf{w}_{i}
        =
        \sum_{l \in \mathcal{L}_{i}} q^{d}_{l}
        && \forall i \in \mathcal{N}
        && [\lambda^{q}]\\
    % Ohm's law
    & \label{eq:SOCOPF:ohm_pf}
        g^{ff}_{e} \mathbf{w}_{i}
        + g^{ft}_{e} \mathbf{w}^{\text{r}}_{e}
        + b^{ft}_{e} \mathbf{w}^{\text{i}}_{e}
        - \mathbf{p}^{\text{f}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{pf}]\\
    & \label{eq:SOCOPF:ohm_qf}
        -b^{ff}_{e} \mathbf{w}_{i}
        - b^{ft}_{e} \mathbf{w}^{\text{r}}_{e}
        + g^{ft}_{e} \mathbf{w}^{\text{i}}_{e}
        - \mathbf{q}^{\text{f}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{qf}]\\
    & \label{eq:SOCOPF:ohm_pt}
        g^{tt}_{e} \mathbf{w}_{j}
        + g^{tf}_{e} \mathbf{w}^{\text{r}}_{e}
        - b^{tf}_{e} \mathbf{w}^{\text{i}}_{e}
        - \mathbf{p}^{\text{t}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{pt}]\\
    & \label{eq:SOCOPF:ohm_qt}
        -b^{tt}_{e} \mathbf{w}_{j}
        - b^{tf}_{e} \mathbf{w}^{\text{r}}_{e}
        - g^{tf}_{e} \mathbf{w}^{\text{i}}_{e}
        - \mathbf{q}^{\text{t}}_{e} = 0
        && \forall e \in \mathcal{E}
        && [\lambda^{qt}]\\
    % Jabr constraints
    & \label{eq:SOCOPF:jabr}
        \left(
            \frac{\mathbf{w}_{i}}{\sqrt{2}},
            \frac{\mathbf{w}_{j}}{\sqrt{2}},
            \mathbf{w}^{\text{r}}_{e},
            \mathbf{w}^{\text{i}}_{e}
        \right)
        \in \mathcal{Q}_{r}^{4}
        && \forall e \in \mathcal{E}
        && [\omega]\\
    % Thermal limits
    & \label{eq:SOCOPF:sm_f}
        (\bar{s}_{e}, \mathbf{p}^{\text{f}}_{e}, \mathbf{q}^{\text{f}}_{e})
        \in \mathcal{Q}^{3}
        && \forall e \in \mathcal{E}
        && [\nu^{f}]\\
    & \label{eq:SOCOPF:sm_t}
        (\bar{s}_{e}, \mathbf{p}^{\text{t}}_{e}, \mathbf{q}^{\text{t}}_{e})
        \in \mathcal{Q}^{3}
        && \forall e \in \mathcal{E}
        && [\nu^{t}]\\
    % Voltage angle deviation
    & \label{eq:SOCOPF:va_diff}
        \underline{\Delta} \theta_{e} \leq \Delta \theta_{e}
        \leq \bar{\Delta} \theta_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{\Delta \theta}]\\
    % Variable bounds
    & \label{eq:SOCOPF:wm_bounds}
        \underline{v}_{i}^{2} \leq \mathbf{w}_{i} \leq \bar{v}_{i}^{2}, 
        && \forall i \in \mathcal{N}
        && [\mu^{w}]\\ 
    & \label{eq:SOCOPF:pg_bounds}
        \underline{p}^{g}_{i} \leq \mathbf{p}^{\text{g}}_{i} \leq \bar{p}^{g}_{i}, 
        && \forall i \in \mathcal{G}
        && [\mu^{pg}]\\
    & \label{eq:SOCOPF:qg_bounds}
        \underline{q}^{g}_{i} \leq \mathbf{q}^{\text{g}}_{i} \leq \bar{q}^{g}_{i},
        && \forall i \in \mathcal{G}
        && [\mu^{qg}]\\
    & \label{eq:SOCOPF:wr_bounds}
        \underline{w}^{\text{r}}_{e} \leq \mathbf{w}^{\text{r}}_{e} \leq \bar{w}^{\text{r}}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{wr}]\\
    & \label{eq:SOCOPF:wi_bounds}
        \underline{w}^{\text{i}}_{e} \leq \mathbf{w}^{\text{i}}_{e} \leq \bar{w}^{\text{i}}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{wi}]\\
    & \label{eq:SOCOPF:pf_bounds}
        -\bar{s}_{e} \leq \mathbf{p}^{\text{f}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{pf}]\\
    & \label{eq:SOCOPF:qf_bounds}
        -\bar{s}_{e} \leq \mathbf{q}^{\text{f}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{qf}]\\
    & \label{eq:SOCOPF:pt_bounds}
        -\bar{s}_{e} \leq \mathbf{p}^{\text{t}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{pt}]\\
    & \label{eq:SOCOPF:qt_bounds}
        -\bar{s}_{e} \leq \mathbf{q}^{\text{t}}_{e} \leq \bar{s}_{e}
        && \forall e \in \mathcal{E}
        && [\mu^{qt}]
\end{align}
```
where ``\Delta \theta_{e} = \theta_{i} - \theta_{j}`` for branch ``e = (i, j) \in \mathcal{E}``.
Similarly, in constraints ``\eqref{eq:SOCOPF:ohm_pf}-\eqref{eq:SOCOPF:jabr}``,
indices ``i`` and ``j`` denote the source and destination of branch ``e \in \mathcal{E}``, respectively.

### Variables

* ``\mathbf{w} \in \mathbb{R}^{N}``: squared nodal voltage magnitude
* ``\mathbf{p}^{\text{g}} \in \mathbb{R}^{G}``: active power dispatch
* ``\mathbf{q}^{\text{g}} \in \mathbb{R}^{G}``: reactive power dispatch
* ``\mathbf{w}^{\text{r}} \in \mathbb{R}^{E}``: real part of voltage product
* ``\mathbf{w}^{\text{i}} \in \mathbb{R}^{E}``: imaginary part of voltage product
* ``\mathbf{p}^{\text{f}} \in \mathbb{R}^{E}``: active power flow "from"
* ``\mathbf{q}^{\text{f}} \in \mathbb{R}^{E}``: reactive power flow "from"
* ``\mathbf{p}^{\text{t}} \in \mathbb{R}^{E}``: active power flow "to"
* ``\mathbf{q}^{\text{t}} \in \mathbb{R}^{E}``: reactive power flow "to"

### Objective

The objective function minimizes the cost of active power generation.

### Constraints

* ``\eqref{eq:SOCOPF:kcl_p}-\eqref{eq:SOCOPF:kcl_q}``:
    Kirchhoff's current law for active and reactive power
* ``\eqref{eq:SOCOPF:ohm_pf}-\eqref{eq:SOCOPF:ohm_qt}``:
    Ohm's law for active/reactive power flows in _from_ and _to_ directions.
* ``\eqref{eq:SOCOPF:jabr}``: Jabr constraint
* ``\eqref{eq:SOCOPF:sm_f}-\eqref{eq:SOCOPF:sm_t}``: Thermal limits
* ``\eqref{eq:SOCOPF:va_diff}``: voltage angle deviation constraints.
* ``\eqref{eq:SOCOPF:wm_bounds}``: nodal voltage angle limits
* ``\eqref{eq:SOCOPF:pg_bounds}-\eqref{eq:SOCOPF:qg_bounds}``: active/reactive generation limits
* ``\eqref{eq:SOCOPF:wr_bounds}-\eqref{eq:SOCOPF:wi_bounds}``: bounds on voltage product variables
* ``\eqref{eq:SOCOPF:pf_bounds}-\eqref{eq:SOCOPF:qt_bounds}``: power flow bounds, derived from thermal limits

!!! info
    Although power flow variable bounds``\eqref{eq:SOCOPF:pf_bounds}-\eqref{eq:SOCOPF:qt_bounds}``
    are redundant with thermal limits ``\eqref{eq:SOCOPF:sm_f}-\eqref{eq:SOCOPF:sm_t}``, 
    their inclusion improves the performance of interior-point solvers like Ipopt.

!!! tip
    When using the `SOCOPF` formulation, all convex quadratic are passed to the solver in conic form.
    To use quadratic form, e.g., when using a solver that doesn't support conic constraints,
    use `SOCOPFQuad`.

## Data format

### Primal solution

| Variable                  | Data | Size  | Description 
|:--------------------------|:-----|:------|:----------------------------------|
| ``\mathbf{w}``            | `w`  | ``N`` | Squared nodal voltage magnitude
| ``\mathbf{p}^{\text{g}}`` | `pg` | ``G`` | Active power generation
| ``\mathbf{q}^{\text{g}}`` | `pg` | ``G`` | Reactive power generation
| ``\mathbf{w}^{\text{r}}`` | `wr` | ``E`` | Voltage product variable (real part)
| ``\mathbf{w}^{\text{i}}`` | `wi` | ``E`` | Voltage product variable (imaginary part)
| ``\mathbf{p}^{\text{f}}`` | `pf` | ``E`` | Active power flow (from)
| ``\mathbf{p}^{\text{t}}`` | `pt` | ``E`` | Active power flow (to)
| ``\mathbf{q}^{\text{f}}`` | `qf` | ``E`` | Reactive power flow (from)
| ``\mathbf{q}^{\text{t}}`` | `qt` | ``E`` | Reactive power flow (to)

### Dual solution


| Constraint                      | Data         | Size  | Domain 
|:--------------------------------|:-------------|:------|:----------------------------------|
| ``\eqref{eq:SOCOPF:kcl_p}``     | `kcl_p`      | ``N`` | Kirchhoff's current law (active power)
| ``\eqref{eq:SOCOPF:kcl_q}``     | `kcl_q`      | ``N`` | Kirchhoff's current law (reactive power)
| ``\eqref{eq:SOCOPF:ohm_pf}``    | `ohm_pf`     | ``E`` | Ohm's law, active power flow (from)
| ``\eqref{eq:SOCOPF:ohm_qf}``    | `ohm_qf`     | ``E`` | Ohm's law, reactive power flow (from)
| ``\eqref{eq:SOCOPF:ohm_pt}``    | `ohm_pt`     | ``E`` | Ohm's law, active power flow (to)
| ``\eqref{eq:SOCOPF:ohm_qt}``    | `ohm_qt`     | ``E`` | Ohm's law, reactive power flow (to)
| ``\eqref{eq:SOCOPF:jabr}``      | `jabr `      | ``E \times 4`` | Jabr constraint
| ``\eqref{eq:SOCOPF:sm_f}``      | `sm_fr`      | ``E \times 3`` | Thermal limit (from)
| ``\eqref{eq:SOCOPF:sm_t}``      | `sm_to`      | ``E \times 3`` | Thermal limit (to)
| ``\eqref{eq:SOCOPF:va_diff}``   | `va_diff`    | ``E`` | Voltage angle deviation
| ``\eqref{eq:SOCOPF:wm_bounds}`` (lower) | `w_lb` | ``N`` | Squared voltage magnitude lower bounds
| ``\eqref{eq:SOCOPF:wm_bounds}`` (upper) | `w_ub` | ``N`` | Squared voltage magnitude upper bounds
| ``\eqref{eq:SOCOPF:pg_bounds}`` (lower) | `pg_lb` | ``G`` | Active power generation lower bound
| ``\eqref{eq:SOCOPF:pg_bounds}`` (upper) | `pg_ub` | ``G`` | Active power generation upper bound
| ``\eqref{eq:SOCOPF:qg_bounds}`` (lower) | `qg_lb` | ``G`` | Reactive power generation lower bound
| ``\eqref{eq:SOCOPF:qg_bounds}`` (upper) | `qg_ub` | ``G`` | Reactive power generation upper bound
| ``\eqref{eq:SOCOPF:wr_bounds}`` (lower) | `wr_lb` | ``E`` | Voltage product (real part) lower bound
| ``\eqref{eq:SOCOPF:wr_bounds}`` (upper) | `wr_ub` | ``E`` | Voltage product (real part) upper bound
| ``\eqref{eq:SOCOPF:wi_bounds}`` (lower) | `wi_lb` | ``E`` | Voltage product (imaginary part) lower bound
| ``\eqref{eq:SOCOPF:wi_bounds}`` (upper) | `wi_ub` | ``E`` | Voltage product (imaginary part) upper bound
| ``\eqref{eq:SOCOPF:pf_bounds}`` (lower) | `pf_lb` | ``E`` | Active power flow (from) lower bound
| ``\eqref{eq:SOCOPF:pf_bounds}`` (upper) | `pf_ub` | ``E`` | Active power flow (from) upper bound
| ``\eqref{eq:SOCOPF:qf_bounds}`` (lower) | `qf_lb` | ``E`` | Reactive power flow (from) lower bound
| ``\eqref{eq:SOCOPF:qf_bounds}`` (upper) | `qf_ub` | ``E`` | Reactive power flow (from) upper bound
| ``\eqref{eq:SOCOPF:pt_bounds}`` (lower) | `pt_lb` | ``E`` | Active power flow (to) lower bound
| ``\eqref{eq:SOCOPF:pt_bounds}`` (upper) | `pt_ub` | ``E`` | Active power flow (to) upper bound
| ``\eqref{eq:SOCOPF:qt_bounds}`` (lower) | `qt_lb` | ``E`` | Reactive power flow (to) lower bound
| ``\eqref{eq:SOCOPF:qt_bounds}`` (upper) | `qt_ub` | ``E`` | Reactive power flow (to) upper bound