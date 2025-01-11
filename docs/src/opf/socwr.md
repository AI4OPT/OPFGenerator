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
\wm_{i} &= \VM_{i}^{2} 
    && \forall i \in \mathcal{N}\\
\wr_{e} &= \VM_{i} \VM_{j} \cos (\theta_{i} - \theta_{j}) 
    && \forall e = (i, j) \in \EDGES\\
\wi_{e} &= \VM_{i} \VM_{j} \sin (\theta_{i} - \theta_{j}) 
    && \forall e = (i, j) \in \EDGES
\end{align*}
```
Note that ``\wr_{e}`` and ``\wi_{e}`` correspond to 
the real and imaginary parts of the complex voltage product ``V_{i}V_{j}^{*}``, respectively.
This transformation implies the non-convex equality
```math
(\wr)^{2} + (\wi)^{2} = \VM_{i}^{2} \times \VM_{j}^{2} = \wm_{i} \times \wm_{j}
```
The Jabr relaxation is obtained by relaxing this into the convex constraint
```math
(\wr)^{2} + (\wi)^{2} \leq \wm_{i} \times \wm_{j}
```

The result SOC-OPF model is presented below.
```math
\begin{align}
    \min \quad 
    & \label{eq:SOCOPF:objective}
        \sum_{i \in \NODES} \sum_{j \in \GENERATORS_{i}} c_j \PG_j + c_0\\
    \text{s.t.} \quad
    & \label{eq:SOCOPF:kcl_p}
        \sum_{j\in\GENERATORS_i}\PG_j 
        - \sum_{e \in \EDGES^{+}_{i}} \PF_{e}
        - \sum_{e \in \EDGES^{-}_{i}} \PT_{e}
        - \GS_i \wm_{i}
        = \sum_{j\in\LOADS_i} \PD_j
        & \forall i & \in \NODES
        % && [\lambda^{p}]
        \\
    & \label{eq:SOCOPF:kcl_q}
        \sum_{g \in \mathcal{G}_{i}} \QG_{g}
        - \sum_{e \in \EDGES^{+}_{i}} \QF_{e}
        - \sum_{e \in \EDGES^{-}_{i}} \QT_{e}
        + \BS_{i} \wm_{i}
        = \sum_{j\in\LOADS_i} \QD_j
        & \forall i & \in \mathcal{N}
        % && [\lambda^{q}]
        \\
    % Ohm's law
    & \label{eq:SOCOPF:ohm_pf}
        \gff_{e} \wm_{i}
        + \gft_{e} \wr_{e}
        + \bft_{e} \wi_{e}
        - \PF_{e} = 0
        & \forall e & \in \EDGES
        % && [\lambda^{pf}]
        \\
    & \label{eq:SOCOPF:ohm_qf}
        -\bff_{e} \wm_{i}
        - \bft_{e} \wr_{e}
        + \gft_{e} \wi_{e}
        - \QF_{e} = 0
        & \forall e & \in \EDGES
        % && [\lambda^{qf}]
        \\
    & \label{eq:SOCOPF:ohm_pt}
        \gtt_{e} \wm_{j}
        + \gtf_{e} \wr_{e}
        - \btf_{e} \wi_{e}
        - \PT_{e} = 0
        & \forall e & \in \EDGES
        % && [\lambda^{pt}]
        \\
    & \label{eq:SOCOPF:ohm_qt}
        -\btt_{e} \wm_{j}
        - \btf_{e} \wr_{e}
        - \gtf_{e} \wi_{e}
        - \QT_{e} = 0
        & \forall e & \in \EDGES
        % && [\lambda^{qt}]
        \\
    % Jabr constraints
    & \label{eq:SOCOPF:jabr}
        \left(
            \frac{\wm_{i}}{\sqrt{2}},
            \frac{\wm_{j}}{\sqrt{2}},
            \wr_{e},
            \wi_{e}
        \right)
        \in \mathcal{Q}_{r}^{4}
        & \forall e = (i,j) & \in \EDGES
        % && [\omega]
        \\
    % Thermal limits
    & \label{eq:SOCOPF:sm_f}
        (\overline{S}_{e}, \PF_{e}, \QF_{e})
        \in \mathcal{Q}^{3}
        & \forall e & \in \EDGES
        % && [\nu^{f}]
        \\
    & \label{eq:SOCOPF:sm_t}
        (\overline{S}_{e}, \PT_{e}, \QT_{e})
        \in \mathcal{Q}^{3}
        & \forall e & \in \EDGES
        % && [\nu^{t}]
        \\
    % Voltage angle deviation
    & \label{eq:SOCOPF:va_diff}
        \tan(\dvamin_{e}) \wr_{e} \leq \wi_{e} \leq \tan(\dvamax_{e}) \wr_{e}
        & \forall e & \in \EDGES
        % && [\mu^{\Delta \theta}]
        \\
    % Variable bounds
    & \label{eq:SOCOPF:pg_bounds}
        \pgmin_{i} \leq \PG_{i} \leq \pgmax_{i}, 
        & \forall i & \in \mathcal{G}
        % && [\mu^{pg}]
        \\
    & \label{eq:SOCOPF:qg_bounds}
        \qgmin_{i} \leq \QG_{i} \leq \qgmax_{i},
        & \forall i & \in \mathcal{G}
        % && [\mu^{qg}]
        \\
    & \label{eq:SOCOPF:wm_bounds}
        \vmmin_{i}^{2} \leq \wm_{i} \leq \vmmax_{i}^{2}, 
        & \forall i & \in \mathcal{N}
        % && [\mu^{w}]
        \\ 
    & \label{eq:SOCOPF:wr_bounds}
        \wrmin_{e} \leq \wr_{e} \leq \wrmax_{e}
        & \forall e & \in \EDGES
        % && [\mu^{wr}]
        \\
    & \label{eq:SOCOPF:wi_bounds}
        \wimin_{e} \leq \wi_{e} \leq \wimax_{e}
        & \forall e & \in \EDGES
        % && [\mu^{wi}]
        \\
    & \label{eq:SOCOPF:pf_bounds}
        {-}\overline{S}_{e} \leq \PF_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % && [\mu^{pf}]
        \\
    & \label{eq:SOCOPF:qf_bounds}
        {-}\overline{S}_{e} \leq \QF_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % && [\mu^{qf}]
        \\
    & \label{eq:SOCOPF:pt_bounds}
        {-}\overline{S}_{e} \leq \PT_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % % && [\mu^{pt}]
        \\
    & \label{eq:SOCOPF:qt_bounds}
        {-}\overline{S}_{e} \leq \QT_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % && [\mu^{qt}]
\end{align}
```

### Variables

* ``\wm \in \mathbb{R}^{N}``: squared nodal voltage magnitude
* ``\PG \in \mathbb{R}^{G}``: active power dispatch
* ``\QG \in \mathbb{R}^{G}``: reactive power dispatch
* ``\wr \in \mathbb{R}^{E}``: real part of voltage product
* ``\wi \in \mathbb{R}^{E}``: imaginary part of voltage product
* ``\PF \in \mathbb{R}^{E}``: active power flow "from"
* ``\QF \in \mathbb{R}^{E}``: reactive power flow "from"
* ``\PT \in \mathbb{R}^{E}``: active power flow "to"
* ``\QT \in \mathbb{R}^{E}``: reactive power flow "to"

### Objective

The objective function ``\eqref{eq:SOCOPF:objective}`` minimizes the cost of active power generation.

!!! todo
    OPFGenerator currently supports only linear cost functions.
    Support for quadratic functions is planned for a later stage; please open an issue if 
    you would like to request this feature.

### Constraints

* ``\eqref{eq:SOCOPF:kcl_p}-\eqref{eq:SOCOPF:kcl_q}``:
    Kirchhoff's current law for active and reactive power
* ``\eqref{eq:SOCOPF:ohm_pf}-\eqref{eq:SOCOPF:ohm_qt}``:
    Ohm's law for active/reactive power flows in _from_ and _to_ directions.
* ``\eqref{eq:SOCOPF:jabr}``: Jabr constraint
* ``\eqref{eq:SOCOPF:sm_f}-\eqref{eq:SOCOPF:sm_t}``: Thermal limits
* ``\eqref{eq:SOCOPF:va_diff}``: voltage angle deviation constraints.
* ``\eqref{eq:SOCOPF:pg_bounds}-\eqref{eq:SOCOPF:qg_bounds}``: active/reactive generation limits
* ``\eqref{eq:SOCOPF:wm_bounds}``: nodal voltage angle limits
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
    !!! warning
        The `SOCOPFQuad` formulation does not support conic duality.

## Data format

### Primal solution

| Variable                  | Data | Size  | Description 
|:--------------------------|:-----|:------|:----------------------------------|
| ``\mathbf{w}``            | `w`  | ``N`` | Squared nodal voltage magnitude
| ``\PG`` | `pg` | ``G`` | Active power generation
| ``\QG`` | `pg` | ``G`` | Reactive power generation
| ``\wr`` | `wr` | ``E`` | Voltage product variable (real part)
| ``\wi`` | `wi` | ``E`` | Voltage product variable (imaginary part)
| ``\PF`` | `pf` | ``E`` | Active power flow (from)
| ``\PT`` | `pt` | ``E`` | Active power flow (to)
| ``\QF`` | `qf` | ``E`` | Reactive power flow (from)
| ``\QT`` | `qt` | ``E`` | Reactive power flow (to)

### Dual solution


| Constraint                              | Data         | Size           | 
|:----------------------------------------|:-------------|:---------------|
| ``\eqref{eq:SOCOPF:kcl_p}``             | `kcl_p`      | ``N``          |
| ``\eqref{eq:SOCOPF:kcl_q}``             | `kcl_q`      | ``N``          |
| ``\eqref{eq:SOCOPF:ohm_pf}``            | `ohm_pf`     | ``E``          |
| ``\eqref{eq:SOCOPF:ohm_qf}``            | `ohm_qf`     | ``E``          |
| ``\eqref{eq:SOCOPF:ohm_pt}``            | `ohm_pt`     | ``E``          |
| ``\eqref{eq:SOCOPF:ohm_qt}``            | `ohm_qt`     | ``E``          |
| ``\eqref{eq:SOCOPF:jabr}``              | `jabr `      | ``E \times 4`` |
| ``\eqref{eq:SOCOPF:sm_f}``              | `sm_fr`      | ``E \times 3`` |
| ``\eqref{eq:SOCOPF:sm_t}``              | `sm_to`      | ``E \times 3`` |
| ``\eqref{eq:SOCOPF:va_diff}``           | `va_diff`    | ``E``          |
| ``\eqref{eq:SOCOPF:pg_bounds}``         | `pg`         | ``G``          |
| ``\eqref{eq:SOCOPF:qg_bounds}``         | `qg`         | ``G``          | 
| ``\eqref{eq:SOCOPF:wm_bounds}``         | `w`          | ``N``          |
| ``\eqref{eq:SOCOPF:wr_bounds}``         | `wr`         | ``E``          |
| ``\eqref{eq:SOCOPF:wi_bounds}``         | `wi`         | ``E``          |
| ``\eqref{eq:SOCOPF:pf_bounds}``         | `pf`         | ``E``          |
| ``\eqref{eq:SOCOPF:qf_bounds}``         | `qf`         | ``E``          |
| ``\eqref{eq:SOCOPF:pt_bounds}``         | `pt`         | ``E``          |
| ``\eqref{eq:SOCOPF:qt_bounds}``         | `qt`         | ``E``          |
