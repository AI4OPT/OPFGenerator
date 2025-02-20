# SDP-OPF

The SDP-OPF model considered in OPFGenerator is presented below.

### Definitions
Let ``\im`` be the imaginary unit. Let ``\cdot^{*}`` be the complex conjugate operator.

The _second order cone_ of order ``n`` is defined as
```math
\mathcal{Q}^{n} = \left\{
    x \in \mathbb{R}^{n}   
\ \middle| \
    x_{1} \geq \sqrt{x_{2}^{2} + ... + x_{n}^{2}}
\right\}.
```


The set of ``n \times n`` symmetric matrices is
```math
\symmat^{n} = \left\{
    X \in \mathbb{R}^{n \times n}
\ \middle| \
    X = X^{\top}
\right\}.
```
The set of ``n \times n`` skew-symmetric matrices is
```math
\skewmat^{n} = \left\{
    X \in \mathbb{R}^{n \times n}
\ \middle| \
    X = -X^{\top}
\right\}.
```
The real _positive semidefinite cone_ of order ``n`` is defined as
```math
\psdsymmat^{n} = \left\{
    X \in \symmat^{n}
\ \middle| \
    X \succeq 0
\right\}.
```
which is the set of positive semidefinite ``n \times n`` symmetric matrices in the real space.

A complex square matrix is Hermitian if it is equal to its conjugate transpose.
The set of ``n \times n`` Hermitian matrices is
```math
\hermat^{n} = \left\{
    X \in \mathbb{C}^{n \times n}
\ \middle| \
    X = X^{H}
\right\}.
```
Note that the diagonal entries of a Hermitian matrix are real.
The complex _positive semidefinite cone_ of order ``n`` is defined as
```math
\psdhermat^{n} = \left\{
    X \in \hermat^{n}
\ \middle| \
    X \succeq 0
\right\}.
```
which is the set of positive semidefinite ``n \times n`` Hermitian matrices in the complex space.

## Mathematical formulation

The SDP-OPF formulation in OPFGenerator is obtained through the change of variable ``\Wsdp = \V \V^{*}``; i.e., ``\Wsdp_{ij} = \V_{i} \V_{j}^{*}``.
This change allows formulating AC-OPF so that constraints involving ``\Wsdp`` are linear in ``\Wsdp``.
It is known that ``\Wsdp = \V \V^{*}`` for some ``\V`` if and only if ``\Wsdp \in \psdhermat^{N}`` and ``rank(\Wsdp) = 1``. By relaxing the rank constraint, a complex semidefinite program (SDP) relaxation of AC-OPF is obtained, where a positive semidefinite constraint requires ``\Wsdp \in \psdhermat^{N}``.

The complex SDP can be converted into a real formulation.
Define ``\Wrsdp \in \symmat^{N}`` and ``\Wisdp \in \skewmat^{N}``.
Then, it is known that
```math
\begin{pmatrix}
    \Wrsdp & \Wisdp \\
    -\Wisdp & \Wrsdp
\end{pmatrix} = \begin{pmatrix}
    \Wrsdp & \Wisdp \\
    {\Wisdp}^{\top} & \Wrsdp
\end{pmatrix} \in \symmat^{2N}
```
is positive semidefinite if and only if ``\Wsdp = \Wrsdp + \im \Wisdp \in \hermat^{N}`` also is positive semidefinite in the complex space.
Note that, by definition, ``diag(\Wsdp) = diag(\Wrsdp)``.

The resulting SDP-OPF model in the real space is presented below.
```math
\begin{align}
    \min \quad
    & \label{eq:SDPOPF:objective}
        \sum_{i \in \NODES} \sum_{j \in \GENERATORS_{i}} c_j \PG_j + c_0
        \\
    \textrm{s.t.} \quad
    & \label{eq:SDPOPF:kcl_p}
        \sum_{j \in \GENERATORS_{i}} \PG_{j}
        - \sum_{e \in \EDGES^{+}_{i}} \PF_{e}
        - \sum_{e \in \EDGES^{-}_{i}} \PT_{e}
        - \GS_i \Wrsdp_{ii}
        = \sum_{j\in\LOADS_i} \PD_j
        & \forall i & \in \NODES
        % && [\lambda^{p}]
        \\
    & \label{eq:SDPOPF:kcl_q}
        \sum_{g \in \mathcal{G}_{i}} \QG_{g}
        - \sum_{e \in \EDGES^{+}_{i}} \QF_{e}
        - \sum_{e \in \EDGES^{-}_{i}} \QT_{e}
        + \BS_{i} \Wrsdp_{ii}
        = \sum_{j\in\LOADS_i} \QD_j
        & \forall i & \in \NODES
        % && [\lambda^{q}]
        \\
    % Ohm's law
    & \label{eq:SDPOPF:ohm_pf}
        \gff_{e} \Wrsdp_{ii}
        + \gft_{e} \Wrsdp_{e}
        + \bft_{e} \Wisdp_{e}
        - \PF_{e} = 0
        & \forall e = (i,j) & \in \EDGES
        % && [\lambdaPf]
        \\
    & \label{eq:SDPOPF:ohm_qf}
        -\bff_{e} \Wrsdp_{ii}
        - \bft_{e} \Wrsdp_{e}
        + \gft_{e} \Wisdp_{e}
        - \QF_{e}
        =0
        & \forall e = (i, j) & \in \EDGES
        % && [\lambdaQf]
        \\
    & \label{eq:SDPOPF:ohm_pt}
        \gtt_{e} \Wrsdp_{jj}
        + \gtf_{e} \Wrsdp_{e}
        - \btf_{e} \Wisdp_{e}
        - \PT_{e}
        = 0
        & \forall e = (i, j) & \in \EDGES
        % && [\lambdaPt]
        \\
    & \label{eq:SDPOPF:ohm_qt}
        -\btt_{e} \Wrsdp_{jj}
        - \btf_{e} \Wrsdp_{e}
        - \gtf_{e} \Wisdp_{e}
        - \QT_{e}
        =0 
        & \forall e = (i, j) & \in \EDGES
        % && [\lambdaQt]
        \\  
    % Thermal limits
    & \label{eq:SDPOPF:sm_f}
        (\bar{s}_{e}, \PF_{e}, \QF_{e}) \in \mathcal{Q}^{3}
        & \forall e & \in \EDGES 
        % && [\nuThermalfr]
        \\
    & \label{eq:SDPOPF:sm_t}
        (\bar{s}_{e}, \PT_{e}, \QT_{e}) \in \mathcal{Q}^{3}
        & \forall e & \in \EDGES 
        % && [\nuThermalto]
        \\
    % Voltage angle deviation
    & \label{eq:SDPOPF:va_diff}
        \tan(\dvamin_{e}) \Wrsdp_{e} \leq \Wisdp_{e} \leq \tan(\dvamax_{e}) \Wrsdp_{e} 
        & \forall e & \in \EDGES
        % && [\muAngleDiff]
        \\
    % Variable bounds
    &  \label{eq:SDPOPF:pg_bounds}
        \pgmin_{i} \leq \PG_{i} \leq \pgmax_{i}
        & \forall i & \in \GENERATORS
        % && [\muPg]
        \\
    & \label{eq:SDPOPF:qg_bounds}
        \qgmin_{i} \leq \QG_{i} \leq \qgmax_{i} 
        & \forall i & \in \GENERATORS
        % && [\muQg]
        \\
    & \label{eq:SDPOPF:wm_bounds}
        \vmmin_{i}^{2} \leq \Wrsdp_{ii} \leq \vmmax_{i}^{2} 
        & \forall i & \in \NODES
        % && [\muWm]
        \\
    & \label{eq:SDPOPF:pf_bounds}
        {-}\overline{S}_{e} \leq \PF_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % && [\mu^{pf}]
        \\
    & \label{eq:SDPOPF:qf_bounds}
        {-}\overline{S}_{e} \leq \QF_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % && [\mu^{qf}]
        \\
    & \label{eq:SDPOPF:pt_bounds}
        {-}\overline{S}_{e} \leq \PT_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % % && [\mu^{pt}]
        \\
    & \label{eq:SDPOPF:qt_bounds}
        {-}\overline{S}_{e} \leq \QT_{e} \leq \overline{S}_{e}
        & \forall e & \in \EDGES
        % && [\mu^{qt}]
        \\
    % PSD constraint
    & \label{eq:SDPOPF:psd}
        \begin{pmatrix}
            \Wrsdp & \Wisdp \\
            -\Wisdp & \Wrsdp
        \end{pmatrix} \in \psdsymmat^{2N}
        &
        % && [\Ssdp]
\end{align}
```

!!! info
    If there are multiple branches ``e_{1}, e_{2}, ...`` from bus ``i`` to ``j``, then ``\Wrsdp_{e_{1}}, \Wrsdp_{e_{2}},...`` will all refer to the same ``\Wrsdp_{ij}`` entry in the matrix. Similarly for ``\Wisdp``.

!!! info
    If branch ``e \in \EDGES`` is connected from bus ``i`` to ``j``, then ``\Wisdp_{e}`` refers to ``\Wisdp_{ij}``.
    If instead branch ``e`` is connected from bus ``j`` to ``i``, then ``\Wisdp_{e}`` refers to ``\Wisdp_{ji}``, which by definition is equal to ``-\Wisdp_{ij}``.


### Variables

* ``\PG \in \mathbb{R}^{G}``: active power dispatch
* ``\QG \in \mathbb{R}^{G}``: reactive power dispatch
* ``\Wrsdp \in \symmat^{N}``: real part of voltage product
* ``\Wisdp \in \skewmat^{N}``: imaginary part of voltage product
* ``\PF \in \mathbb{R}^{E}``: active power flow "from"
* ``\QF \in \mathbb{R}^{E}``: reactive power flow "from"
* ``\PT \in \mathbb{R}^{E}``: active power flow "to"
* ``\QT \in \mathbb{R}^{E}``: reactive power flow "to"

### Objective

The objective function ``\eqref{eq:SDPOPF:objective}`` minimizes the cost of active power generation.

!!! todo
    OPFGenerator currently supports only linear cost functions.
    Support for quadratic functions is planned for a later stage; please open an issue if 
    you would like to request this feature.

### Constraints

* ``\eqref{eq:SDPOPF:kcl_p}-\eqref{eq:SDPOPF:kcl_q}``:
    Kirchhoff's current law for active and reactive power
* ``\eqref{eq:SDPOPF:ohm_pf}-\eqref{eq:SDPOPF:ohm_qt}``:
    Ohm's law for active/reactive power flows in _from_ and _to_ directions
* ``\eqref{eq:SDPOPF:sm_f}-\eqref{eq:SDPOPF:sm_t}``: thermal limits
* ``\eqref{eq:SDPOPF:va_diff}``: voltage angle deviation constraints
* ``\eqref{eq:SDPOPF:pg_bounds}-\eqref{eq:SDPOPF:qg_bounds}``: active/reactive generation limits
* ``\eqref{eq:SDPOPF:wm_bounds}``: bounds on squared voltage magnitude
* ``\eqref{eq:SDPOPF:pf_bounds}-\eqref{eq:SDPOPF:qt_bounds}``: power flow bounds, derived from thermal limits
* ``\eqref{eq:SDPOPF:psd}``: positive semidefinite constraint

!!! info
    Although power flow variable bounds``\eqref{eq:SDPOPF:pf_bounds}-\eqref{eq:SDPOPF:qt_bounds}``
    are redundant with thermal limits ``\eqref{eq:SDPOPF:sm_f}-\eqref{eq:SDPOPF:sm_t}``, 
    their inclusion improves the performance of interior-point solvers like Ipopt.


## Data format

### Primal solution

| Variable                  | Data | Size  | Description 
|:--------------------------|:-----|:------|:----------------------------------|
| ``diag(\Wrsdp)``          | `w`  | ``N`` | Squared nodal voltage magnitude, the diagonal entries of ``\Wrsdp``
| ``\PG`` | `pg` | ``G`` | Active power generation
| ``\QG`` | `pg` | ``G`` | Reactive power generation
| ``\Wrsdp_{\EDGES}`` | `wr` | ``E`` | Voltage product variable (real part) corresponding to branches
| ``\Wisdp_{\EDGES}`` | `wi` | ``E`` | Voltage product variable (imaginary part) corresponding to branches
| ``\PF`` | `pf` | ``E`` | Active power flow (from)
| ``\PT`` | `pt` | ``E`` | Active power flow (to)
| ``\QF`` | `qf` | ``E`` | Reactive power flow (from)
| ``\QT`` | `qt` | ``E`` | Reactive power flow (to)

Since ``\Wrsdp`` and ``\Wisdp`` can be dense, only the off-diagonal entries that correspond to branches are extracted to save space.
Other entries only appear in constraints the positive semidefinite constraint ``\eqref{eq:SDPOPF:psd}`` and not other constraints. Therefore, values of these entries can be obtained by solving a positive semidefinite matrix completion problem to recover a full feasible matrix.

!!! info
    If there are multiple parallel branches ``e_{1}, e_{2}, ...`` from bus ``i`` to ``j``, the same entries ``\Wrsdp_{ij}`` and ``\Wisdp_{ij}`` are extracted for each of the branches.


### Dual solution


| Constraint                              | Data         | Size           | 
|:----------------------------------------|:-------------|:---------------|
| ``\eqref{eq:SDPOPF:kcl_p}``             | `kcl_p`      | ``N``          |
| ``\eqref{eq:SDPOPF:kcl_q}``             | `kcl_q`      | ``N``          |
| ``\eqref{eq:SDPOPF:ohm_pf}``            | `ohm_pf`     | ``E``          |
| ``\eqref{eq:SDPOPF:ohm_qf}``            | `ohm_qf`     | ``E``          |
| ``\eqref{eq:SDPOPF:ohm_pt}``            | `ohm_pt`     | ``E``          |
| ``\eqref{eq:SDPOPF:ohm_qt}``            | `ohm_qt`     | ``E``          |
| ``\eqref{eq:SDPOPF:sm_f}``              | `sm_fr`      | ``E \times 3`` |
| ``\eqref{eq:SDPOPF:sm_t}``              | `sm_to`      | ``E \times 3`` |
| ``\eqref{eq:SDPOPF:va_diff}``           | `va_diff`    | ``E``          |
| ``\eqref{eq:SDPOPF:pg_bounds}``         | `pg`         | ``G``          |
| ``\eqref{eq:SDPOPF:qg_bounds}``         | `qg`         | ``G``          | 
| ``\eqref{eq:SDPOPF:wm_bounds}``         | `w`          | ``N``          |
| ``\eqref{eq:SDPOPF:pf_bounds}``         | `pf`         | ``E``          |
| ``\eqref{eq:SDPOPF:qf_bounds}``         | `qf`         | ``E``          |
| ``\eqref{eq:SDPOPF:pt_bounds}``         | `pt`         | ``E``          |
| ``\eqref{eq:SDPOPF:qt_bounds}``         | `qt`         | ``E``          |
| ``\eqref{eq:SDPOPF:psd}``               | `s`, `sr`, `si`| ``N``, ``E``, ``E`` |

The dual variable corresponding to the positive semidefinite constraint ``\eqref{eq:SDPOPF:psd}`` of the real SDP-OPF is ``\Ssdp`` which is constrained to be in ``\psdsymmat^{2N}``.

By deriving the dual of the complex SDP formulation, it can be shown that any feasible ``\Ssdp`` in the real formulation has the same 4-block structure as the primal matrix, and each block has a sparsity structure corresponding to ``\EDGES``.
Without considering the complex formulation, the dual problem of the real formulation also implies a transformation to obtain such a feasible sparse 4-block matrix from any feasible ``\Ssdp`` easily.
Therefore, without loss of generality,
```math
\Ssdp =
\begin{pmatrix}
    \Srsdp & \Sisdp \\
    -\Sisdp & \Srsdp
\end{pmatrix}
```
where the off-diagonal entries of ``\Srsdp \in \symmat^{N}`` and ``\Sisdp \in \skewmat^{N}`` that don't correspond to any branch are 0.
!!! info
    Currently, there is no necessity to compute the aforementioned transformation, since interior point solvers and splitting cone solvers are observed to produce ``\Ssdp`` with the sparse 4-block structure when solving SDP-OPF, even though sparsity is not explicitly exploited.

`s`, `sr` and `si` are the diagonal entries of ``\Srsdp``, the off-diagonal entries of ``\Srsdp`` that correspond to branches, and the off-diagonal entries of ``\Sisdp`` that correspond to branches, respectively.
Unlike in the primal case, these values are sufficient to directly reconstruct the original ``\Ssdp``.