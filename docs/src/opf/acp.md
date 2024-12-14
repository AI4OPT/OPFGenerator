# AC-OPF

The ACOPF model considered in OPFGenerator is presented below.

```math
\begin{align}
    \min \quad 
    & \sum_{g \in \mathcal{G}} c_{g} \mathbf{p}^{\text{g}}_{g} + c^{0}_{g}\\
    \text{s.t.} \quad
    & \theta_{i_{s}} = 0\\
    & \sum_{g \in \mathcal{G}_{i}} \mathbf{p}^{\text{g}}_{g}
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{p}^{\text{f}}_{e}
        - \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{p}^{\text{t}}_{e}
        - g^{s}_{i} \mathbf{v}_{i}^{2}
        =
        \sum_{l \in \mathcal{L}_{i}} p^{d}_{l}
        & \forall i \in \mathcal{N}
        &&& [\lambda^{p}]\\
    & \sum_{g \in \mathcal{G}_{i}} \mathbf{q}^{\text{g}}_{g}
        - \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{q}^{\text{f}}_{e}
        - \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{q}^{\text{t}}_{e}
        + b^{s}_{i} \mathbf{v}_{i}^{2}
        =
        \sum_{l \in \mathcal{L}_{i}} q^{d}_{l}
        & \forall i \in \mathcal{N}
        &&& [\lambda^{q}]\\
    % Ohm's law
    & g^{ff}_{e} \mathbf{v}_{i}^{2}
        + g^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        + b^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{p}^{\text{f}}_{e} = 0
        & \forall e \in \mathcal{E}
        &&& [\lambda^{pf}]\\
    & -b^{ff}_{e} \mathbf{v}_{i}^{2}
        - b^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        + g^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{q}^{\text{f}}_{e} = 0
        & \forall e \in \mathcal{E}
        &&& [\lambda^{qf}]\\
    & g^{tt}_{e} \mathbf{v}_{j}^{2}
        + g^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        - b^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{p}^{\text{t}}_{e} = 0
        & \forall e \in \mathcal{E}
        &&& [\lambda^{pt}]\\
    & -b^{tt}_{e} \mathbf{v}_{j}^{2}
        - b^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
        - g^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
        - \mathbf{q}^{\text{t}}_{e} = 0
        & \forall e \in \mathcal{E}
        &&& [\lambda^{qt}]\\
    % Thermal limits
    & (\mathbf{p}^{\text{f}}_{e})^{2} + (\mathbf{q}^{\text{f}}_{e})^{2} \leq \bar{s}_{e}^{2}
        & \forall e \in \mathcal{E}
        &&& [\nu^{f}]\\
    & (\mathbf{p}^{\text{t}}_{e})^{2} + (\mathbf{q}^{\text{t}}_{e})^{2} \leq \bar{s}_{e}^{2}
        & \forall e \in \mathcal{E}
        &&& [\nu^{t}]\\
    % Voltage angle deviation
    & \underline{\Delta} \theta_{e} \leq \Delta \theta_{e}
        \leq \bar{\Delta} \theta_{e}
        & \forall e \in \mathcal{E}
        &&& [\mu^{\Delta \theta}]\\
    % Variable bounds
    & \underline{v}_{i} \leq \mathbf{v}_{i} \leq \bar{v}_{i}, 
        & \forall i \in \mathcal{N}
        &&& [\mu^{v}]\\ 
    & \underline{p}^{g}_{i} \leq \mathbf{p}^{\text{g}}_{i} \leq \bar{p}^{g}_{i}, 
        & \forall i \in \mathcal{G}
        &&& [\mu^{pg}]\\
    & \underline{q}^{g}_{i} \leq p^{q}_{i} \leq \bar{q}^{g}_{i},
        & \forall i \in \mathcal{G}
        &&& [\mu^{qg}]\\
    & -\bar{s}_{e} \leq \mathbf{p}^{\text{f}}_{e}, \mathbf{q}^{\text{f}}_{e}, \mathbf{p}^{\text{t}}_{e}, \mathbf{q}^{\text{t}}_{e} \leq \bar{s}_{e}
        & \forall e \in \mathcal{E}
        &&& [\mu^{pf,qf,pt,qt}]
\end{align}
```
where ``\Delta \theta_{e} = \theta_{i} - \theta_{j}`` for branch ``e = (i, j) \in \mathcal{E}``.


## Variables

* ``v \in \mathbb{R}^{N}``: nodal voltage magnitude
* ``\theta \in \mathbb{R}^{N}``: nodal voltage angle
* ``\mathbf{p}^{\text{g}} \in \mathbb{R}^{G}``: active power dispatch
* ``\mathbf{q}^{\text{g}} \in \mathbb{R}^{G}``: reactive power dispatch
* ``\mathbf{p}^{\text{f}} \in \mathbb{R}^{E}``: active power flow "from"
* ``\mathbf{q}^{\text{f}} \in \mathbb{R}^{E}``: reactive power flow "from"
* ``\mathbf{p}^{\text{t}} \in \mathbb{R}^{E}``: active power flow "to"
* ``\mathbf{q}^{\text{t}} \in \mathbb{R}^{E}``: reactive power flow "to"

## Objective

The objective function minimizes the cost of active power generation.
OPFGenerator currently supports only linear cost functions.

## Constraints


### Slack bus reference angle

The slack bus reference angle constraint sets the voltage angle to zero at the slack bus.
```math
\theta_{i_{s}} = 0
```
where ``i_{s}`` is the index of the slack bus.

### Kirchhoff current law

Kirchhoff's current law (also known as "nodal power balance") for bus ``i \in \mathcal{N}``
reads, for active power
```math
\sum_{g \in \mathcal{G}_{i}} \mathbf{p}^{\text{g}}_{g}
- \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{p}^{\text{f}}_{e}
- \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{p}^{\text{t}}_{e}
- g^{s}_{i} \mathbf{v}_{i}^{2}
=
\sum_{l \in \mathcal{L}_{i}} p^{d}_{l}
```
and for reactive power
```math
\sum_{g \in \mathcal{G}_{i}} \mathbf{q}^{\text{g}}_{g}
- \sum_{e \in \mathcal{E}^{+}_{i}} \mathbf{q}^{\text{f}}_{e}
- \sum_{e \in \mathcal{E}^{-}_{i}} \mathbf{q}^{\text{t}}_{e}
+ b^{s}_{i} \mathbf{v}_{i}^{2}
=
\sum_{l \in \mathcal{L}_{i}} q^{d}_{l}
```

### Ohm's law

Active and reactive power flows on branch ``e \in \mathcal{E}``
    are obtained via Ohm's law as follows
```math
\begin{align*}
g^{ff}_{e} \mathbf{v}_{i}^{2}
    + g^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
    + b^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
    - \mathbf{p}^{\text{f}}_{e} &= 0\\
-b^{ff}_{e} \mathbf{v}_{i}^{2}
    - b^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
    + g^{ft}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
    - \mathbf{q}^{\text{f}}_{e} &= 0\\
g^{tt}_{e} \mathbf{v}_{j}^{2}
    + g^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
    - b^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
    - \mathbf{p}^{\text{t}}_{e} &= 0\\
-b^{tt}_{e} \mathbf{v}_{j}^{2}
    - b^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \cos (\Delta \theta_{e})
    - g^{tf}_{e} \mathbf{v}_{i} \mathbf{v}_{j} \sin (\Delta \theta_{e})
    - \mathbf{q}^{\text{t}}_{e} &= 0
\end{align*}
```
where ``\Delta \theta_{e} = \theta_{i} - \theta_{j}`` and branch ``e`` is from bus ``i`` to bus ``j``.

### Thermal limits

### Voltage angle difference limits

### Variable bounds

**Nodal voltage magnitude**
```math
\underline{v}_{i} \leq \mathbf{v}_{i} \leq \bar{v}_{i}, \quad \forall i \in \mathcal{N}
```

* Active/reactive power dispatch
```math
\underline{p}^{g}_{i} \leq \mathbf{p}^{\text{g}}_{i} \leq \bar{p}^{g}_{i}, \quad \forall i \in \mathcal{G}\\
\underline{q}^{g}_{i} \leq p^{q}_{i} \leq \bar{q    }^{g}_{i}, \quad \forall i \in \mathcal{G}
```
  Minimum and/or maximum limits for active/reactive power dispatch may take negative values.

* Bounds on power flows are derived from thermal limits
```math
-S_{e} \leq \mathbf{p}^{\text{f}}_{e}, \mathbf{q}^{\text{f}}_{e}, \mathbf{p}^{\text{t}}_{e}, \mathbf{q}^{\text{t}}_{e} \leq S_{e}, \quad \forall e \in \mathcal{E}
```

!!! info
    Although these lower/upper bounds on active/reactive power flows
     are redundant with thermal constraints, their presence improves the performance
    of interior-point solvers like Ipopt.

## Data format



