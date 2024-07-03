# Notations

## Input parameters

Sets:
* ``\mathcal{N} = \{1, ..., N\}``: Set of buses
* ``\mathcal{E} = \{1, ..., E\}``: Set of branches
* ``\mathcal{E}^{+}_{i}``: Set of branches _leaving_ bus ``i \in \mathcal{N}``
* ``\mathcal{E}^{-}_{i}``: Set of branches _entering_ bus ``i \in \mathcal{N}``
* ``\mathcal{G} = \{1, ..., G\}``: Set of generators
* ``\mathcal{G}_{i} = \{1, ..., G_{i}\}``: Set of generators at bus ``i \in \mathcal{N}``
* ``\mathcal{L} = \{1, ..., L\}``: Set of loads at bus ``i \in \mathcal{N}``
* ``\mathcal{S} = \{1, ..., S\}``: Set of shunts
* ``\mathcal{S} = \{1, ..., S_{i}\}``: Set of shunts at bus ``i \in \mathcal{N}``.



Network data:
* ``g^{s}_{s} + \mathbf{j} \, b^{s}_{s}``: complex admittance of shunt ``s \in \mathcal{S}``
* ``g_{e} + \mathbf{j} \, b_{e}``: complex admittance of branch ``e \in \mathcal{E}``
