Early on, we've been using `FICTURE` and `Banksy` for clustering. We noticed
that `FICTURE` was better able to validate known biology (like white matter
tracts), but perhaps for downstream analysis it's better to work with spatially
resolved cells, as we've done with `Banksy`. This directory contains an approach
to try to combine the benefits from both `Banksy` and `FICTURE` at the cell
level. `FICTURE` is run on bins nearby to each cell, and eventually we end up
with cells clustered by `Banksy` and a distribution across `FICTURE` clusters
describing the nearby microenvironment for each cell.
