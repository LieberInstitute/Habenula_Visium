From `code/10_HD_bin_level/cell_environment/22_nearby_transcription.R` at k = 10
we get the suspicious result that cluster 0 strongly colocalizes with habenula
but is fairly negatively correlated transcriptionally with all habenula cell
types. There are several other similar cases of transcriptional disagreement that
seem unlikely. This directory attempts to approach the data in a different way
to either verify or find an issue with the current results.

This is the idea:

To test whether that weird result with cluster 0 at extracellular FICTURE k = 10 is correct (specifically that transcription is anti-correlated with habenula), we can register the k = 10 results against a "cellular FICTURE" version of the data. That is:

- For each cell, label that cell with the FICTURE cluster most dominant around that cell. Then take the cellular (not extracellular) transcriptional data and calculate markers for each FICTURE cluster, followed by registration

Given that k = 10 cluster 0 seems to reliably locate in the habenula, if the weird result is correct, we should see anticorrelation between extracellular and cellular cluster 0.
