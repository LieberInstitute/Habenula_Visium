Despite nuclear segmentations not looking terribly wrong at first, it was
discovered that H1-XQQD7C7_A1_8518 and H1-XQQD7C7_D1_8518 had many more supposed
cells than the other 3 samples, and closer inspection revealed a significant
false-positive rate in those samples. This directory contains bin2cell tests to
try to control false positives, particularly in nuclear (primary) segmentations
where inspection by eye is possible.

Later, this directory was used for separate tests, trying to improve secondary
segmentations. These segmentations had huge numbers of false positives, and in
downstream clustering, clusters of ambiguous cell type had disproportionately
many secondary cells (suggesting poor segmentation/ not getting real cells).

