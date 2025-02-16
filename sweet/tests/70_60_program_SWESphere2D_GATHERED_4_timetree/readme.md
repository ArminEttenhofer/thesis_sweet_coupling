#### Directory for ETDSDC benchmarking

Scripts and sample plots for benchmarking SDC (ETD/IMEX) and other exponential methods. Part of master's thesis, see [here](https://www.martin-schreiber.info/data/student_projects/2024_MA_elizaveta_boriskova.pdf).

Included tests:
- Galewsky test case
- Williamson's test case
- Strong scaling for LRZ linux cluster of SCCS ("tiny" partition)
- order verification (on Galewsky, for any SDC `order=min(sweeps+1, substeps)`)

`etdsdc_wil_15_M128` folser also contains copies of total geopotential plotting scripts (by @x-JCalda). For linux cluster access please ask @x-KGadda.
