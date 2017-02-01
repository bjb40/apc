*dev Stata 14
*checks analysis
*Bryce Bartlett

import delimited using "H:/projects/apc/output/sim_data/tdat.csv", clear

gen cfactor = c+20

reg y2 i.p i.cfactor

*margins i.p
*marginsplot

margins i.cfactor
marginsplot

margins i.p
marginsplot


import delimited using "H:/projects/apc/output/sim_data/cordat.csv", clear

gen xfactor = round(xcor+20)

reg y2 i.p i.xfactor

margins i.xfactor
marginsplot
