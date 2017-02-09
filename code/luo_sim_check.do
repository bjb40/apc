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



import delimited using "H:/projects/apc/output/sim_data/cordat.csv", clear

gen xfactor = round(xcor+20)

reg y2 i.p i.xfactor

margins i.xfactor
marginsplot

gen p2 = p*p
gen xcor2 = xcor*xcor

reg y2 p xcor p2 xcor2

margins, at(xcor=(-20(1)20))
marginsplot
