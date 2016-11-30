*dev Stata 14
*APC analyees -- data simulated in R
*Bryce Bartlett


import delimited using "H:/projects/APC/output/simdat.csv", clear

regress y age age2 if age>18
