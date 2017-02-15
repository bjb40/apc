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

margins, at(xcor2=(0(10)400))
marginsplot

reg y2 c.p##c.p c.xcor##c.xcor

margins, dydx(xcor) at(p=(1(1)20))
marginsplot

margins, at(xcor=(-20(1)20)) 
marginsplot

margins, at(xcor=(-20(1)20) p=(0,1,10.5,20))
marginsplot

margins, at(p=(0(1)20))
marginsplot

reg y2 c.a##c.a c.p##c.p c.xcor##c.xcor 

margins, at(xcor=(-20(1)20) p=(10.5) a=(10.5))
marginsplot

*@@@@@@@@@@@@@
*Checking M1
*@@@@@@@@@@@@@@

import delimited using "H:/projects/apc/output/sim_data/tdat.csv", clear

gen c = p-a

gen newc = c+20

egen cfac = cut(newc), at(-20(2)20)

reg y1 i.a i.p i.cfac

margins i.a /*why does this look wrong!?*/
marginsplot

/*this works from the yellow paper calculations*/
reg y1 c.a##c.a c.p##c.p c.a#c.p
margins, dydx(a) at(a=(1(1)20)) atmeans
marginsplot

/*specialized measurement model*/
mixed y1 a p c || c:

reg y4 c.a##c.a c.p##c.p
margins, dydx(a) at(a=(1(1)20)) atmeans
marginsplot
