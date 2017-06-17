###equivalents
#written on paper reducing secenario 1 to a and p only

age.x=1:20
p=10.5 #(held at mean for interaction


plot(age.x,-0.05*age.x-0.0115*age.x^2-0.003*age.x*p,
     type='l',lty=2,ylim=c(-7,4))
lines(age.x,.3*age.x-0.01*age.x^2)
