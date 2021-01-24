library(VGAM)
nboot<- 1000
set.seed(12345)

#Application 1:
###############

df1<- read.table("c:/David Documents/Research/Bias correction/Topp-Leone/2021/SC16.txt", header=TRUE)
x<- df1$X
x
n<- length(x)
xbar<- mean(x)
# OBTAIN THE MLE AND MOM ESTIMATOR OF "NU"
ml<- -n/(sum(log(x))+sum(log(2-x)))
mom<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
cs<- (n-1)*ml/n
ml
cs
mom
n

# BOOTSTRAP THE STD. ERRORS FOR MLE & MOM
# ---------------------------------------
mlboot<- c()
momboot<- c()
csboot<- c()

for(ii in 1:nboot) {

xx<- sample(x,n,replace=T)
xxbar<- mean(xx)
vhatboot<- -n/(sum(log(xx))+sum(log(2-xx)))
vhatcsboot<- (n-1)*vhatboot/n
vtildboot<- uniroot(function(z) gamma(2+2*z)*(xxbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
mlboot<- c(mlboot,vhatboot)
momboot<- c(momboot,vtildboot)
csboot<- c(csboot, vhatcsboot)

}

# End of Bootstrap Loop

seml<- sd(mlboot)
semom<- sd(momboot)
secs<- sd(csboot)
c(seml,secs, semom) # Bootstrapped std. errors of MLE's and MOM

# CALCULATE THE BOOTSTRAP BIASES

mlbootbias<- mean(mlboot)-ml
mombootbias<- mean(momboot)-mom

# COMPUTE THE BIAS-ADJUSTED ESTIMATORS
mlba<- ml-mlbootbias
momba<- mom-mombootbias

c(mlba,momba)     # The bootstrap bias-adjusted MLE and MOM point estimates

# AND THEIR BOOTSTRAPPED STD. ERRORS:

# Set up the "sample" of bootstrap-bias-adjusted estimates
mlbootad<- mlboot-mlbootbias
mombootad<- momboot-mombootbias

semlboot<- sd(mlbootad)
semomboot<- sd(mombootad)
c(semlboot, semomboot)

# Plot the data and fitted models:
denml <-  function(z){ 2*ml*(1-z)*(z*(2-z))^(ml-1)
}
dencs <-  function(z){ 2*cs*(1-z)*(z*(2-z))^(cs-1)
}
denmom <-  function(z){ 2*mom*(1-z)*(z*(2-z))^(mom-1)
}
denmlba <-  function(z){ 2*mlba*(1-z)*(z*(2-z))^(mlba-1)
}
denmomba <-  function(z){ 2*momba*(1-z)*(z*(2-z))^(momba-1)
}
hist(x, prob=TRUE, main="Fig. 1: SC16 Data", xlab="Capacity Factor", breaks=12)  # Altun & Hamedani
curve(denml(x), add=TRUE, col="red",lwd=2)
curve(dencs(x), add=TRUE, col="blue",lwd=2)
curve(denmom(x), add=TRUE, col="green",lwd=2)
curve(denmlba(x), add=TRUE, col="orange",lwd=2)
curve(denmomba(x), add=TRUE, col="purple",lwd=2)

legend(0.6, 4, legend=c("Actual (Hist)","ML", "C-S", "MOM","ML-Boot", "MOM-Boot"),
       col=c("black","red", "blue", "green", "orange", "purple"), lty=c(1,1,1,1,1,1), cex=0.8)


#========================

# Goodness-of-fit tests:
# ---------------------

#Use Dn test - it has the best power
# Re-order the data
x<- sort(x)
Fml<- x^ml*(2-x)^ml
Fcs<- x^cs*(2-x)^cs
Fmom<- x^mom*(2-x)^mom
Fmlba<- x^mlba*(2-x)^mlba
Fmomba<- x^momba*(2-x)^momba
plot(x,Fml)              # Change "Fml" to "Fcs" and to "Fmom" here & below
				 # And do the same thing to get the K-S test results based on the bootstrapped estimates
ii<- rep(1:n)

lower<- ii/n-Fml
upper<- Fml-(ii-1)/n

m<- pmax(lower,upper)
#lower
#upper
#m
Dn<- max(m)
Dn
e<- ii/n
plot(x, Fml,col="blue", type="l", xlab="Capacity Factor",ylab="c.d.f.", main="Theoretical & Empirical c.d.f.'s")
lines(x,e, col="red", lty=3, lwd=2)
lines(x, Fcs, col="black", type="l", lty=2)
lines(x, Fmom, col="green", typ="l")
legend(0.05,0.97,legend=c("Theoretical (ML)", "Theoretical (C-S)", "Theoretical (MOM)","Empirical"),
       col=c("blue", "black", "green","red"), lty=c(1,2,1,3),cex=0.8,lwd=c(1,1,1,2))


####################################################################################
#Application 2:
# ============

df2<- read.table("c:/David Documents/Research/Bias correction/Topp-Leone/2021/electronics.txt", header=TRUE)
x<- df2$X/1000
x
n<- length(x)
xbar<- mean(x)
# OBTAIN THE MLE AND MOM ESTIMATOR OF "NU"
ml<- -n/(sum(log(x))+sum(log(2-x)))
mom<- uniroot(function(z) gamma(2+2*z)*(xbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
cs<- (n-1)*ml/n
ml
cs
mom
n

# BOOTSTRAP THE STD. ERRORS
# -------------------------

# BOOTSTRAP THE STD. ERRORS FOR MLE & MOM
# ---------------------------------------
mlboot<- c()
momboot<- c()
csboot<- c()
for(ii in 1:nboot) {

xx<- sample(x,n,replace=T)
xxbar<- mean(xx)
vhatboot<- -n/(sum(log(xx))+sum(log(2-xx)))
vhatcsboot<- (n-1)*vhatboot/n
vtildboot<- uniroot(function(z) gamma(2+2*z)*(xxbar-1)+4^z*(gamma(1+z))^2, lower = 0, upper = 70,
            tol = 0.0001)$root
mlboot<- c(mlboot, vhatboot)
momboot<- c(momboot,vtildboot)
csboot<- c(csboot, vhatcsboot)

}
# End of Bootstrap Loop

seml<- sd(mlboot)
semom<- sd(momboot)
secs<- sd(csboot)
c(seml,secs, semom) # Bootstrapped std. errors of MLE's and MOM

# CALCULATE THE BOOTSTRAP BIASES
mlbootbias<- mean(mlboot)-ml
mombootbias<- mean(momboot)-mom
# COMPUTE THE BIAS-ADJUSTED ESTIMATORS
mlba<- ml-mlbootbias
momba<- mom-mombootbias

c(mlba,momba)     # The bootstrap bias-adjusted MLE and MOM point estimates

# AND THEIR BOOTSTRAPPED STD. ERRORS:

# Set up the "sample" of bootstrap-bias-adjusted estimates
mlbootad<- mlboot-mlbootbias
mombootad<- momboot-mombootbias

semlboot<- sd(mlbootad)
semomboot<- sd(mombootad)
c(semlboot, semomboot)

# Plot the data and fitted models:

denml <-  function(z){ 2*ml*(1-z)*(z*(2-z))^(ml-1)
}
dencs <-  function(z){ 2*cs*(1-z)*(z*(2-z))^(cs-1)
}
denmom <-  function(z){ 2*mom*(1-z)*(z*(2-z))^(mom-1)
}
denmlba <-  function(z){ 2*mlba*(1-z)*(z*(2-z))^(mlba-1)
}
denmomba <-  function(z){ 2*momba*(1-z)*(z*(2-z))^(momba-1)
}
hist(x, prob=TRUE, main="Fig. 2: Electronic Device Data", xlab="Life-Span (/1000)", breaks=12)  # Altun & Hamedani
curve(denml(x), add=TRUE, col="red",lwd=2)
curve(dencs(x), add=TRUE, col="blue",lwd=2)
curve(denmom(x), add=TRUE, col="green",lwd=2)
curve(denmlba(x), add=TRUE, col="orange",lwd=2)
curve(denmomba(x), add=TRUE, col="purple",lwd=2)

legend(0.3, 5, legend=c("Actual (Hist)","ML", "C-S", "MOM","ML-Boot", "MOM-Boot"),
       col=c("black","red", "blue", "green", "orange", "purple"), lty=c(1,1,1,1,1,1), cex=0.8)

# Goodness-of-fit tests:
# ---------------------

#Use Dn test - it has the best power
# Re-order the data
x<- sort(x)
Fml<- x^ml*(2-x)^ml
Fcs<- x^cs*(2-x)^cs
Fmom<- x^mom*(2-x)^mom
Fmlba<- x^mlba*(2-x)^mlba
Fmomba<- x^momba*(2-x)^momba

plot(x,Fml)              # Change "Fml" to "Fcs" and to "Fmom" here & below
				 # And do the same thing to get the K-S test results based on the bootstrapped estimates

ii<- rep(1:n)

lower<- ii/n-Fml
upper<- Fml-(ii-1)/n

m<- pmax(lower,upper)
#lower
#upper
#m
Dn<- max(m)
Dn
e<- ii/n
plot(x, Fml,col="blue", type="l", xlab="Life-Span (/1000)",ylab="c.d.f.", main="Theoretical & Empirical c.d.f.'s")
lines(x,e, col="red", lty=3, lwd=2)
lines(x, Fcs, col="black", type="l", lty=2)
lines(x, Fmom, col="green", typ="l")
legend(0.25,0.25,legend=c("Theoretical (ML)", "Theoretical (C-S)", "Theoretical (MOM)","Empirical"),
       col=c("blue", "black", "green","red"), lty=c(1,2,1,3),cex=0.8,lwd=c(1,1,1,2))


