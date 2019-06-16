
skaters <- as.data.frame(matrix(scan("fivektenk-data",skip=5),ncol=8,byrow=T))
skate = skaters
skate$fivekmin = skaters[,3]*60+skate[,4]+skate[,5]/100
skate$tenkmin = skaters[,6]*60+skate[,7]+skate[,8]/100
skate = skate[,-c(3,4,5,6,7,8)]
skate$yy = skate$tenkmin/skate$fivekmin

year <- skaters[ ,1]
fivek <- 60*skaters[ ,3] + skaters[ ,4] + skaters[ ,5]/100
tenk <- 60*skaters[ ,6] + skaters[ ,7] + skaters[ ,8]/100
yy <- tenk/fivek
events <- 2000:2019
mm <- 0*events
mean5k <- 0*events
upper5k <- 0*events
lower5k <- 0*events
mean10k <- 0*events
upper10k <- 0*events
lower10k <- 0*events
meanyy <- 0*events
upperyy <- 0*events
loweryy <- 0*events

yyear = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,
          2016,2017,2018,2019)

kk <- length(events)



for (j in 1:kk) {
  now5k <- fivek[year==events[j]]
  now10k <- tenk[year==events[j]]
  nowyy <- now10k/now5k
  
  mean5k[j] <- mean(now5k)
  mean10k[j] <- mean(now10k)
  meanyy[j] <- mean(nowyy)
  
  
  upper10k[j] <- mean(now10k)+1.96*sd(now10k)
  lower10k[j] <- mean(now10k)-1.96*sd(now10k)
  upperyy[j] <- mean(nowyy)+1.96*sd(nowyy)
  loweryy[j] <- mean(nowyy)-1.96*sd(nowyy)
  
  mm[j] <- length(now5k)
  }


plot(skate$V1, skate$tenkmin,ylim=c(725,925), xlab="Verdensmesterskap", ylab="10000 m")
#

#plot(events,mean5k,ylim=c(350,450))
matlines(events,mean10k)
matlines(events,upper10k, col="red")
matlines(events,lower10k, col="red")


plot(skate$V1, skate$tenkmin/skate$fivekmin,ylim=c(2,2.2), xlab="Verdensmesterskap", ylab="ratio 10k/5k")
#

#plot(events,mean5k,ylim=c(350,450))
matlines(events,meanyy)
matlines(events,upperyy, col="red")
matlines(events,loweryy, col="red")


mu6=ML6[21:40]
muse6=se6[21:40]
lowermu6=mu6-1.645*muse6
uppermu6=mu6+1.645*muse6

plot(yyear, mu6,ylim=c(2,2.2), xlab="Verdensmesterskap", ylab="xi")
#

#plot(events,mean5k,ylim=c(350,450))
matlines(events,mu6)
matlines(events,lowermu6, col="red")
matlines(events,uppermu6, col="red")

sigma6=ML6[1:20]
sigmase6=se6[1:20]
lowersigma6=sigma6-1.645*sigmase6
uppersigma6=sigma6+1.645*sigmase6

plot(yyear, sigma6,ylim=c(0.006,0.06), xlab="Verdensmesterskap", ylab="sigma")
#

#plot(events,mean5k,ylim=c(350,450))
matlines(events,sigma6)
matlines(events,lowersigma6, col="red")
matlines(events,uppersigma6, col="red")


# for each year

y2019 = skate[skate$V1 == 2019,]$yy
y2018 = skate[skate$V1 == 2018,]$yy
y2017 = skate[skate$V1 == 2017,]$yy
y2016 = skate[skate$V1 == 2016,]$yy
y2015 = skate[skate$V1 == 2015,]$yy
y2014 = skate[skate$V1 == 2014,]$yy
y2013 = skate[skate$V1 == 2013,]$yy
y2012 = skate[skate$V1 == 2012,]$yy
y2011 = skate[skate$V1 == 2011,]$yy
y2010 = skate[skate$V1 == 2010,]$yy
y2009 = skate[skate$V1 == 2009,]$yy
y2008 = skate[skate$V1 == 2008,]$yy
y2007 = skate[skate$V1 == 2007,]$yy
y2006 = skate[skate$V1 == 2006,]$yy
y2005 = skate[skate$V1 == 2005,]$yy
y2004 = skate[skate$V1 == 2004,]$yy
y2003 = skate[skate$V1 == 2003,]$yy
y2002 = skate[skate$V1 == 2002,]$yy
y2001 = skate[skate$V1 == 2001,]$yy
y2000 = skate[skate$V1 == 2000,]$yy



yyy = list(y2000, y2001, y2002, y2003, y2004, y2005, y2006, y2007, y2008, y2009,
           y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019)


  # MODEL 1
  
minuslogL1 <- function(para)
{
  mu <- para[1]
  sigma <- para[2]
    
  sum ( 0.5*(yy - mu)^2/sigma^2 + log(sigma) )
}
logL1 <- function(para) {-minuslogL1(para)}

nils1 <- nlm(minuslogL1,c(2,1),hessian=T)

ML1 <- nils1$estimate
Jhat1 <- nils1$hessian
se1 <- sqrt(diag(solve(Jhat1)))
showme1a <- cbind(ML1,se1)
print(round(showme1a,4))

aic1 <- 2*logL1(ML1) -2*2
aic1





  # MODEL 2
  
  
minuslogL2 <- function(para) {
  
  out =0
  for (j in 1:length(yyy)) {
    sigma <- para[1]
    y = unlist(yyy[j], use.names=FALSE)
      
    out = out+ sum( 0.5*(y - para[j+1])^2/sigma^2 + log(sigma) )
    
  }
  return(out)
}
logL2 <- function(para) {-minuslogL2(para)}
  
nils2 <- nlm(minuslogL2,numeric(21)+1,hessian=T)

ML2 <- nils2$estimate
Jhat2 <- nils2$hessian
se2 <- sqrt(diag(solve(Jhat2)))
showme2a <- cbind(ML2,se2)
print(round(showme2a,4))

aic2 <- 2*logL2(ML2) -2*21
aic2



# MODEL 3


minuslogL3 <- function(para) {
  
  out =0
  for (j in 1:length(yyy)) {
    mu = para[1]
    y = unlist(yyy[j], use.names=FALSE)
    
    out = out+ sum( 0.5*(y - mu)^2/para[j+1]^2 + log(para[j+1]) )
    
  }
  return(out)
}
logL3 <- function(para) {-minuslogL3(para)}

nils3 <- nlm(minuslogL3,numeric(21)+1,hessian=T)

ML3 <- nils3$estimate
Jhat3 <- nils3$hessian
se3 <- sqrt(diag(solve(Jhat3)))
showme3a <- cbind(ML3,se3)
print(round(showme3a,4))

aic3 <- 2*logL3(ML3) -2*21
aic3


# MODEL 4

minuslogL4 <- function(para) {
  
  out =0
  sigma <- para[1]
  a = para[2]
  b = para[3]
  
  
  for (j in 1:length(yyy)) {
    
    y = unlist(yyy[j], use.names=FALSE)
    w = (j-9)/5.627
    mu = a+w*b
    
    out = out+ sum( 0.5*(y - mu)^2/sigma^2 + log(sigma) )
    
  }
  return(out)
}
logL4 <- function(para) {-minuslogL4(para)}

nils4 <- nlm(minuslogL4,c(1,1,1),hessian=T)

ML4 <- nils4$estimate
Jhat4 <- nils4$hessian
se4 <- sqrt(diag(solve(Jhat4)))
showme4a <- cbind(ML4,se4)
print(round(showme4a,4))

aic4 <- 2*logL4(ML4) -2*3
aic4


# MODEL 5

minuslogL5 <- function(para) {
  
  out =0
  sigma0 <- para[1]
  a = para[2]
  b = para[3]
  g = para[4]
  
  for (j in 1:length(yyy)) {
    
    y = unlist(yyy[j], use.names=FALSE)
    w = (j-9)/5.627
    mu = a+w*b
    sigma = sigma0*exp(g*w)
    
    out = out+ sum( 0.5*(y - mu)^2/sigma^2 + log(sigma) )
    
  }
  return(out)
}
logL5 <- function(para) {-minuslogL5(para)}

nils5 <- nlm(minuslogL5,c(1,1,1,1),hessian=T)

ML5 <- nils5$estimate
Jhat5 <- nils5$hessian
se5 <- sqrt(diag(solve(Jhat5)))
showme5a <- cbind(ML5,se5)
print(round(showme5a,4))

aic5 <- 2*logL5(ML5) -2*4
aic5



# MODEL 6


minuslogL6 <- function(para) {
  
  out =0
  for (j in 1:length(yyy)) {
    y = unlist(yyy[j], use.names=FALSE)
    
    out = out+ sum( 0.5*(y - para[j+20])^2/para[j]^2 + log(para[j]) )
    
  }
  return(out)
}
logL6 <- function(para) {-minuslogL6(para)}

nils6 <- nlm(minuslogL6,numeric(40)+1,hessian=T)

ML6 <- nils6$estimate
Jhat6 <- nils6$hessian
se6 <- sqrt(diag(solve(Jhat6)))
showme6a <- cbind(ML6,se6)
print(round(showme6a,4))

aic6 <- 2*logL6(ML6) -2*40
aic6


# MULTINORMAL



minuslogL7 <- function(para) {
  mu = para[1]
  tau = para[2]
  sigma = para[3]
  
  out =0
  for (j in 1:length(yyy)) {
    y = unlist(yyy[j], use.names=FALSE)
    m = length(y)
    sig = sigmaMatrix(c(m, tau, sigma))
    
    out = out+ m/2*log(2*pi)+1/2*log(determinant(sig, logarithm = FALSE)$modulus)+1/2*t(y-mu)%*%solve(sig)%*%(y-mu)
    
    
  }
  return(as.double(out))
}
logL7 <- function(para) {-minuslogL7(para)}

nils7 <- nlm(minuslogL7,c(2.1,0.01,0.02),hessian=T)

ML7 <- nils7$estimate
Jhat7 <- nils7$hessian
se7 <- sqrt(diag(solve(Jhat7)))
showme7a <- cbind(ML7,se7)
print(round(showme7a,4))

aic7 <- 2*logL7(ML7) -2*3
aic7


sigmaMatrix <- function(para) {
  m = para[1]
  tau = para[2]
  sigma = para[3]
  
  s = matrix(0,m,m)
  
  for (i in 1:m) {
    for (j in 1:m) {
      
      if (i == j) {
        s[i,j] = tau^2 + sigma^2
      }
      else {
        s[i,j] = tau^2
      }
      
    }
  }
  s
}





