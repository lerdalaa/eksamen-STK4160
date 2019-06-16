
accidents <- as.data.frame(matrix(scan("sweden-accidents",skip=11),ncol=5,byrow=T))
colnames(accidents) <- c("year","year195x","cars","gas", "y")

accidents$y_small = accidents$y/100

y = accidents[,5]
y_small = accidents[,6]
x1 = accidents[,2]
x2 = accidents[,3]
x3 = accidents[,4]
accidents$cars_norm = (x2-mean(x2))/sd(x2)

fit1 <- glm(y ~ year195x + cars + gas, family=poisson, data=accidents)
fit1_small <- glm(y_small ~ year195x + cars + gas, family=poisson, data=accidents)

fit_lm <- lm(y ~ year195x + cars + gas, data=accidents)

summary( fit1 )



plot(accidents$year,y,xlab="år", ylab="dødsfall i trafikken")
lines(accidents$year,fit1$fitted.values)

residuals = (y-fit1$fitted.values)/sqrt(fit1$fitted.values)

plot(fitted(fit1),residuals, xlab="y values", ylab="standardized residuals")
plot(x1,residuals, xlab="�r etter 1950", ylab="standardized residuals")

plot(x2,residuals)
plot(x1,residuals)
plot(x1,fit1$residuals)


one = numeric(56)+1
scorem1 = cbind(one,x1,x2,x3)

JM1 = t(scorem1)%*%scorem1





lambda = exp(fit1$coefficients[1]+x1*fit1$coefficients[2]+x2*fit1$coefficients[3]
             +x3*fit1$coefficients[4])
res_Kmatrix = matrix(0,4,4)
res_Jmatrix = matrix(0,4,4)
for (i in 1:56) {
  
  xmatrix = cbind(c(1, x1[i],x2[i],x3[i]),
  c(x1[i],x1[i]*x1[i],x1[i]*x2[i],x1[i]*x3[i]),
  c(x2[i],x2[i]*x1[i],x2[i]*x2[i],x2[i]*x3[i]),
  c(x3[i],x3[i]*x1[i],x3[i]*x2[i],x3[i]*x3[i]))
  
  Kmatrix = (y[i]-lambda[i])^2*xmatrix
  Jmatrix = lambda[i]*xmatrix
  
  res_Kmatrix = res_Kmatrix + Kmatrix
  res_Jmatrix = res_Jmatrix + Jmatrix
}
Kmatrix = res_Kmatrix/56
Jmatrix = res_Jmatrix/56

Jinverse = solve(Jmatrix)
sum(diag(Jinverse*Kmatrix))



fit2 <- glm(y ~ year195x + cars + gas + I(cars_norm^2), family=poisson, data=accidents)
fit3 <- glm(y ~ year195x + cars + gas + I(cars_norm^2) + I(cars_norm^3), family=poisson, data=accidents)
fit4 <- glm(y ~ year195x + cars + gas + I(cars_norm^2) + I(cars_norm^3) + I(cars_norm^4), family=poisson, data=accidents)

sum1 =0
sum2 =0
sum3 =0
sum4 =0
for (i in 1:56) {
  push = accidents[c(i),]
  
  temp <- accidents[-c(i),]
  temp1 <- glm(y ~ year195x + cars + gas, family=poisson, data=temp)
  temp2 <- glm(y ~ year195x + cars + gas, I(cars_norm^2), family=poisson, data=temp)
  temp3 <- glm(y ~ year195x + cars + gas + I(cars_norm^2) + I(cars_norm^3), family=poisson, data=temp)
  temp4 <- glm(y ~ year195x + cars + gas + I(cars_norm^2) + I(cars_norm^3) + I(cars_norm^4), family=poisson, data=temp)
  
  # predict for each
  sum1 = sum1 + abs(exp(predict(temp1, newdata=push))-push$y)
  sum2 = sum2 + abs(exp(predict(temp2, newdata=push))-push$y)
  sum3 = sum3 + abs(exp(predict(temp3, newdata=push))-push$y)
  sum4 = sum4 + abs(exp(predict(temp4, newdata=push))-push$y)

}
sum1=sum1/56
sum2=sum2/56
sum3=sum3/56
sum4=sum4/56

sum21 =0
sum22 =0
sum23 =0
sum24 =0
for (i in 1:56) {
  push = accidents[c(i),]
  
  temp <- accidents[-c(i),]
  temp1 <- glm(y ~ year195x + cars + gas, family=poisson, data=temp)
  temp2 <- glm(y ~ year195x + cars + gas, I(cars_norm^2), family=poisson, data=temp)
  temp3 <- glm(y ~ year195x + cars + gas + I(cars_norm^2) + I(cars_norm^3), family=poisson, data=temp)
  temp4 <- glm(y ~ year195x + cars + gas + I(cars_norm^2) + I(cars_norm^3) + I(cars_norm^4), family=poisson, data=temp)
  
  lambda1I = exp(temp1$coefficients[1]+x1[i]*temp1$coefficients[2]
                 +x2[i]*temp1$coefficients[3]+x3[i]*temp1$coefficients[4])
  lambda2I = exp(temp2$coefficients[1]+x1[i]*temp1$coefficients[2]
                 +x2[i]*temp1$coefficients[3]+x3[i]*temp1$coefficients[4])
  lambda3I = exp(temp3$coefficients[1]+x1[i]*temp1$coefficients[2]
                 +x2[i]*temp1$coefficients[3]+x3[i]*temp1$coefficients[4])
  lambda4I = exp(temp4$coefficients[1]+x1[i]*temp1$coefficients[2]
                 +x2[i]*temp1$coefficients[3]+x3[i]*temp1$coefficients[4])
  # predict for each
  sum21 = sum21 + log(dpois(y[i],lambda1I))
  sum22 = sum22 + log(dpois(y[i],lambda2I))
  sum23 = sum23 + log(dpois(y[i],lambda3I))
  sum24 = sum24 + log(dpois(y[i],lambda4I))
  
}
sum21 = sum21/56
sum22 = sum22/56
sum23 = sum23/56
sum24 = sum24/56


# predict y_2011
pred = accidents[c(i),]
pred$year195x = 61
exp(predict(fit1, newdata=pred))
exp(predict(fit2, newdata=pred))
exp(predict(fit3, newdata=pred))
exp(predict(fit4, newdata=pred))





logLG <- function(para) {
  beta0 = para[1]
  beta1 = para[2]
  beta2 = para[3]
  beta3 = para[4]
  c = para[5]
  
  #lambdaG = dgamma(y_small,c*exp(beta0+x1*beta1+x2*beta2+x3*beta3),c)
  #sum(-lambdaG +y_small*log(lambdaG)-factorial(y_small))
  
  sum = 0
  for (i in 1:56) {
    lambdaG = dgamma(y_small[i],c*exp(beta0+x1[i]*beta1+x2[i]*beta2+x3[i]*beta3),c)
    sum = sum -lambdaG +y_small[i]*log(lambdaG)-log(factorial(y_small[i]))
  }
  sum
  
}

logLG <- function(para) {
  beta0 = para[1]
  beta1 = para[2]
  beta2 = para[3]
  beta3 = para[4]
  c = para[5]
  
  #lambdaG = dgamma(y_small,c*exp(beta0+x1*beta1+x2*beta2+x3*beta3),c)
  #sum(-lambdaG +y_small*log(lambdaG)-factorial(y_small))
  
  sum = 0
  for (i in 1:56) {
    lambdaG = dgamma(y_small[i],c*exp(beta0+x1[i]*beta1+x2[i]*beta2+x3[i]*beta3),c)
    sum = sum -1/2*log(2*pi)-1/2*log(lambdaG)-1/(2*lambdaG^2)*(y[i]-lambdaG)^2
  }
  sum
  
}

minuslogLG <- function(para) {-logLG(para)}

nilsG <- nlm(minuslogLG,c(2.01,-0.086,0.93,-0.0198,1),hessian=T)

MLG <- nilsG$estimate
JhatG <- nilsG$hessian
seG <- sqrt(diag(solve(JhatG)))
showmeGa <- cbind(MLG,seG)
print(round(showmeGa,4))

aicG <- 2*logLG(MLG) -2*3
aicG




