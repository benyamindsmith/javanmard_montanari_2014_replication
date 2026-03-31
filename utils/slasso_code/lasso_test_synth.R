library(Matrix);
library(glmnet);
library(expm);
library(flare);
source("lasso_inference.R");

p <- 200;
n <- 150;
s0 <- 15;
b <- 1;
b0 <- 0;
sigma <- 0.5;
set.seed('1')

X <- rbinom(p*n,1,prob=0.15);
dim(X) <- c(n,p);
X <- X %*% diag(1+9*runif(p))

theta0 <- c(rep(b,s0),rep(0,p-s0));
w <- sigma*rnorm(n);
y <- (b0+X%*%theta0+w);

th <- SSLasso(X,y,verbose = TRUE)
dev.new()
plot(th$coef,ylim=c(-1.5*b,1.5*b), main='Confidence Intervals based on de-biased LASSO', ylab='', xlab = 'Coefficients');
points(th$unb.coef,col="blue");
points(theta0,col="green");
lines(th$up.lim,col="red");
lines(th$low.lim,col="red");
legend('topright', legend=c('LASSO','de-biased LASSO','Ground-truth','Confidence Intervals'), col=c('black', 'blue','green','red'), pch=c(1,1,1,NA_integer_), lty = c(0,0,0,1))


print(paste("Estimated Noise Sdev = ",th$noise.sd));
print(paste("Estimated Norm0 = ",th$norm0));
nc <- sum(as.numeric(th$low.lim>theta0)+as.numeric(th$up.lim<theta0));
print(paste("Coverage = ",(p-nc)/p));
print("P-values:")
print(th$pvals)

