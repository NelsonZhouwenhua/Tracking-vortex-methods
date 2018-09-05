setwd("~/Desktop/MSc_Project_Vortex_Methods")
# install required packages:
install.packages("jpeg") # load image data
install.packages("imager") # blur image to improve boundary
install.packages("pracma") # Gauss-Newton method to fit circle
install.packages("astsa") # fit AR ARMA time series models
install.packages("rugarch") # fit ARCH GARCH time series models
install.packages("vars") # fit VAR VARMA time series models
install.packages("rmgarch") # fit MGARCH time series models

# library required packages
library(jpeg)
library(imager)
library(pracma)
library(astsa)
library(rugarch)
library(vars)
library(rmgarch)
library(grid)

# FUNCTIONS

# write function for image processing to isolated sandstorms only
imgPro <- function(k){
  for (i in 1:600){
    for (j in 1:800){
      if (abs(k[i,j,1]-1)>40/255 || abs(k[i,j,3]-1)>40/255 || abs(k[i,j,2]-0)>180/255){
        k[i,j,1] <- 1;
        k[i,j,2] <- 1;
        k[i,j,3] <- 1;
      }
    }
  }
  return(k)
}

# write function to classify image by crude thresholds 
imgCT <- function(k){
  a <- matrix(0,nrow=dim(k)[1],ncol=dim(k)[2])
  for (i in 1:600){
    for(j in 1:800){
      if (k[i,j,1]==1 && k[i,j,2]==1 && k[i,j,3]==1){
        a[i,j] = 1
      }
    }
  }
  return(a)
}

# naive fit circle
fitcircle <- function(k){
  z <- matrix(1,nrow=dim(k)[1],ncol=dim(k)[2])
  a <- sum(k==0)
  t <- which(k==0,arr.ind=TRUE)
  x <- median(t[,1])
  y <- median(t[,2])
  r <- round(mean(sqrt(a/pi)))
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt(abs(i-x)^2+abs(j-y)^2)<r){
        z[i,j] <- 0
      }
    }
  }
  return(z)
}

# find boundaries for regular shapes and irregular shapes
fb <- function(k){
  a <- matrix(1,nrow=dim(k)[1],ncol=dim(k)[2])
  a[,1] <- k[,1]
  a[1,] <- k[1,]
  a[,800] <- k[,800]
  a[600,] <- k[600,]
  for (i in 2:599){
    for(j in 2:799){
      if (k[i,j-1]==1 && k[i,j+1]==1 && k[i-1,j]==1 && k[i+1,j]==1 && k[i,j]==0){
        k[i,j] = 1
      }
      if (k[i,j-1]==1 && k[i,j]==0){
        a[i,j] = 0
      }
      if (k[i-1,j]==1 && k[i,j]==0){
        a[i,j] = 0
      }
      if (k[i,j+1]==1 && k[i,j]==0){
        a[i,j] = 0
      }
      if (k[i+1,j]==1 && k[i,j]==0){
        a[i,j] = 0
      }
    }
  }
  return(a)
}

fb1 <- function(k){
  a <- matrix(1,nrow=dim(k)[1],ncol=dim(k)[2])
  img <- as.cimg(t(k))
  img1 <- (isoblur(img,4)  < .5 )
  x <- matrix(1,nrow=dim(k)[1],ncol=dim(k)[2])
  t <- which(img1==FALSE,arr.ind = TRUE)
  for (i in 1:length(t[,1])){
    x[t[i,2],t[i,1]] <- 0
  }
  a[,1] <- k[,1]
  a[1,] <- k[1,]
  a[,800] <- k[,800]
  a[600,] <- k[600,]
  for (i in 2:599){
    for(j in 2:799){
      if (x[i,j-1]==1 && x[i,j]==0){
        a[i,j] = 0
      }
      if (x[i-1,j]==1 && x[i,j]==0){
        a[i,j] = 0
      }
      if (x[i,j+1]==1 && x[i,j]==0){
        a[i,j] = 0
      }
      if (x[i+1,j]==1 && x[i,j]==0){
        a[i,j] = 0
      }
    }
  }
  return(a)
}

# find the unique longest chain of boundary points
find_neighbours <- function(pool,current){
  within_x <- lapply(c(-1,0,1),function(r){return(which( pool[,1]==(current[1]+r) ))})
  within_y <- lapply(c(-1,0,1),function(r){return(which( pool[,2]==(current[2]+r) ))})
  neighbours_id <- c(intersect(within_x[[1]],within_y[[1]]),
                     intersect(within_x[[1]],within_y[[2]]),
                     intersect(within_x[[1]],within_y[[3]]),
                     intersect(within_x[[2]],within_y[[1]]), #(2,2) is `current`
                     intersect(within_x[[2]],within_y[[3]]),
                     intersect(within_x[[3]],within_y[[1]]),
                     intersect(within_x[[3]],within_y[[2]]),
                     intersect(within_x[[3]],within_y[[3]]))
  #return row numbers of the neignbours in the pool
  return(neighbours_id)
}

tiny2 <- function(pool){ #pool is a Nx2 matrix
  chain_list <- list()
  while(length(as.array(pool)) != 0){
    # initialize current chain and queue
    current_chain <- matrix(0,nrow(pool),2)
    current_chain_length <- 0
    current_queue <- matrix(0,nrow(pool),2)
    current_queue_length <- 0
    
    current_queue[1,] <- pool[1,]
    current_queue_length <- 1
    current_chain[1,] <- pool[1,]
    current_chain_length <- 1
    
    pool <- pool[-1,]
    if(length(as.array(pool))<=2){break}
    
    while(current_queue_length!=0){
      current <- current_queue[1,]
      
      current_queue <- current_queue[-1,] #remove from queue
      current_queue_length <- current_queue_length - 1
      
      # current_chain[current_chain_length+1,] <- current #add to chain
      # current_chain_length <- current_chain_length + 1
      
      neighbours_id <- find_neighbours(pool,current)
      
      if(length(neighbours_id)==0){next}
      old_queue_length <- current_queue_length #update queue length
      current_queue_length <- old_queue_length + length(neighbours_id)
      current_queue[(old_queue_length+1):current_queue_length,] <- 
        pool[neighbours_id,] #add to queue
      
      old_chain_length <- current_chain_length #update queue length
      current_chain_length <- old_chain_length + length(neighbours_id)
      current_chain[(old_chain_length+1):current_chain_length,] <- 
        pool[neighbours_id,] #add to chain
      
      pool <- pool[-neighbours_id,] #remove from pool
      if(length(as.array(pool))<=2){break}
    }
    current_chain <- current_chain[which(rowSums(current_chain)!=0),] #remove 0 rows
    chain_list <- c(chain_list,list(current_chain))
    if(length(as.array(pool))<=2){break}
  }
  chain_lengths <- sapply(1:length(chain_list),function(i){if (length(chain_list[[i]])>2){
    return( nrow(chain_list[[i]]) )
  }
    else{
      return(1)
    }  
  })
  return(chain_list[[which.max(chain_lengths)]])
}

# fit circles using boundary
fitSS <- function(xy,
                  a0=mean(xy[,1]),
                  b0=mean(xy[,2]),
                  r0 = mean(sqrt((xy[,1]-a0)^2 + (xy[,2]-b0)^2))){
  SS <- function(abr){
    (abr[3] - sqrt((xy[,1]-abr[1])^2 + (xy[,2]-abr[2])^2))^2
  }
  gaussNewton(c(a0,b0,r0), SS)
}

# fit ellipse using boundary
fit.ellipse <- function (x, y = NULL) {
  EPS <- 1.0e-8 
  dat <- xy.coords(x, y) 
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y) 
  D2 <- cbind(dat$x, dat$y, 1) 
  S1 <- t(D1) %*% D1 
  S2 <- t(D1) %*% D2 
  S3 <- t(D2) %*% D2 
  T <- -solve(S3) %*% t(S2) 
  M <- S1 + S2 %*% T 
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2) 
  evec <- eigen(M)$vec 
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2 
  a1 <- evec[, which(cond > 0)] 
  f <- c(a1, T %*% a1) 
  names(f) <- letters[1:6] 
  
  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b
  
  b2 <- f[2]^2 / 4
  
  center <- c(soln[1], soln[2]) 
  names(center) <- c("x", "y") 
  
  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6]) 
  den1 <- (b2 - f[1]*f[3]) 
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2) 
  den3 <- f[1] + f[3] 
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) )) 
  
  # calculate the angle of rotation 
  term <- (f[1] - f[3]) / f[2] 
  angle <- atan(1 / term) / 2 
  
  list(coef=f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle)) 
}

# fit GLM
imgGLM <- function(k){  
  f <- matrix(k,ncol = 1,nrow = 480000)
  i <- rep(seq(1,600),800)
  t <- rep(seq(1,800),each = 600)
  yx <- data.frame(cbind(f,i,t,i^2,t^2,i*t))#
  names(yx)[1]<- "y"
  fit <- glm(y~.,data = yx,family=binomial)
  new <- predict(fit,yx)
  fit$coefficients[1] <- fit$coefficients[1] - min(new)
  new <- predict(fit,yx)
  new <- range01(new)#
  new[new>sum(k==0)/480000] <- 1
  new[new<=sum(k==0)/480000] <- 0
  newk <- matrix(new,ncol = 800,nrow = 600)
  return(newk)
}

# fit a sequence of images in circles
cir <- function(a){
  x <- numeric(length(a))
  y <- numeric(length(a))
  r <- numeric(length(a))
  for (i in 1:length(a)){
    img <- readJPEG(a[i])
    img1 <- imgPro(img)
    img2 <- imgCT(img1)
    img3 <- fb1(img2)
    t <- which(img3==0,arr.ind=TRUE)
    if (length(t) == 0){
      next
    }
    k <- tiny2(t)
    fit <- fitSS(k)
    x[i] <- fit$xs[1]
    y[i] <- fit$xs[2]
    r[i] <- fit$xs[3]
  }
  return(cbind(x,y,r))
}

# fit a sequence of images in ellipses
ell <- function(k){
  x <- numeric(length(k))
  y <- numeric(length(k))
  maj <- numeric(length(k))
  min <- numeric(length(k))
  ang <- numeric(length(k))
  for (i in 1:length(k)){
    img <- readJPEG(k[i])
    img1 <- imgPro(img)
    img2 <- imgCT(img1)
    img3 <- fb1(img2)
    a <- which(img3==0,arr.ind=TRUE)
    if (length(a) == 0){
      next
    }
    b <- tiny2(a)
    imgell <- fit.ellipse(b)
    x[i] <- imgell$center[1]
    y[i] <- imgell$center[2]
    maj[i] <- imgell$major
    min[i] <- imgell$minor
    ang[i] <- imgell$angle
  }
  return(cbind(x,y,maj,min,ang))
}

# using glm fit image sequencen to ellipse
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
GLM <- function(k){
  x <- numeric(length(k))
  y <- numeric(length(k))
  maj <- numeric(length(k))
  min <- numeric(length(k))
  ang <- numeric(length(k))
  for (i in 1:length(k)){
    img <- readJPEG(k[i])
    img1 <- imgPro(img)
    img2 <- imgCT(img1)
    img3 <- imgGLM(img2)
    img4 <- fb(img3)
    t <- which(img4==0,arr.ind=TRUE)
    x[i] <- round(mean(t[,1]))
    y[i] <- round(mean(t[,2]))
    imgell <- fit.ellipse(t)
    maj[i] <- imgell$major
    min[i] <- imgell$minor
    ang[i] <- imgell$angle
  }
  return(cbind(x,y,maj,min,ang))
}


# fit images for the first and second sequence

# loi1 <- list.files(path="~/Desktop/MSc_Project_Vortex_Methods/sandstorm1",pattern="*.jpg")
# setwd("~/Desktop/MSc_Project_Vortex_Methods/sandstorm1")
# t <- Sys.time()
# CIR1 <- cir(loi1)
# Sys.time() -t
# t <- Sys.time()
# ELL1 <- ell(loi1)
# Sys.time() -t
# t <- Sys.time()
# GLM1 <- GLM(loi1)
# Sys.time() -t

# loi2 <- list.files(path="~/Desktop/MSc_Project_Vortex_Methods/sandstorm4",pattern="*.jpg")
# setwd("~/Desktop/MSc_Project_Vortex_Methods/sandstorm4")
# CIR2 <- cir(loi2[1:50])
# ELL2 <- ell(loi2[1:50])
# GLM2 <- GLM(loi2[1:50])
# CIR3 <- cir(loi2[97:146])
# ELL3 <- ell(loi2[97:146])
# GLM3 <- GLM(loi2[97:146])
# CIRm <- cir(loi2)
# ELLm <- ell(loi2)
# GLMm <- GLM(loi2)

# Autoregression time series model ARMA(p,q)

# fit ARMA for circles and predict
cirARMA <- function(k,p,q){
  train <- round(nrow(k)*0.8)
  cir1 <- arima(k[1:train,1], order = c(p, 0L, q),method = "CSS-ML")
  cir2 <- arima(k[1:train,2], order = c(p, 0L, q),method = "CSS-ML")
  cir3 <- arima(k[1:train,3], order = c(p, 0L, q),method = "CSS-ML")
  pred1 <- predict(cir1, n.ahead = nrow(k)-train)
  pred2 <- predict(cir2, n.ahead = nrow(k)-train)
  pred3 <- predict(cir3, n.ahead = nrow(k)-train)
  return(cbind(pred1$pred,pred2$pred,pred3$pred))
}
cirARMA1 <- function(k,p,q){
  train <- round(nrow(k)*0.8)
  cir1 <- arima(k[1:train,1], order = c(p, 0L, q),method = "ML")
  cir2 <- arima(k[1:train,2], order = c(p, 0L, q),method = "ML")
  cir3 <- arima(k[1:train,3], order = c(p, 0L, q),method = "ML")
  pred1 <- predict(cir1, n.ahead = nrow(k)-train)
  pred2 <- predict(cir2, n.ahead = nrow(k)-train)
  pred3 <- predict(cir3, n.ahead = nrow(k)-train)
  return(cbind(pred1$pred,pred2$pred,pred3$pred))
}
# fit ARMA for Ellipses and GLM and predict
ellARMA <- function(k,p,q){
  train <- round(nrow(k)*0.8)
  ell1 <- arima(k[1:train,1], order = c(p, 0, q),method = "CSS-ML")
  ell2 <- arima(k[1:train,2], order = c(p, 0, q),method = "CSS-ML")
  ell3 <- arima(k[1:train,3], order = c(p, 0, q),method = "CSS-ML")
  ell4 <- arima(k[1:train,4], order = c(p, 0, q),method = "CSS-ML")
  ell5 <- arima(k[1:train,5], order = c(p, 0, q),method = "CSS-ML")
  pred1 <- predict(ell1, n.ahead = nrow(k)-train)
  pred2 <- predict(ell2, n.ahead = nrow(k)-train)
  pred3 <- predict(ell3, n.ahead = nrow(k)-train)
  pred4 <- predict(ell4, n.ahead = nrow(k)-train)
  pred5 <- predict(ell5, n.ahead = nrow(k)-train)
  return(cbind(pred1$pred,pred2$pred,pred3$pred,pred4$pred,pred5$pred))
}
ellARMA1 <- function(k,p,q){
  train <- round(nrow(k)*0.8)
  ell1 <- arima(k[1:train,1], order = c(p, 0, q),method = "ML")
  ell2 <- arima(k[1:train,2], order = c(p, 0, q),method = "ML")
  ell3 <- arima(k[1:train,3], order = c(p, 0, q),method = "ML")
  ell4 <- arima(k[1:train,4], order = c(p, 0, q),method = "ML")
  ell5 <- arima(k[1:train,5], order = c(p, 0, q),method = "ML")
  pred1 <- predict(ell1, n.ahead = nrow(k)-train)
  pred2 <- predict(ell2, n.ahead = nrow(k)-train)
  pred3 <- predict(ell3, n.ahead = nrow(k)-train)
  pred4 <- predict(ell4, n.ahead = nrow(k)-train)
  pred5 <- predict(ell5, n.ahead = nrow(k)-train)
  return(cbind(pred1$pred,pred2$pred,pred3$pred,pred4$pred,pred5$pred))
}
# ARMA(1,1) predictions
# PREDcir1 <- cirARMA1(CIR1,1,1)
# PREDell1 <- ellARMA1(ELL1,1,1)
# PREDGLM1 <- ellARMA1(GLM1,1,1)
# AR(1) predictions
# PREDcir11 <- cirARMA1(CIRm,1,0)
# PREDell11 <- ellARMA1(ELLm,1,0)
# PREDGLM11 <- ellARMA1(GLMm,1,0)

# function help to plot 
cirCOM <- function(k,t,img){
  a <- matrix(1,nrow = 600,ncol = 800)
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((k[1]-i)^2+(k[2]-j)^2)<k[3]){
        a[i,j] <- 0
      }
    }
  }
  b1 <- fb(a)
  a <- matrix(1,nrow = 600,ncol = 800)
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((t[1]-i)^2+(t[2]-j)^2)<t[3]){
        a[i,j] <- 0
      }
    }
  }
  b2 <- fb(a)
  c <- img
  for (i in 1:length(which(b1==0,arr.ind = TRUE)[,1]))
    c[which(b1==0,arr.ind = TRUE)[i,1],which(b1==0,arr.ind = TRUE)[i,2],] <- c(0,0,0)
  for (i in 1:length(which(b2==0,arr.ind = TRUE)[,1]))
    c[which(b2==0,arr.ind = TRUE)[i,1],which(b2==0,arr.ind = TRUE)[i,2],] <- c(0,0,1)
  return(c)
}
cirCOM1 <- function(k,t,d,img){
  a <- matrix(1,nrow = 600,ncol = 800)
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((k[1]-i)^2+(k[2]-j)^2)<k[3]){
        a[i,j] <- 0
      }
    }
  }
  b1 <- fb(a)
  a <- matrix(1,nrow = 600,ncol = 800)
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((t[1]-i)^2+(t[2]-j)^2)<t[3]){
        a[i,j] <- 0
      }
    }
  }
  b2 <- fb(a)
  a <- matrix(1,nrow = 600,ncol = 800)
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((d[1]-i)^2+(d[2]-j)^2)<d[3]){
        a[i,j] <- 0
      }
    }
  }
  b3 <- fb(a)
  c <- img
  for (i in 1:length(which(b1==0,arr.ind = TRUE)[,1]))
    c[which(b1==0,arr.ind = TRUE)[i,1],which(b1==0,arr.ind = TRUE)[i,2],] <- c(0,0,0)
  for (i in 1:length(which(b2==0,arr.ind = TRUE)[,1]))
    c[which(b2==0,arr.ind = TRUE)[i,1],which(b2==0,arr.ind = TRUE)[i,2],] <- c(0,0,1)
  for (i in 1:length(which(b3==0,arr.ind = TRUE)[,1]))
    c[which(b3==0,arr.ind = TRUE)[i,1],which(b3==0,arr.ind = TRUE)[i,2],] <- c(1,0,0)
  return(c)
}

ellCOM <- function(k,t,img){
  a <- matrix(1,nrow = 600,ncol = 800)
  r <- sqrt(k[3]^2-k[4]^2)
  f1 <- c(k[1]- r*sin(k[5]),k[2]+r*cos(k[5])) 
  f2 <- c(k[1]+ r*sin(k[5]),k[2]-r*cos(k[5]))
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((i-f1[1])^2+(j-f1[2])^2)+sqrt((i-f2[1])^2+(j-f2[2])^2)<2*k[3]){
        a[i,j] <- 0
      }
    }
  }
  b1 <- fb(a)
  a <- matrix(1,nrow = 600,ncol = 800)
  r <- sqrt(t[3]^2-t[4]^2)
  f1 <- c(t[1]- r*sin(t[5]),t[2]+r*cos(t[5])) 
  f2 <- c(t[1]+ r*sin(t[5]),t[2]-r*cos(t[5]))
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((i-f1[1])^2+(j-f1[2])^2)+sqrt((i-f2[1])^2+(j-f2[2])^2)<2*t[3]){
        a[i,j] <- 0
      }
    }
  }
  b2 <- fb(a)
  c <- img
  for (i in 1:length(which(b1==0,arr.ind = TRUE)[,1]))
    c[which(b1==0,arr.ind = TRUE)[i,1],which(b1==0,arr.ind = TRUE)[i,2],] <- c(0,0,0)
  for (i in 1:length(which(b2==0,arr.ind = TRUE)[,1]))
    c[which(b2==0,arr.ind = TRUE)[i,1],which(b2==0,arr.ind = TRUE)[i,2],] <- c(0,0,1)
  return(c)
}
ellCOM1 <- function(k,t,d,img){
  a <- matrix(1,nrow = 600,ncol = 800)
  r <- sqrt(k[3]^2-k[4]^2)
  f1 <- c(k[1]- r*sin(k[5]),k[2]+r*cos(k[5])) 
  f2 <- c(k[1]+ r*sin(k[5]),k[2]-r*cos(k[5]))
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((i-f1[1])^2+(j-f1[2])^2)+sqrt((i-f2[1])^2+(j-f2[2])^2)<2*k[3]){
        a[i,j] <- 0
      }
    }
  }
  b1 <- fb(a)
  a <- matrix(1,nrow = 600,ncol = 800)
  r <- sqrt(t[3]^2-t[4]^2)
  f1 <- c(t[1]- r*sin(t[5]),t[2]+r*cos(t[5])) 
  f2 <- c(t[1]+ r*sin(t[5]),t[2]-r*cos(t[5]))
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((i-f1[1])^2+(j-f1[2])^2)+sqrt((i-f2[1])^2+(j-f2[2])^2)<2*t[3]){
        a[i,j] <- 0
      }
    }
  }
  b2 <- fb(a)
  a <- matrix(1,nrow = 600,ncol = 800)
  r <- sqrt(d[3]^2-d[4]^2)
  f1 <- c(d[1]- r*sin(d[5]),d[2]+r*cos(d[5])) 
  f2 <- c(d[1]+ r*sin(d[5]),d[2]-r*cos(d[5]))
  for (i in 1:600){
    for (j in 1:800){
      if (sqrt((i-f1[1])^2+(j-f1[2])^2)+sqrt((i-f2[1])^2+(j-f2[2])^2)<2*d[3]){
        a[i,j] <- 0
      }
    }
  }
  b3 <- fb(a)
  c <- img
  for (i in 1:length(which(b1==0,arr.ind = TRUE)[,1]))
    c[which(b1==0,arr.ind = TRUE)[i,1],which(b1==0,arr.ind = TRUE)[i,2],] <- c(0,0,0)
  for (i in 1:length(which(b2==0,arr.ind = TRUE)[,1]))
    c[which(b2==0,arr.ind = TRUE)[i,1],which(b2==0,arr.ind = TRUE)[i,2],] <- c(0,0,1)
  for (i in 1:length(which(b3==0,arr.ind = TRUE)[,1]))
    c[which(b3==0,arr.ind = TRUE)[i,1],which(b3==0,arr.ind = TRUE)[i,2],] <- c(1,0,0)
  return(c)
}
# img <- readJPEG(loi2[127])
# img <- readJPEG(loi2[142])
# img <- readJPEG(loi2[157])
# img <- imgPro(img)
# img1 <- imgPro(img)
# img1 <- imgCT(img1)
# img1 <- fb1(img1)
# k <- which(img1==0,arr.ind = TRUE)
# k <- tiny2(k)
# for (i in 1:length(k[,1]))
# img[k[i,1],k[i,2],] <- c(0,0,0)
# img1 <- cirCOM( CIRm[157,],PREDcir11[31,],img)
# img1 <- ellCOM( ELLm[157,],PREDell11[31,],img)
# img1 <- ellCOM( GLMm[157,],PREDGLM11[31,],img)
# img1 <- cirCOM( CIRm[157,],cirvar[31,],img)
# img1 <- ellCOM( ELLm[127,],ellvar[1,],img)
# img1 <- ellCOM( GLMm[157,],glmvar[31,],img)
# img1 <- cirCOM1( CIRm[157,],PREDcir11[31,],cirvar[31,],img)
# img1 <- ellCOM1( ELLm[157,],PREDell11[31,],ellvar[31,],img)
# img1 <- ellCOM1( GLMm[157,],PREDGLM11[31,],glmvar[31,],img)
# img1 <- ellCOM1( ELLm[127,],ellvar[1,],DCCPell[1,],img)
# img1 <- ellCOM1( ELLm[157,],ellvar[31,],DCCPell[31,],img)
# img1 <- ellCOM1( ELLm[157,],DCCPell[31,],CGARCHSell,img)
# img1 <- ellCOM1( ELLm[127,],DCCPell[1,],CGARCHSell,img)
# a <- rgb(img1[,,1],img1[,,2],img1[,,3])
# dim(a) <- dim(img1[,,1])
# grid.raster(a)

# ARMA rmse
cirRMSE <- function(k,t){
  test <- round(nrow(t)*0.2)
  k <- unclass(k)
  v1 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),1]-k[,1])^2))
  v2 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),2]-k[,2])^2))
  v3 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),3]-k[,3])^2))
  z <- cbind(v1,v2,v3)
  colnames(z) <- c("x","y","r")
  return(z)
}

ellRMSE <- function(k,t){
  test <- round(nrow(t)*0.2)
  k <- unclass(k)
  v1 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),1]-k[,1])^2))
  v2 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),2]-k[,2])^2))
  v3 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),3]-k[,3])^2))
  v4 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),4]-k[,4])^2))
  v5 <- sqrt(mean((t[(nrow(t)-test+1):nrow(t),5]-k[,5])^2))
  z <- cbind(v1,v2,v3,v4,v5)
  colnames(z) <- c("x","y","major","minor","angle")
  return(z)
}

# ARMAcirRMSE <- cirRMSE(PREDcir1,CIR1)
# ARMAellRMSE <- ellRMSE(PREDell1,ELL1)
# ARMAGLMRMSE <- ellRMSE(PREDGLM1,GLM1)
# ARcirRMSE <- cirRMSE(PREDcir11,CIRm)
# ARellRMSE <- ellRMSE(PREDell11,ELLm)
# ARGLMRMSE <- ellRMSE(PREDGLM11,GLMm)

# ARMA MAPE
cirMAPE <- function(k,t){
  test <- round(nrow(t)*0.2)
  k <- unclass(k)
  v1 <- mean(abs((t[(nrow(t)-test+1):nrow(t),1]-k[,1])/t[(nrow(t)-test+1):nrow(t),1]))
  v2 <- mean(abs((t[(nrow(t)-test+1):nrow(t),2]-k[,2])/t[(nrow(t)-test+1):nrow(t),2]))
  v3 <- mean(abs((t[(nrow(t)-test+1):nrow(t),3]-k[,3])/t[(nrow(t)-test+1):nrow(t),3]))
  z <- cbind(v1,v2,v3)
  colnames(z) <- c("x","y","r")
  return(z)
}

ellMAPE <- function(k,t){
  test <- round(nrow(t)*0.2)
  k <- unclass(k)
  v1 <- mean(abs((t[(nrow(t)-test+1):nrow(t),1]-k[,1])/t[(nrow(t)-test+1):nrow(t),1]))
  v2 <- mean(abs((t[(nrow(t)-test+1):nrow(t),2]-k[,2])/t[(nrow(t)-test+1):nrow(t),2]))
  v3 <- mean(abs((t[(nrow(t)-test+1):nrow(t),3]-k[,3])/t[(nrow(t)-test+1):nrow(t),3]))
  v4 <- mean(abs((t[(nrow(t)-test+1):nrow(t),4]-k[,4])/t[(nrow(t)-test+1):nrow(t),4]))
  v5 <- mean(abs((t[(nrow(t)-test+1):nrow(t),5]-k[,5])/t[(nrow(t)-test+1):nrow(t),5]))
  z <- cbind(v1,v2,v3,v4,v5)
  colnames(z) <- c("x","y","major","minor","angle")
  return(z)
}

# ARMAcirMAPE <- cirMAPE(PREDcir1,CIR1)
# ARMAellMAPE <- ellMAPE(PREDell1,ELL1)
# ARMAGLMMAPE <- ellMAPE(PREDGLM1,GLM1)

# ARCH(p) and GARCH(p,q) using ARMA(1,1)
cirGARCH <- function(k,p,q){
  spec <- ugarchspec(variance.model = list(garchOrder = c(p, q)),mean.model=list(
    armaOrder=c(1,1)))
  train <- round(nrow(k)*0.8)
  f1 <- ugarchfit(spec, k[1:train,1])
  f2 <- ugarchfit(spec, k[1:train,2])
  f3 <- ugarchfit(spec, k[1:train,3])
  cir1 <- ugarchforecast(f1)
  cir2 <- ugarchforecast(f2)
  cir3 <- ugarchforecast(f3)
  return(cbind(fitted(cir1),fitted(cir2),fitted(cir3)))
}

ellGARCH <- function(k,p,q){
  spec <- ugarchspec(variance.model = list(garchOrder = c(p, q)),mean.model=list(
    armaOrder=c(1,1)))
  train <- round(nrow(k)*0.8)
  f1 <- ugarchfit(spec, k[1:train,1])
  f2 <- ugarchfit(spec, k[1:train,2])
  f3 <- ugarchfit(spec, k[1:train,3])
  f4 <- ugarchfit(spec, k[1:train,4])
  f5 <- ugarchfit(spec, k[1:train,5])
  ell1 <- ugarchforecast(f1)
  ell2 <- ugarchforecast(f2)
  ell3 <- ugarchforecast(f3)
  ell4 <- ugarchforecast(f4)
  ell5 <- ugarchforecast(f5)
  return(cbind(fitted(ell1),fitted(ell2),fitted(ell3),fitted(ell4),fitted(ell5)))
}

# GARCH(1,1) predictions
# cirGARCH112 <- cirGARCH(CIR2,1,1)
# ellGARCH112 <- ellGARCH(ELL2,1,1)
# GLMGARCH112 <- ellGARCH(GLM2,1,1)
# cirGARCH113 <- cirGARCH(CIR3,1,1)
# ellGARCH113 <- ellGARCH(ELL3,1,1)
# GLMGARCH113 <- ellGARCH(GLM3,1,1)

# ARCH(1) predictions
# cirGARCH12 <- cirGARCH(CIR2,1,0)
# ellGARCH12 <- ellGARCH(ELL2,1,0)
# GLMGARCH12 <- ellGARCH(GLM2,1,0)
# cirGARCH13 <- cirGARCH(CIR3,1,0)
# ellGARCH13 <- ellGARCH(ELL3,1,0)
# GLMGARCH13 <- ellGARCH(GLM3,1,0)

# img <- readJPEG(loi2[41])
# img <- imgPro(img)
# img1 <- cirCOM( CIR2[41,],cirGARCH112[1,],img)
# img1 <- ellCOM( ELL2[41,],ellGARCH112[1,],img)
# img1 <- ellCOM( GLM2[41,],GLMGARCH112[1,],img)
# img1 <- cirCOM( CIR2[41,],cirGARCH12[1,],img)
# img1 <- ellCOM( ELL2[41,],ellGARCH12[1,],img)
# img1 <- ellCOM( GLM2[41,],GLMGARCH12[1,],img)
# a <- rgb(img1[,,1],img1[,,2],img1[,,3])
# dim(a) <- dim(img1[,,1])
# grid.raster(a)

# GARCHcirRMSE <- cirRMSE(cirGARCH112,CIR2)
# GARCHellRMSE <- ellRMSE(ellGARCH112,ELL2)
# GARCHGLMRMSE <- ellRMSE(GLMGARCH112,GLM2)

# ARCHcirRMSE <- cirRMSE(cirGARCH12,CIR2)
# ARCHellRMSE <- ellRMSE(ellGARCH12,ELL2)
# ARCHGLMRMSE <- ellRMSE(GLMGARCH12,GLM2)

# ARMAcirRMSE4 <- cirRMSE(PREDcir4,CIR4)
# ARMAellRMSE4 <- ellRMSE(PREDell4,ELL4)
# ARMAGLMRMSE4 <- ellRMSE(PREDGLM4,GLM4)

# VAR circles and ellipses
var1 <- VAR(CIRm[1:126,])
prd1 <- predict(var1,n.ahead = 31)
cirvar <- cbind(prd1$fcst$x[,1],prd1$fcst$y[,1],prd1$fcst$r[,1])
var2 <- VAR(ELLm[1:126,])
prd2 <- predict(var2,n.ahead = 31, ci = 0.95, dumvar = NULL)
ellvar <- cbind(prd2$fcst$x[,1],prd2$fcst$y[,1],prd2$fcst$maj[,1],prd2$fcst$min[,1],prd2$fcst$ang[,1])
var3 <- VAR(GLMm[1:126,])
prd3 <- predict(var3,n.ahead = 31, ci = 0.95, dumvar = NULL)
glmvar <- cbind(prd3$fcst$x[,1],prd3$fcst$y[,1],prd3$fcst$maj[,1],prd3$fcst$min[,1],prd3$fcst$ang[,1])


# DCC-GARCH and COPULA-GARCH circles
xspec = ugarchspec(mean.model = list(armaOrder = c(1, 1)), variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), distribution.model = 'norm')
uspec = multispec(replicate(3, xspec))
dspec = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
cspec = cgarchspec(uspec = uspec, dccOrder = c(1, 1))
cgarchcir <- cgarchfit(cspec,CIRm[1:126,])
dcccir <- dccfit(dspec,CIRm[1:126,])
DCCPC <-dccforecast(dcccir, n.ahead = 31)
DCCPcir <- cbind(fitted(DCCPC)[,1,],fitted(DCCPC)[,2,],fitted(DCCPC)[,3,])
CGARCHSC <- cgarchsim(cgarchcir, n.sim = 1000,m.sim = 31)###
CGARCHScir <- c(mean(fitted(CGARCHSC)[101:1000,1]),mean(fitted(CGARCHSC)[101:1000,2]),mean(fitted(CGARCHSC)[101:1000,3]))


# DCC-GARCH and COPULA-GARCH ellipses
uspec = multispec(replicate(5, xspec))
dspec = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
cspec = cgarchspec(uspec = uspec, dccOrder = c(1, 1))
cgarchell <- cgarchfit(cspec,ELLm[1:126,])
dccell <- dccfit(dspec,ELLm[1:126,])
DCCPE <-dccforecast(dccell, n.ahead = 31)
DCCPell <- cbind(fitted(DCCPE)[,1,],fitted(DCCPE)[,2,],fitted(DCCPE)[,3,],fitted(DCCPE)[,4,],fitted(DCCPE)[,5,])
CGARCHSE <- cgarchsim(cgarchell, n.sim = 1000,m.sim = 1)###
CGARCHSell <- c(mean(fitted(CGARCHSE)[101:1000,1]),mean(fitted(CGARCHSE)[101:1000,2]),mean(fitted(CGARCHSE)[101:1000,3]),mean(fitted(CGARCHSE)[101:1000,4]),mean(fitted(CGARCHSE)[101:1000,5]))

dccglm <- dccfit(dspec,GLMm[1:126,])
DCCPG <-dccforecast(dccglm, n.ahead = 31)
DCCPGLM <- cbind(fitted(DCCPG)[,1,],fitted(DCCPG)[,2,],fitted(DCCPG)[,3,],fitted(DCCPG)[,4,],fitted(DCCPG)[,5,])



