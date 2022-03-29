library("rmutil")

print("DONT MIX UP SD AND VARIANCE. YOU ARE WARNED")
samp_z <- function(n=1, h = 0.2, dof = 1, dist = "normal")
{
   ret <- NULL
   if(is.matrix(h) && dim(h)[1] != 1 && dim(h)[2] != 1)
   {
       print(h)
       stop("h is matrix, samp_z is written for scalar and vector h. take care. might produce wrong samples.")
   }
   if(dist == "normal")
   {
       ret <- rnorm(n, sd = h)
   }
   else if (dist == "bim_normal") {
      if(any(h > 1))
      {
          print(h)
          stop("choose h less than 1. mean is not defined.")
      }

      ret <- rnorm(n, mean = sqrt(1 - h^2), sd = h)
      # invert it to generate from bimodal distribution, possibly for each target.
      invert_or_not <- (runif(n) < 0.5)
      weight <- 2 * invert_or_not - 1
      ret <- ret * weight
   }
   else if (dist == "laplace") {
       ret <- rlaplace(n, s = h)
   }
   else if (is.null(ret))
   {
      stop("something's wrong, I can feel it. ret is null")
   }
   return(ret)
}

log_proposal_dist <- function(x, mean = 0, h = 1, dof = 1, dist = "normal")
{
    if(dist == "normal")
    {
        ret <- dnorm(x, mean = mean, sd = h, log = TRUE)
    }
    else if(dist == "bim_normal")
    {
        if(h > 1)
        {
            print("h > 1, choose correct h.")
            quit()
        }
        ret <- log(0.5 * ( dnorm(x, mean = mean - sqrt(1 - h^2), sd = h) + dnorm(x, mean = mean + sqrt(1 - h^2), sd = h) ))
    }
    else if (dist == "laplace") {
        ret <- dlaplace(x, mean = mean, s = h, log = TRUE)
    }
    return(ret)
}

ksd <- function(x, gradlogp, c = 1, beta = 0.5)
{
  c2 = c^2
  num_points = nrow(x)
  dim_x = ncol(x)
  imq_ksd_sum = 0

  #Calculate KSD
  for(i in 1:num_points){
    for(j in i:num_points){
      x1 = x[i,]
      x2 = x[j,]
      gradlogp1 = gradlogp[i,]
      gradlogp2 = gradlogp[j,]

      diff = x1-x2
      diff2 = sum(diff^2)

      base = diff2 + c2
      base_beta = base^(-beta)
      base_beta1 = base_beta/base

      kterm_sum = sum(gradlogp1*gradlogp2)*base_beta
      coeffgrad = -2.0 * beta * base_beta1
      gradx1term_sum = sum(gradlogp1*(-diff))*coeffgrad
      gradx2term_sum = sum(gradlogp2*diff)*coeffgrad
      gradx1x2term_sum = (-dim_x + 2*(beta+1)*diff2/base)*coeffgrad
      m <- 1 + 1*(i!=j)
      imq_ksd_sum = imq_ksd_sum + m*(kterm_sum + gradx1term_sum + gradx2term_sum + gradx1x2term_sum)
    }
  }
  imq_ksd = sqrt(imq_ksd_sum)/num_points
  return(imq_ksd)
}