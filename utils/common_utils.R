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