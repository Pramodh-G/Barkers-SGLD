samp_z <- function(n=1, h = 0.2, dist = "normal")
{
   ret <- rnorm(n)
   if(dist == "normal")
   {
       ret <- rnorm(n, sd = h)
   }
   else if (dist == "bim_normal") {
      if(h > 1)
      {
          print("choose h less than 1. mean is not defined.")
          quit()
      }

      ret <- rnorm(n, mean = sqrt(1 - h^2), sd = h)
      # invert it to generate from bimodal distribution, possibly for each target.
      invert_or_not <- (runif(n) < 0.5)
      weight <- 2 * invert_or_not - 1
      ret <- ret * weight
   }
   return(ret)
}

log_proposal_dist <- function(x, mean = 0, h = 1, dist = "normal")
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
    return(ret)
}