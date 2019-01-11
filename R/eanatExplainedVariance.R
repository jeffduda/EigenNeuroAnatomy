eanatExplainedVariance <- function( reference, basis, by.component=FALSE, r.squared=FALSE ) {

  if ( dim(reference)[2] != dim(basis)[2] ) {
    stop("Incompatible size of inputs")
  }

  nvec = dim(basis)[1]
  outvec = rep(NA,nvec)
  y = t(basis)

  if (r.squared) {
    sst = sum( (reference - mean(reference))^2 )

    for ( j in 1:nvec ) {

      linmod = NA
      if (by.component==TRUE) {
        linmod = lm( reference ~  ( reference %*% y[,j] ))
      }
      else {
        submat = as.matrix(y[,1:j])
        linmod = lm( reference ~  ( reference %*% submat ) )
      }
      linmodResid = resid(linmod)
      ssr = sum(linmodResid*linmodResid)

      outvec[j] = 1-(ssr/sst)
    }

  }
  else {
    for ( j in 1:nvec ) {
      linmod=NA
      if (by.component==TRUE) {
        linmod = lm( reference ~  ( reference %*% y[,j] ))
      }
      else {
        submat = as.matrix(y[,1:j])
        linmod = lm( reference ~  ( reference %*% submat ) )
      }

      reconmat = predict( linmod )

      meanFit = 0
      for ( i in 1:ncol(reference) )
        {
        corValue = cor(reference[,i],reconmat[,i])
        meanFit = meanFit + corValue*corValue
        }

      outvec[j] = meanFit / ncol( reference )

    }
  }

  return(outvec)

}
