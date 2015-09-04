\name{newdata}
\alias{newdata.RBFgenerator}
\alias{newdata}
\alias{newdata.TreeEnsemble}

\title{ Generate semi-artificial data using a generator}
\description{
   Using a generator build with \code{\link{rbfDataGen}} or \code{\link{treeEnsemble}} the method generates \code{size} new instances.
    }
\usage{
\method{newdata}{RBFgenerator}(object, size, var=c("estimated","Silverman"), 
                               classProb=NULL, defaultSpread=0.05, ... )
\method{newdata}{TreeEnsemble}(object, size=1, ...) 

}

\arguments{
  \item{object}{ An object of class RBFgenerator or TreeEnsemble containing a generator structure as returned by \code{\link{rbfDataGen}} or \code{\link{treeEnsemble}}, respectively.  }
  \item{size}{ A number of instances to generate.}
  \item{var}{For RBFgenerator type parameter \code{var} determines the method of kernel width (variance) estimation. Supported options are \code{"estimated"} and \code{"Silverman"}.}
  \item{classProb}{A vector of desired class value probability distribution. Default value \code{classProb=NULL} uses probability distribution of the generator training instances. }
  \item{defaultSpread}{For RBFgenerator type parameter is a numeric value replacing zero spread in case \code{var="estimated"} is used. The value \code{defaultSpread=NULL} keeps zero spread values.  }
  
  \item{... }{ Additional parameters passed to density estimation functions kde, logspline, and quantile.} 
}
\details{
   The function uses the \code{object} structure as returned by \code{\link{rbfDataGen}} or \code{\link{treeEnsemble}}. 
   In case of \code{RBFgenerator} the object contains descriptions of the Gaussian kernels, which model the original data. 
   The kernels are used to generate  a required number of new instances.
    The kernel width of provided kernels can be set in two ways. By setting \code{var="estimated"} the estimated spread of the
    training instances that have the maximal activation value for the particular kernel is used. 
    Using \code{var="Silverman"} width is set by the generalization of Silverman's rule of thumb to multivariate
     case (unreliable for larger dimensions). 
     
  In case of TreeEnsemble generator no additional parameters are needed, except for the number of generated instances.
 }
\value{
  The method returns a \code{data.frame} object  with required number of instances.
  }

\author{ Marko Robnik-Sikonja}

\seealso{
\code{\link{rbfDataGen}}.
}
\examples{
# inspect properties of the iris data set
plot(iris, col=iris$Species)
summary(iris)

# create RBF generator
irisRBF<- rbfDataGen(Species~.,iris)
# create treesemble  generator
irisEnsemble<- treeEnsemble(Species~.,iris,noTrees=10)


# use the generator to create new data with both generators
irisNewRBF <- newdata(irisRBF, size=150)
irisNewEns <- newdata(irisEnsemble, size=150)

#inspect properties of the new data
plot(irisNewRBF, col = irisNewRBF$Species) #plot generated data
summary(irisNewRBF)
plot(irisNewEns, col = irisNewEns$Species) #plot generated data
summary(irisNewEns)
}


\keyword{ datagen }
\keyword{ multivariate }