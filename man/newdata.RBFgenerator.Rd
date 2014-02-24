\name{newdata.RBFgenerator}
\alias{newdata.RBFgenerator}
\alias{newdata}

\title{ Generate semi-artificial data with RBF generator}
\description{
   Using a generator build with \code{\link{rbfDataGen}} the method generates \code{size} new instances.
    }
\usage{
\method{newdata}{RBFgenerator}(object, size, var=c("estimated","Silverman"), ... )
}

\arguments{
  \item{object}{ An object of type RBFgenerator containing a generator structure as returned by \code{\link{rbfDataGen}}. }
  \item{size}{ A number of instances to generate.}
  \item{var}{The method of kernel width (variance) estimation. Supported options are \code{"estimated"} and \code{"Silverman"}.}
  \item{... }{ Additional parameters passed to future extended versions of the function. } 
}
\details{
   The function uses the \code{object} structure as returned by \code{\link{rbfDataGen}} containing descriptions of 
   the Gaussian kernels which model the original data. The kernels are used to generate  a required number of new instances.
   
   The kernel width of provided kernels can be set in two ways. By setting \code{var="estimated"} the estimated spread of the
    training instances that have the maximal activation value for the particular kernel is used. 
    Using \code{var="Silverman"} width is set by the generalization of Silverman's rule of thumb to multivariate
     case (unreliable for larger dimensions). 
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
irisGenerator<- rbfDataGen(Species~.,iris)

# use the generator to create new data
irisNew <- newdata(irisGenerator, size=200)

#inspect properties of the new data
plot(irisNew, col = irisNew$Species) #plot generated data
summary(irisNew)
}


\keyword{ datagen }
\keyword{ multivariate }