\name{ariCompare}
\alias{ariCompare}


\title{ Evaluate clustering similarity of two data sets}
\description{
   Similarity of two data sets is compared with Adjusted Rand Index (ARI) and Fowlkes-Mallows index(FM).
    }
\usage{
ariCompare(data1, data2, ARI=TRUE, FM=TRUE) 
}

\arguments{
  \item{data1}{ A \code{data.frame} containing the reference data.}
  \item{data2}{ A \code{data.frame} with the same number and names of columns as \code{data1}.}
  \item{ARI}{ A logical value indicating whether ARI shall be reported.}
  \item{FM}{ A logical value indicating whether FM index shall be reported.}
}
\details{
   The function compares data stored in \code{data1} with \code{data2} by first performing partitioning around medoids (PAM) clustering on \code{data1}. 
   Instances from \code{data2} are than assigned to the cluster with the closest medoid. 
   In second step PAM clustering is performed on \code{data2} and instances from \code{data1} are assigned to the clusters with closest medoids. 
   The procedure gives us two clusterings on the same instances which we can compare using either ARI or FM. 
   The higher the value of ARI/FM the more similar are the two data sets.
   
   For random clustering ARI returns a value around zero (negative values are possible) and for perfectly matching clustering ARI is 1. FM values are strictly in [0, 1].
   }
\value{
  The method returns a value of a list containing ARI and/or FM, depending on the parameters.
  }

\author{ Marko Robnik-Sikonja}

\seealso{
\code{\link{newdata.RBFgenerator}}.
}
\examples{
# use iris data set

# create RBF generator
irisGenerator<- rbfDataGen(Species~.,iris)

# use the generator to create new data
irisNew <- newdata(irisGenerator, size=200)

# compare ARI computed on clustering with original and new data
ariCompare(iris, irisNew)

}


\keyword{ datagen }
\keyword{ multivariate }