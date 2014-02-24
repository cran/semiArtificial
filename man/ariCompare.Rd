\name{ariCompare}
\alias{ariCompare}


\title{ Evaluate clustering similarity of two data sets with ARI}
\description{
   Similarity of two data sets is compared with Adjusted Rand Index (ARI).
    }
\usage{
ariCompare(data1, data2)
}

\arguments{
  \item{data1}{ A \code{data.frame} containing the reference data.}
  \item{data2}{ A \code{data.frame} with the same number and names of columns as \code{data1}.}
}
\details{
   The function compares data stored in \code{data1} with \code{data2} by first performing partitioning
    around medoids (PAM) clustering on \code{data1}. Instances from \code{data2} are than assigned to the
     cluster with the closest medoid. 
   In second step PAM clustering is performed on \code{data2} and instances from \code{data1} are assigned to
    the clusters with
   closest medoids. The procedure gives us two clusterings on the same instances which we can compare
    with ARI. The higher the value of ARI the more similar the two data sets.
   }
\value{
  The method returns a value of ARI.
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