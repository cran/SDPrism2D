\name{sdprism2d}
\alias{sdprism2d}
\title{Visualizing the Standard Deviation as the Size of a Prism}
\usage{
sdprism2d(data, hlim = NULL, xyscale = NULL)
}
\arguments{
  \item{data}{The data that a user inputs, usually a vector of values.}
  \item{hlim}{Optional, 4 by default. The height limit for the plot of step 2, step3, and step 4.}
  \item{xyscale}{Optional, 4 by default. The ratio of scales between the x-axis and the y-axis.}
}
\value{
No return value, the function will open a new window and display the graphs of the 4 steps of visualizing the standard deviation.
}
\description{
We visualize the standard deviation of a data set as the size of a prism whose volume equals the total volume of several prisms made from the Empirical Cumulative Distribution Function.
}
\examples{
  sdprism2d(c(10,18,23,30,36),4,4)
}
