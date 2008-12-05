\name{design}
\alias{mpi96_measure}
\alias{mpi96_control}
\alias{mpi96_dilution}
\alias{mpi384_measure}
\alias{mpi384_control}
\alias{mpi384_dilution}
\alias{default96_measure}
\alias{default96_control}
\alias{default96_dilution}
\alias{default384_measure}
\alias{default384_control}
\alias{default384_dilution}
\docType{data}
\title{Configuration files specifying the design of a compound screen on
  84 NSCLC cell lines}
\description{
  The experimental setup in a compound screen on 84 non-small cell lung
  cancer (NSCLC) cell lines was designed according to the arrangement
  specified in these data sets. However, the functions for screen
  evaluation actually expect this information as files on the local
  harddisk such as \code{"mpi384_measure.txt"}, \code{"mpi384_control.txt"} and
  \code{"mpi384_dilution.txt"} which are installed together with the
  package. Default examples of these configuration files are included as
  \code{"default96_measure.txt"}, \code{"default384_measure.txt"} etc.

  In the \code{"measure"} file, there must be one row specified for each
  measurement series. This consists of the respective compound name
  followed by tab-delimited information on the wells where the
  measurements are located. Each of these must be given as a comma-delimited
  pair of coordinates. The same is expected for the control wells
  in the \code{"control"} file. Here, a particular well can be specified
  multiple times and will then be repeatedly used as a control well for
  the signal intensity without a compound being applied. If the
  \code{normalize} argument to the functions is specified as \code{"single"}, there must be
  one control well for each single measurement; if on the other hand,
  \code{"mean"} is selected, an arbitrary number of wells can be
  specified and the mean of those values is used for
  normalization. Finally, the number of rows and the row names in the
  \code{"dilution"} file must equal those in the preceding two
  configuration files. Each row contains the compound name followed by a
  tab-delimited list of the concentrations used for the respective
  measurement series. It should by obvious that the number of concentrations in one row must equal the
  number of wells in the \code{"measure"} file for each row. However,
  the number of control wells can be distinct from these if the
  \code{normalize} argument is set to \code{"mean"} such that the mean
  of the respective control row is taken.

  Importantly, the number of rows must be equal in all three files as
  well as the row names, where case-sensitivity and literal equality has
  to be carefully verified. The easiest way to create the configuration
  files is to simply start the GUI using \code{ic50()} which
  automatically creates a default version to be modified by the
  user. After having saved this configuration, it can be repeatedly used
  for screen evaluations as long as the experimental setup is not
  changed.

  A step-by-step tutorial document describing how to prepare the data and
  configuration  is included in the \code{ic50} package.
}
\references{
  Frommolt P, Thomas RK (2008): Standardized high-throughput evaluation
  of cell-based compound screens. BMC Bioinformatics, 9(1): 475
  
  Sos ML, Michel K, Zander T, Frommolt P, Weiss J, et al. (2008): Predicting drug 
  activity in non-small cell lung cancer based on genetic lesions. Submitted
}
\keyword{design}
\keyword{datasets}
