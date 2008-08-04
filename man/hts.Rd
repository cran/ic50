\name{hts}
\alias{hts.96}
\alias{hts.384}
\title{High-throughput evaluation of compound screens}
\description{
  Simultaneous evaluation of a large number of compound screens on 96- and 384-well plates.
}
\usage{
hts.96(indir=".",plates=2,measure=NULL,control=NULL,dilution=NULL,inhib=NULL,
       normalize="mean",graphics="mean",outdir="./results")
hts.384(indir=".",plates=2,measure=NULL,control=NULL,dilution=NULL,inhib=NULL,
        normalize="single",graphics="mean",outdir="./results")
}
\arguments{
  \item{indir}{A character specifying the directory which contains the raw
    data files.}
  \item{plates}{Number of plates used for each experiment.}
  \item{measure}{Configuration file for the locations of the measurement wells.}
  \item{control}{Configuration file for the locations of the control wells.}
  \item{dilution}{Configuration for the concentrations in each measurement.}
  \item{inhib}{Vector of real numbers between 0 and 1 specifying the percentage of
    inhibition to compute concentrations for. Defaults to 0.5 for all compounds.}
  \item{normalize}{Method to normalize the measurement by the
    controls. If \code{"mean"}, the mean of the controls specified by
    \code{control} is used; \code{"single"} requires one individual
    control well per measurement well.}
  \item{graphics}{A character specifying the plotting method. For
    \code{"mean"}, a killcurve of the mean values of the measurement
    series is given, whereas one curve is plotted for each if \code{"single"} is
    specified. For \code{"fitted"}, a sigmoid-shaped derivation of the
    logistic model is fitted to the data.}  
  \item{outdir}{The directory where the results will be written.}
}
\details{
  In cytotoxicity screens of chemical compounds, biological activity is
  typically quantified by the concentration for which a particular
  fraction (typically 0.5) of the cells has been killed after a predefined
  treatment period. For this purpose, all concentrations are plotted
  against the percentages of cells still being alive under this
  treatment, forming a kill-curve under which the preimage of the 0.5
  point is defined as the half-maximum inhibitory concentration
  (IC50). For high-throughput screens (HTS), in particular, the
  evaluation of the data needs to be performed in an automatic fashion.
  
  The \code{hts.96} and \code{hts.384} functions provide a powerful tool
  to simultaneously evaluate all data in the specified input directory
  \code{indir}. The data files are handled in groups of the size specified
  by \code{plates} and the file names should be arranged in a way that
  two plates with replicates for the same measurements
  are displayed one below the other in a file browser. The data are
  expected to be arranged in tab-delimited text files as the typical
  output of appropriate signal readers. Just as for the
  evaluation of a single measurement, the design must be specified by
  tab-delimited files for \code{measure}, \code{control} and
  \code{dilution}. Details on these are given in the manual of the
  \code{\link{default384_measure}} and \code{\link{default384_control}} files.

  For each compound in the screen and each group of data files, a
  graphics output is given in the file \code{"killcurves.pdf"} in the
  current workspace directory. In addition, the text file
  \code{"ic50.txt"} contains a tab-delimted table with the same
  evaluation as for the \code{\link{ic50.96}} and \code{\link{ic50.384}}
  functions but for all experiments one below the other.
}
\value{
  A matrix with the following columns:
  \item{firstfile}{Filename of the first input file.}
  \item{ic}{The inhibitory concentrations for the respective compounds.}
  \item{cl}{Lower 0.95 confidence limits for the IC values.}
  \item{cu}{Upper 0.95 confidence limits for the IC values.}
  \item{inhib}{The inhibition percentages for the concentrations \code{ic}.}
  \item{maxsd}{Maximum of the standard deviations at the measured
    concentrations as determined from the single replicates.}
  \item{cv}{Coefficient of variation of the IC values as determined from the single replicates.}
  \item{area}{Area between the measured curve and an ideal kill-curve.}
}
\references{
Zhang JH, Chung TD, Oldenburg KR (1999): A simple statistical parameter
for use in evaluation and validation of high-throughput screening
assays. Journal of Biomolecular Screening 4, pp. 67-73
 
Tallarida RJ (2000): Drug synergism and dose-effect data
analysis. Chapman & Hall, New York
%Michel K, Zander T, Frommolt P et al. (2008): An integrated functional
%genomics approach for prediction of drug activity in non-small cell lung
%cancer.
}
\author{
  Peter Frommolt, University of Cologne \email{peter.frommolt@medizin.uni-koeln.de}\cr
  \url{http://www.medizin.uni-koeln.de/kai/imsie/Peter.Frommolt/}  
}
\examples{
#Example from a non-small cell lung cancer (NSCLC) cell line screen. In
#total, 84 samples were screened. The evaluation is exemplarily shown for
#the cell lines A549, Calu1, H322 and HCC2429.

data(A549_1,A549_2,Calu1_1,Calu1_2,H322_1,H322_2,HCC2429_1,HCC2429_2)
dir.create("NSCLC_screen")
write.table(A549_1,file="NSCLC_screen/A549_1.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(A549_2,file="NSCLC_screen/A549_2.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(Calu1_1,file="NSCLC_screen/Calu1_1.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(Calu1_2,file="NSCLC_screen/Calu1_2.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(H322_1,file="NSCLC_screen/H322_1.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(H322_2,file="NSCLC_screen/H322_2.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(HCC2429_1,file="NSCLC_screen/HCC2429_1.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(HCC2429_2,file="NSCLC_screen/HCC2429_2.txt",row.names=FALSE,col.names=FALSE,sep="\t")

data(mpi384_measure,mpi384_control,mpi384_dilution)
write.table(mpi384_measure,file="mpi384_measure.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(mpi384_control,file="mpi384_control.txt",row.names=FALSE,col.names=FALSE,sep="\t")
write.table(mpi384_dilution,file="mpi384_dilution.txt",row.names=FALSE,col.names=FALSE,sep="\t")

print(hts.384(indir="NSCLC_screen",
              measure="mpi384_measure.txt",control="mpi384_control.txt",dilution="mpi384_dilution.txt",
              inhib=rep(0.5,7),outdir="NSCLC_results",normalize="mean"))
}
\keyword{design}
\keyword{hplot}
\keyword{dplot}
\keyword{nonlinear}
\keyword{smooth}
