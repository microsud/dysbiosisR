% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combinedShannonJSD.R
\name{combinedShannonJSD}
\alias{combinedShannonJSD}
\title{Combined Alpha Beta Diveristy Based Score}
\arguments{
\item{x}{A phyloseq object}

\item{reference_samples}{Vector of samples to use as reference.}
}
\value{
A data frame with combinedShannonJSD score and sample information
}
\description{
Combined Alpha Beta Diveristy Based Score
}
\details{
Calculates a dysbiosis score that combined Shannon
         diversity and Jensen–Shannon divergence.

The combined alpha-beta diversity approach was used by
\code{Santiago M E et al. 2019}. This approach uses Shannon diversity as the
alpha diversity measure and Jensen–Shannon divergence as the beta diversity
measure. The score is mean difference of Shannon diversity between test sample
and all references samples multiplied by the mean JSD of the test sample to all
reference samples. When calculating this score for reference samples, the
sample being used is excluded from calculating means for alpha and beta diversity.
}
\examples{
data("WirbelJ_2018")
library(phyloseq)
ps <- WirbelJ_2018
# Define controls as reference samples
ref.samples <- sample_names(subset_samples(WirbelJ_2018,
                                           disease == "healthy"))
alpha.beta.dysbiosis <- combinedShannonJSD(ps,
                                           reference_samples = ref.samples)
head(alpha.beta.dysbiosis)

}
\references{
\itemize{
\item{}{Santiago M et al. (2019).
Microbiome predictors of dysbiosis and VRE decolonization in patients with
recurrent C. difficile infections in a multi-center retrospective study.
\emph{AIMS Microbiol}, 5:1–18.}
}
}
\author{
Sudarshan A. Shetty
}
