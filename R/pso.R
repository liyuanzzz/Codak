#' Psoriasis vulgaris data
#'
#' Data from a psoriasis vulgaris disease study to identify DE genes between two groups of samples.
#'     The results can be summarized as 18151 observations of signed p-values.
#'     This dataset also contains two other useful covariates for this problem.
#'
#' @docType data
#'
#' @format A data frame with the following three variables:
#' \describe{
#'     \item{signp}{the signed p-values. Corresponding p-values are obtained using the limma-voom method.}
#'     \item{len_gene}{the lengths of the gene coding regions}
#'     \item{tval_mic}{the test statistics from a related microarray study in Gudjonsson et al. (2010)}
#' }
#'
#' @keywords datasets
#'
#' @references Jabbari et al. (2014) Journal of Investigative Dermatology 132: 246-249
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/21850022}{PubMed})
#'
#' Gudjonsson et al. (2010) Journal of Investigative Dermatology 130: 1829-1840
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/20220767}{PubMed})
#'
#' @examples
#' data(pso)
#' x=data.frame(x1=pso$len_gene,x2=ifelse(is.na(pso$tval_mic),0,pso$tval_mic))
#' FDR=codak_group(pso$signp,x,c("s(x1,x2)","s(x1)"),group=1+as.numeric(is.na(pso$tval_mic)))
#' rejResult=(FDR<=0.1)
"pso"
