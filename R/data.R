#' Gene expression profiling of chronic lung disease for the Lung Genomics Research Consortium (GSE47460)
#' Microarray data from total RNA extracted from whole lung homogenate from subjects undergoing thoracic surgery.
#' There are two platforms in this GEO dataset and this data contains the platform GPL14550. We substract COPD and ILD patients in this dataset.
#' @format a list of four components.
#' \itemize{
#' \item{outcome: }{outcome (FEV1) of 331 samples}
#' \item{Expression: }{Expression matrix with 331 samples(Cols) and 12958 genes(Rows)}
#' \item{Covariates: }{age and gender for 331 samples}
#' \item{label: }{Diagosis for 331 samples, either COPD or ILD}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460}
"GSE47460_GPL14550"

#' Gene expression profiling of chronic lung disease for the Lung Genomics Research Consortium (GSE47460)
#' Microarray data from total RNA extracted from whole lung homogenate from subjects undergoing thoracic surgery.
#' There are two platforms in this GEO dataset and this data contains the platform GPL6480. We substract COPD and ILD patients in this dataset.
#' @format a list of four components.
#' \itemize{
#' \item{outcome: }{outcome (FEV1) of 136 samples}
#' \item{Expression: }{Expression matrix with 136 samples(Cols) and 12958 genes(Rows)}
#' \item{Covariates: }{age and gender for 136 samples}
#' \item{label: }{Diagosis for 136 samples, either COPD or ILD}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460}
"GSE47460_GPL6480"

#' Gene expression profiling of METARBIC breast cancer dataset, 275 triple negative breast cancer samples are picked.
#' @format a list of five components.
#' \itemize{
#' \item{Expression: }{Expression matrix with 275 samples(Rows) and 13964  genes(Cols)}
#' \item{OS: }{overall survival}
#' \item{OS.event: }{censoring indicators where 1 means events and 0 means censoring}
#' \item{covariate: }{covariates matrix: age, tumor size, number of lymph nodes and tumor grade}
#' \item{PAM50: }{PAM50 labels}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460}
"Data.Metabric"

#' Gene expression profiling of SCAN-B breast cancer dataset, 151 triple negative breast cancer samples are picked.
#' @format a list of five components.
#' \itemize{
#' \item{Expression: }{Expression matrix with 151 samples(Rows) and 18964  genes(Cols)}
#' \item{OS: }{overall survival}
#' \item{OS.event: }{censoring indicators where 1 means events and 0 means censoring}
#' \item{covariate: }{covariates matrix: age, tumor size, number of lymph nodes and tumor grade}
#' \item{PAM50: }{PAM50 labels}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60789}
"Data.ScanB"
