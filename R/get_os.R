#' Internal function
#'
#' Detect the user's operating system (OS).
#'
#' See this link for further details:
#' https://www.r-bloggers.com/identifying-the-os-from-r/
#'
#' @examples \dontrun{
#' get_os()
#' }
#' @keywords internal
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
