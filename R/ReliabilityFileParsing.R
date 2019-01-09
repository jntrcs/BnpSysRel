

#A set of numbers that is marked as "needsBuilt"--everything that appears after a colon
#Any numbers not after colon are "ready". Once a merge is done move the number after the colon from needsBuilt to ready


#' Estimates system reliability for a specified system.
#'
#' @param file A string with the file path containing the reliability diagram (see details)
#' @param priorList A list with named elements containing bsp object representing the prior for each component in the reliability diagram
#' @param dataList A list with named elements for each component or subsystem for which reliability was measured
#'
#' @return A named list of posterior bsp objects, one for each named component in the system.
#' @export
#'
#' @section Details:
#' The file is a text file specifying your reliability block diagram. Parallel and series subsystems can be specified as in the following example:
#'
#' S(Engine, GasDelivery):GasPropulsion
#'
#' S(Motor, Batteries, Controller, Belt):ElectricPropulsion
#'
#' P(GasPropulsion, ElectricPropulsion):Propulsion
#'
#' S and P represent series and parallel relationships respectively.
#' The first line of this diagram can be read as
#' "Engine and GasDelivery are components in series
#' of a subsystem named GasPropulsion." For each component on the right hand
#' side of an equation, the prior of that component will be determined
#' by the posterior of its subcomponents. Data, if available, can be used to update its
#' posterior by adding the data to dataList (e.g. assigning a matrix to dataList$GasPropulsion).
#' For each component on the left hand side of the colon, priors can be defined
#' by adding a bsp object to priorList. If no prior is provided, a non-informative
#' prior will be used, with a warning.
#'
#' The following syntax conventions for the file will be enforced:
#'\itemize{
#' \item Each line begins with S( or P( indicating series or parallel, one relationship per line
#'
#'\item Components can be uppercase, lowercase, and include numbers, but no special characters.
#'
#'\item All relationships must be named by adding a colon and a label after the end of the relationship
#'
#' \item Two or more elements can be specified in a relationship
#'
#' \item A check will be done to ensure there are no circular dependencies among the named components.
#' For example, Propulsion cannot be used as an input to GasPropulsion because GasPropulsion
#' is an input to Propulsion
#'}
#' @examples
#'
estimateSystemReliability<-function(file, priorList, dataList){

}


#' Creates a non-informative prior
#'
#' @param data The data matrix, which is used to ensure only the points necessary appear on the support
#'
#' @return A bsp object with 0 precision.
#' @export
#'
createUninformativePrior<-function(data){
  #mindata<-min(data[,1])
  maxdata<-max(data[,1])
  bsp(c(maxdata), centeringMeasure = c(.999), precision = 0)

}
