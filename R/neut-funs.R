
#' Calculate neutralization given predicted titer 50
#' @param pt50 positive value, ratio of concentration of IC80
#' @return scalar or numerical vector
convert_pt50_neut = function(pt50, hill = 1) {1 - 1/(1+ pt50^hill)}

#' Calculate neutralization given predicted titer 80
#' @param pt80 positive value, ratio of concentration of IC80
#' @return scalar or numerical vector
convert_pt80_neut = function(pt80, hill = 1) {1 - 1/(1+ (pt80 * 4 ^(1/hill))^hill)}


#' Calculate predicted titer 80 (pt80) given neutralization
#' @param neut neutralization between 0 and 1
#' @return scalar or numerical vector
convert_neut_pt80 = function(neut, hill = 1) {4^(-1/hill) * neut/(1 - neut)}

#' Calculate predicted titer 50 (pt50) given neutralization
#' @param neut neutralization between 0 and 1
#' @return scalar or numerical vector
convert_neut_pt50 = function(neut, hill = 1) {neut/(1 - neut)}

#' Calculate IIP given neutralization
#' @param neut neutralization between 0 and 1
#' @return scalar or numerical vector
convert_neut_iip = function(neut) pmin(-log10(1 - neut), -log10(.Machine$double.eps))

#' Calculate IIP given PT80 (for input into AMP PE curve)
#' @param pt80 positive value, ratio of concentration of IC80
#' @return scalar or numerical vector
convert_pt80_iip = function(pt, hill = 1) {
  pmin(-log10(1 - convert_pt80_neut(pt80, hill = hill)), -log10(.Machine$double.eps))


# these arent ready
  #' #' Calculate general predicted titer given neutralization
  #' #' @param neut neutralization between 0 and 1
  #' #' @param titer_level between 1 and 100
  #' #' @return scalar or numerical vector
  #' convert_neut_pt = function(neut, titer_level = 50, hill = 1) {
  #'   4^(-1/hill*xx) * neut/(1 - neut)
  #' }
  #'

  #' #' Calculate neutralization given generalized titer
  #' #' @param pt positive value, ratio of concentration of given IC
  #' #' @param titer_input level for PT input (ex., 50 -> PT50)
  #' #' @param titer_target level for PT output (ex., 80 -> PT80)
  #' #' @param .titer_as_pct default T, assumed input is percent, converts values < 1
  #' #' @return scalar or numerical vector
  #' convert_titer_base = function(pt,
  #'                               titer_input = 50, titer_target = 80,
  #'                               hill = 1, .titer_as_pct = T) {
  #'   stopifnot(titer_input > 0 & titer_input < 100)
  #'   if(.titer_as_pct & titer_input < 1) titer_input = titer_input * 100
  #'   stopifnot(titer_target > 0 & titer_target < 100)
  #'   if(.titer_as_pct & titer_target < 1) titer_target = titer_target * 100
  #'
  #' }
  #'
  #' convert_pt50_pt80 = function(pt50, hill = 1) convert_titer_base(pt50, 50, 80, hill = hill)
  #' convert_pt80_pt50 = function(pt80, hill = 1) convert_titer_base(pt80, 80, 50, hill = hill)

  #' #' Calculate neutralization given generalized titer
  #' #' @param pt positive value, ratio of concentration of given IC
  #' #' @param titer_level between 1 and 100
  #' #' @param .titer_as_pct default T, assumed input is percent, converts values < 1
  #' #' @return scalar or numerical vector
  #' convert_pt_neut = function(pt, titer_level = 50, hill = 1, .titer_as_pct = T) {
  #'   stopifnot(titer_level > 0 & titer_level < 100)
  #'   if(.titer_as_pct & titer_level < 1) titer_level = titer_level * 100
  #'
  #'   1 - 1 / (1 + (pt * xxx ^ (1 / hill)) ^ hill)
  #' }
}
