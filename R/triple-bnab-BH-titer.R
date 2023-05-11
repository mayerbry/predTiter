# BH titer: common factor applied to each individual Ab concentration that reduces
# neutralization to 50%


#' Theoretical calculation for experimental (pooled sera) neutralization titer
#' under BH interaction for 3 antibodies
#'
#' @param pt50_1 ratio of concentration to IC50 for first Ab
#' @param pt50_2 ratio of concentration to IC50 for second Ab
#' @param pt50_3 ratio of concentration to IC50 for third Ab
#' @param titer_target goal combination titer endpoint: 0.5 or 50 = PT50, 0.8 or 80 = PT80. Pct equivalence fails for tiny titer levels (<PT1): if you want PT<1, must input < 0.01
#' @param min_titer_return smallest triple BH titer returned (default 0.01)
#' @param NaN_as_NA arbitrarily large titers returned as missing
#'
#' This function is vectorized for use over PK
#'
#' @return
#' @export
#'
#' @examples
calc_3bnab_BHtiter = function(pt50_1, pt50_2, pt50_3, titer_target = 50,
                              min_titer_return = 1e-2,NaN_as_NA = F){
  stopifnot(titer_target > 0 & titer_target < 100)

  # under .titer_as_pct, titer_target < 1 taken as is, > 1 converted
  if(titer_target > 1) titer_target = titer_target / 100

  # if a row is all zeros, this is tracked then 0 is returned
  all_zeros = which(pt50_1 < sqrt(.Machine$double.eps) &
                      pt50_2 < sqrt(.Machine$double.eps) &
                      pt50_3 < sqrt(.Machine$double.eps))
  # any zeros has a closed-form solution
  any_zeros = which(pt50_1 < sqrt(.Machine$double.eps) |
                      pt50_2 < sqrt(.Machine$double.eps) |
                      pt50_3 < sqrt(.Machine$double.eps))
  # managing numerical precision
  pt50_1 = pmax(pt50_1, sqrt(.Machine$double.eps))
  pt50_2 = pmax(pt50_2, sqrt(.Machine$double.eps))
  pt50_3 = pmax(pt50_3, sqrt(.Machine$double.eps))

  tmp_pt50_dat = data.frame(a = pt50_1, b = pt50_2, c = pt50_3)

  # the two sufficient values for 2-bnab calculations (as needed)
  tmp_pt50_dat$median_row = apply(tmp_pt50_dat, 1, median)
  tmp_pt50_dat$max_row = pmax(tmp_pt50_dat$a, tmp_pt50_dat$b, tmp_pt50_dat$c)

  tmp_pt50_dat$PT = NA_real_

  # if there are any zeros, the information only exists in at most two of the input titers
  if(length(any_zeros) > 0){

    tmp_pt50_dat$PT[any_zeros] = .calc_2bnab_BHPTxx(
      pt50_1 = tmp_pt50_dat$median_row[any_zeros],
      pt50_2 = tmp_pt50_dat$max_row[any_zeros],
      titer_target = titer_target,
      NaN_as_NA = NaN_as_NA)

    if(length(tmp_pt50_dat[-any_zeros,] > 0)){
      tmp_pt50_dat$PT[-any_zeros] = .calc_3bnab_BHPTxx(
        tmp_pt50_dat[-any_zeros,],
        titer_target = titer_target,
        min_titer_return = min_titer_return,
        NaN_as_NA = NaN_as_NA)
    }

  } else{
    tmp_pt50_dat$PT = .calc_3bnab_BHPTxx(
      tmp_pt50_dat,
      titer_target = titer_target,
      min_titer_return = min_titer_return,
      NaN_as_NA = NaN_as_NA)
  }

  tmp_pt50_dat$PT[all_zeros] = 0
  tmp_pt50_dat$PT
}

#' Theoretical calculation for experimental (pooled sera) neutralization titer
#' under BH interaction for 2 antibodies (closed form solution)
#'
#' @param PT50_1 ratio of concentration to IC50 for first Ab
#' @param PT50_2 ratio of concentration to IC50 for second Ab
#' @param titer_target goal combination titer endpoint: 0.5 or 50 = PT50, 0.8 or 80 = PT80. Pct equivalence fails for tiny titer levels (<PT1): if you want PT<1, must input < 0.01
#' @param NaN_as_NA arbitrarily large titers returned as missing
#'
#' This function is vectorized for use over PK
#'
#' @return
#' @export
#'
#' @examples
calc_2bnab_BHtiter = function(pt50_1, pt50_2, titer_target = 0.5, NaN_as_NA = F){

  # managing numerical precision
  all_zeros = which(pt50_1 < sqrt(.Machine$double.eps) & pt50_2 < sqrt(.Machine$double.eps))

  pt50_1 = pmax(pt50_1, sqrt(.Machine$double.eps))
  pt50_2 = pmax(pt50_2, sqrt(.Machine$double.eps))
  out = .calc_2bnab_BHPTxx(pt50_1, pt50_2, titer_target, NaN_as_NA)
  if(length(all_zeros > 0)) out[all_zeros] = 0
  out
}

.calc_2bnab_BHPTxx = function(pt50_1, pt50_2, titer_target, NaN_as_NA) {
  conj = (pt50_1 + pt50_2) ^ 2 + 4 * pt50_1 * pt50_2 * (titer_target / (1 - titer_target))
  x = (-(pt50_1 + pt50_2) + sqrt(conj)) / (2 * pt50_1 * pt50_2)
  if (any(is.nan(x))) {
    if (!NaN_as_NA)
      stop(paste("BH PT50 nAn, input:", pt50_1, pt50_2))
    x <- NA_real_
  }
  out = 1 / x

  # this handles some cluster issues with high values
  infinites = which(is.infinite(out))
  if (length(infinites > 0))
    out[infinites] = pmax(pt50_1, pt50_2)[infinites]

  out

}

.BH_cubic_fun = function(x, pt50_1, pt50_2, pt50_3, titer_target = 0.5){
  termA = pt50_1 * pt50_2 * pt50_3
  termB =  pt50_1 * pt50_2 +  pt50_1 * pt50_3 +  pt50_2 * pt50_3
  termC =  pt50_1 + pt50_2 + pt50_3
  termD = titer_target/(1 - titer_target)
  x^3 * termA + x^2 * termB + x * termC - termD
}

.calc_3bnab_BHPTxx = function(pt50dat, titer_target, min_titer_return, NaN_as_NA){

  # returns a fixed value when root outside of the interval ( purrr's tryCatch)
  safe_uniroot = purrr::safely(uniroot,
                               otherwise = list(root = 1/min_titer_return))

  x = purrr::map_dbl(1:nrow(pt50dat), function(i){

    res = safe_uniroot(.BH_cubic_fun, interval = c(0, 1/min_titer_return),
                  pt50_1 = pt50dat$a[i], pt50_2 = pt50dat$b[i], pt50_3 = pt50dat$c[i],
                  titer_target = titer_target, tol = .Machine$double.eps ^ 0.5)
    res$result$root

  })


  if(any(is.nan(x))) {
    if(!NaN_as_NA) stop(paste("BH PT50 nAn, input:", pt50_1, pt50_2))
    x <- NA_real_
  }
  out = 1/x

  # this handles some cluster issues with high values
  infinites = which(is.infinite(out))
  if(length(infinites > 0)) out[infinites] = pmax(pt50_1, pt50_2, pt50_3)[infinites]

  out

}
