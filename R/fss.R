### Fractions Skill Score
###   OBSOLETE "pure-R" implementation
### Fast algorithm based on
### N. Faggian et al. (2015)

### fss calculates for a list of scales and thresholds
### 1. get the 'summed_area_table' for a given threshold
### 2. use this for fast calculation of FSS at different scales

### TO DO: 1. check dimensioning, C-code compatibility ...
##         2. sat.fc is a vector (from C code), not a matrix!
##            so either reformat or change indexing -> what is faster?
##         3. return the 3 components in stead of fss? Then you can accumulate over time.


#' @param fcfield A data matrix
#' @param obfield A data matrix that must have exactly the same dimensions
#' @param windows Window sizes in x and y direction as two column
#'   All values must be odd!
#' @param thresholds A list of numerical thresholds
fss_old <- function(obfield,
  fcfield,
  windows = rbind(c(3, 3), c(5, 5)),
  thresholds = c(1, 5, 15),
  ...) {
  if (any(dim(obfield) != dim(fcfield))) {
    stop("FC and OBS must have same dimensions.")
  }
  nx <- dim(fcfield)[1]
  ny <- dim(fcfield)[2]
  nscales <- nrow(windows)
  result <- NULL

  for (thr in thresholds) {
    #    sat.fc <- .C("summed_area_table",
    #                indat=as.integer(forecast>=thr),
    #                nx=as.integer(nx),
    #                ny=as.integer(ny))$indat
    #    sat.fc <- matrix(sat.fc, nrow=nx)
    sat.fc <- apply(apply(fcfield >= thr, 1, cumsum), 1, cumsum)
    #    sat.obs <- .C("summed_area_table",
    #                indat=as.integer(obs>=thr),
    #                nx=as.integer(nx),
    #                ny=as.integer(ny))$indat
    #    sat.obs <- matrix(sat.obs, nrow=nx)
    sat.obs <- apply(apply(obfield >= thr, 1, cumsum), 1, cumsum)

    res <- data.frame(
      threshold = rep(thr, nscales),
      nx = rep(nx, nscales),
      ny = rep(ny, nscales),
      fss = rep(NA, nscales),
      fss1 = rep(NA, nscales),
      fss2 = rep(NA, nscales)
    )

    for (i in 1:nrow(windows)) {
      # works best with odd numbers (symmetric around center)
      sx2 <- floor(windows[i, 1] / 2)
      sy2 <- floor(windows[i, 2] / 2)

      x0 <- 1:nx - sx2
      x0[x0<=0] <- NA
      X0 <- rep(x0, ny)
      #      X0 <- rep(pmax(1:nx - sx2, 1), ny)
      X1 <- rep(pmin(1:nx + sx2, nx), ny)
      Y0 <- rep(pmax(1:ny - sy2, 1),  each = nx)
      Y1 <- rep(pmin(1:ny + sy2, ny), each = nx)

      frac.fc <- (sat.fc[rbind(X0, Y0)] + sat.fc[rbind(X1, Y1)]
        - sat.fc[rbind(X0, Y1)] - sat.fc[rbind(X0, Y1)]) / (sx*sy)
      frac.obs <- (sat.obs[rbind(X0, Y0)] + sat.obs[rbind(X1, Y1)]
        - sat.obs[rbind(X0, Y1)] - sat.obs[rbind(X0, Y1)]) / (sx*sy)
      res$fss1[i] <- sum( (frac.fc - frac.obs)^2 ) # /(nx*ny)
      res$fss2[i] <- sum(frac.fc^2  + frac.obs^2) #/(nx*ny)

    }
    res$fss <- 1 - res$fss1 / res$fss2
    if (is.null(result)) result <- res
    else result <- rbind(result, res)
  }
  result
}

#' Neighbourhood verification for ensembles
#'
#' @param .data A \code{"harp_ens_grid_df"} data frame or a \code{geolist}.
#' @param obs_col The name of the observations column in the
#'   \code{"harp_ens_grid_df"}, or if \code{.data} is a \code{geolist}, a
#'   \code{geofield}.
#' @param threshold A numeric vector of thresholds to for which to compute the
#'   neighbourhood scores. Can also be a character vector with numeric values
#'   between 0 and 1 preceded by a "q" to denote quantile thresholds.
#' @param radius The radius, in pixels of the neighbourhoods for which to
#'   compute the scores.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nbhd_verify <- function(
  .fcst, obs, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, ...
) {
  UseMethod("nbhd_verify")
}

#' @export
nbhd_verify.geofield <- function(
  .fcst, obs, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, ...
) {

  quantile_thresh <- FALSE
  thresholds <- harpCore::parse_thresholds(threshold)
  if (thresholds[["quantiles"]]) {
    quantile_thresh <- TRUE
    threshold <- quantile(obs, thresholds[["thresholds"]])
    thresholds_df <- tibble::tibble(
      quantile  = thresholds[["thresholds"]],
      threshold = threshold
    )
  }

  result <- cpp_neighborhood_scores(
    .fcst, obs, threshold, radius, comparator, include_low, include_high
  )

  if (quantile_thresh) {
    result <- dplyr::relocate(
      suppressMessages(dplyr::inner_join(result, thresholds_df)),
      dplyr::all_of("quantile"),
      .before = "threshold"
    )
  }

  result

}

#' @export
nbhd_verify.harp_det_grid_df <- function(
  .fcst, obs, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, ...
) {

  dplyr::mutate(
    dplyr::rowwise(.fcst),
    dplyr::across(
      {{obs}},
      ~list(
        nbhd_verify(
          .data[["fcst"]],
          .x,
          threshold,
          radius,
          comparator,
          include_low,
          include_high
        )
      ),
      .names = "nbhd_scores"
    )
  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::all_of(c("fcst")), -{{obs}}) %>%
    tidyr::unnest(dplyr::all_of("nbhd_scores"))

}

#' @export
nbhd_verify.harp_ens_grid_df <- function(
    .fcst, obs, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1, ...
) {
  result <- dplyr::bind_rows(
    ens_fss(
      .fcst, {{obs}}, threshold, radius,
      comparator, include_low, include_high
    ),
    ens_efss(
      .fcst, {{obs}}, threshold, radius,
      comparator, include_low, include_high, num_cores
    )
  )

  dplyr::bind_rows(
    result,
    dplyr::inner_join(
      ens_dfss(
        .fcst, harpCore::unique_col(result, "threshold"), radius,
        comparator, include_low, include_high, num_cores
      ),
      dplyr::distinct(
        result,
        dplyr::across(intersect(c("quantile", "threshold"), colnames(result)))
      )
    )
  )
}


#' Title
#'
#' @param x
#' @param threshold
#' @param radius
#' @param num_cores
#'
#' @return
#' @export
#'
#' @examples
ens_dfss <- function(
  x, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1
) {
  UseMethod("ens_dfss")
}

#' @export
ens_dfss.geolist <- function(
  x, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1
) {
  member_pairs <- unique_pairs(seq_along(x))
  if (num_cores > 1) {

    num_cores <- min(num_cores, parallel::detectCores(), length(member_pairs))
    result <- dplyr::bind_rows(
      parallel::mclapply(
        1:nrow(member_pairs),
        function(i) nbhd_verify(
          x[[member_pairs[i, 1]]], x[[member_pairs[i, 2]]],
          threshold, radius, comparator, include_low, include_high
        ),
        mc.cores = num_cores
      )
    )

  } else {

    result <- dplyr::bind_rows(
      lapply(
        1:nrow(member_pairs),
        function(i) nbhd_verify(
          x[[member_pairs[i, 1]]], x[[member_pairs[i, 2]]],
          threshold, radius, comparator, include_low, include_high
        )
      )
    )

  }

  summarise_fss(
    result, type = "dispersion", dx = attr(x[[1]], "domain")[["dx"]]
  )

}

#' @export
ens_dfss.harp_ens_grid_df <- function(
  x, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1
) {
  dplyr::mutate(
    dplyr::rowwise(x),
    dfss = list(
      ens_dfss(
        as_geolist(as.list(dplyr::pick(dplyr::matches("_mbr[[:digit:]]+")))),
        threshold,
        radius,
        comparator,
        include_low,
        include_high,
        num_cores
      )
    )
  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::where(is_geolist)) %>%
    tidyr::unnest(dplyr::all_of("dfss"))
}

#' @export
ens_efss <- function(
  x, y, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1
) {
  UseMethod("ens_efss")
}

#' @export
ens_efss.geolist <- function(
  x, y, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1
) {
  stopifnot(meteogrid::is.geofield(y))
  stopifnot(
    meteogrid::compare.geodomain(
      attr(x[[1]], "domain"),
      meteogrid::as.geodomain(y)
    )
  )

  quantile_thresh <- FALSE
  thresholds <- harpCore::parse_thresholds(threshold)
  if (thresholds[["quantiles"]]) {
    quantile_thresh <- TRUE
    threshold <- quantile(y, thresholds[["thresholds"]])
    thresholds_df <- tibble::tibble(
      quantile  = thresholds[["thresholds"]],
      threshold = threshold
    )
  }

  if (num_cores > 1) {

    num_cores <- min(num_cores, parallel::detectCores(), length(x))
    result <- dplyr::bind_rows(
      parallel::mclapply(
        x, nbhd_verify, y, threshold, radius,
        comparator, include_low, include_high,
        mc.cores = num_cores
      )
    )

  } else {

    result <- dplyr::bind_rows(
      lapply(
        x, nbhd_verify, y, threshold, radius,
        comparator, include_low, include_high
      )
    )

  }

  if (quantile_thresh) {
    result <- dplyr::relocate(
      suppressMessages(dplyr::inner_join(result, thresholds_df)),
      dplyr::all_of("quantile"),
      .before = "threshold"
    )
  }

  summarise_fss(
    result, type = "error", dx = attr(x[[1]], "domain")[["dx"]]
  )

}

#' @export
ens_efss.harp_ens_grid_df <- function(
  x, y, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, num_cores = 1
) {
  dplyr::mutate(
    dplyr::rowwise(x),
    dplyr::across(
      {{y}},
      ~list(
        ens_efss(
          as_geolist(as.list(dplyr::pick(dplyr::matches("_mbr[[:digit:]]+")))),
          .x,
          threshold,
          radius,
          comparator,
          include_low,
          include_high,
          num_cores
        )
      ),
      .names = "efss"
    )
  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::where(is_geolist)) %>%
    tidyr::unnest(dplyr::all_of("efss"))
}

#' @export
ens_fss <- function(
  x, y, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, ...
) {
  UseMethod("ens_fss")
}

#' @export
ens_fss.geolist <- function(
  x, y, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, ...
) {
  stopifnot(meteogrid::is.geofield(y))
  stopifnot(
    meteogrid::compare.geodomain(
      attr(x[[1]], "domain"),
      meteogrid::as.geodomain(y)
    )
  )

  quantile_thresh <- FALSE
  thresholds <- harpCore::parse_thresholds(threshold)
  if (thresholds[["quantiles"]]) {
    quantile_thresh <- TRUE
    threshold <- quantile(y, thresholds[["thresholds"]])
    thresholds_df <- tibble::tibble(
      quantile  = thresholds[["thresholds"]],
      threshold = threshold
    )
  }

  result <- dplyr::bind_rows(
    lapply(
      threshold,
      function(thresh) {
        dplyr::mutate(
          cpp_ens_fss(
            x, y, thresh, radius,
            comparator, include_low, include_high
          ),
          threshold = thresh
        )
      }
    )
  )
  if (quantile_thresh) {
    result <- dplyr::relocate(
      suppressMessages(dplyr::inner_join(result, thresholds_df)),
      dplyr::all_of("quantile"),
      .before = "threshold"
    )
  }

  summarise_fss(
    result, type = "ensemble", dx = attr(x[[1]], "domain")[["dx"]]
  )

}

#' @export
ens_fss.harp_ens_grid_df <- function(
  x, y, threshold, radius, comparator = "ge",
  include_low = TRUE, include_high = TRUE, ...
) {
  dplyr::mutate(
    dplyr::rowwise(x),
    dplyr::across(
      {{y}},
      ~list(
        ens_fss(
          as_geolist(as.list(dplyr::pick(dplyr::matches("_mbr[[:digit:]]+")))),
          .x,
          threshold,
          radius,
          comparator,
          include_low,
          include_high
        )
      ),
      .names = "fss"
    )
  ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::where(is_geolist)) %>%
    tidyr::unnest(dplyr::all_of("fss"))
}

unique_pairs <- function(x) {
  x <- unique(x)
  g <- function(i) {
    z <- setdiff(x, x[seq_len(i)])
    if (length(z)) {
      cbind(x[i], z, deparse.level = 0)
    }
  }
  do.call(rbind, lapply(seq_along(x), g))
}

summarise_fss <- function(result, type, dx) {
  by = intersect(c("quantile", "threshold", "scale"), colnames(result))
  dplyr::summarise(
    result,
    fbs     = sum(.data[["fbs"]]),
    fbs_ref = sum(.data[["fbs_ref"]]),
    fss     = mean(.data[["fss"]]),
    count   = dplyr::n(),
    .by     = by
  ) %>%
    dplyr::mutate(
      grid_length = dx,
      nbhd_length = (.data[["scale"]] * 2 + 1) * dx,
      type        = type
    ) %>%
    dplyr::rename(
      nbhd_radius = dplyr::all_of("scale")
    ) %>%
    dplyr::relocate(
      dplyr::all_of(c("grid_length", "nbhd_length", "type")),
      .after = "nbhd_radius"
    ) %>%
    tibble::as_tibble()
}
