#' Run spatial verification on a (for now) deterministic forecast
#'
#' @param start_date Date of the first forecast to read.
#' @param end_date Date of the last forecast to read.
#' @param model The name of the (deterministic or EPS) model.
#' @param parameter The parameters to read as a character vector.
#' @param lead_time The lead times to read as a numeric vector.
#'   Should be in the units that are also used in fc_file_template.
#' @param lt_unit The unit used for lead_time. Can be "h" (hours), "m" (minutes), "s" (seconds)
#' @param by The time between forecasts. Should be a string of a number followed
#'   by a letter, where the letter gives the units - may be "d" for days, "h" for
#'   hours or "m" for minutes.
#' @param members The (numbers of the) ensemble members to read. While Netcdf and grib2 files 
#'   can contain multiple members, for other formats we assume they are in separate files
#'   (see also fc_file_template)
#' @param fc_file_path The top level path for the forecast files to read.
#' @param fc_file_template The file type to generate the template for. Can be
#'   "harmoneps_grib", "harmeoneps_grib_fp", "harmoneps_grib_sfx", "meps_met",
#'   "harmonie_grib", "harmonie_grib_fp", "harmone_grib_sfx", "vfld", "vobs", or
#'   "fctable". If anything else is passed, it is returned unmodified. In this
#'   case substitutions can be used. Available substitutions are \{YYYY\} for
#'   year, \{MM\} for 2 digit month with leading zero, \{M\} for month with no
#'   leading zero, and similarly \{DD\} or \{D\} for day, \{HH\} or \{H\} for
#'   hour, \{mm\} or \{m\} for minute. Also \{LDTx\} for lead time and \{MBRx\}
#'   for ensemble member where x is the length of the string including leading
#'   zeros - can be omitted or 2, 3 or 4. Note that the full path to the file
#'   will always be file_path/template.
#' @param fc_file_format The format of the files to read. Can be e.g. "fa" or "grib".
#' @param fc_file_options A list with format-specific options for the reader function.
#' @param fc_domain The forecast domain. If provided, the fc reading can be made faster
#'    by not extracting domain information (format option meta).
#' @param fc_interp_method Interpolation method to be used when transforming a forecast
#'   field to the verification grid.
#' @param fc_accumulation The accumulation type of the forecast. This is only used for
#'   accumulated parameters (e.g. precipitation). NULL signifies that the field is accumulated
#'   from the start of the model run. Otherwise this should be a string containing a numerical value
#'   and a time unit, e.g. "15m" or "1h".
#' @param ob_file_path The top level path for the forecast files to read.
#' @param ob_file_template The file type to generate the template for. Can be
#'   "harmoneps_grib", "harmeoneps_grib_fp", "harmoneps_grib_sfx", "meps_met",
#'   "harmonie_grib", "harmonie_grib_fp", "harmone_grib_sfx", "vfld", "vobs", or
#'   "fctable". If anything else is passed, it is returned unmodified. In this
#'   case substitutions can be used. Available substitutions are {YYYY} for
#'   year, \{MM\} for 2 digit month with leading zero, \{M\} for month with no
#'   leading zero, and similarly \{DD\} or \{D\} for day, \{HH\} or \{H\} for
#'   hour, \{mm\} or \{m\} for minute. Also \{LDTx\} for lead time and \{MBRx\}
#'   for ensemble member where x is the length of the string including leading
#'   zeros - can be omitted or 2, 3 or 4. Note that the full path to the file
#'   will always be file_path/template.
#' @param ob_file_format The format of the files to read. Can be e.g. "hdf5" or "grib".
#' @param ob_file_options A list with format-specific options for the reader function.
#' @param ob_domain The observation domain. If provided, the obs reading can be made faster
#'    by not extracting domain information (format option meta).
#' @param ob_interp_method Interpolation method to be used when transforming a forecast
#'   field to the verification grid.
#' @param ob_accumulation The accumulation type of the observation (or reference). This is only used for
#'   accumulated parameters (e.g. precipitation). NULL signifies that the field is accumulated
#'   from the start of the model run. That is probably rare for observations. 
#'   Otherwise this should be a string containing a numerical value
#'   and a time unit, e.g. "15m" or "1h". "0" means an instantaneous value.
#' @param verif_domain A \code{geodomain} that defines the common verification grid.
#' @param return_data          = TRUE,
#' @param thresholds Thresholds used for FSS, ...
#' @param window_sizes Scales used for fuzzy methods like FSS. A vector of box sizes.
#'   All values must be odd integers (so the central point is really in the center of a box).
#' @param sqlite_path If specified, SQLite files are generated and written to
#'   this directory.
#' @param sqlite_file Name of SQLite file.
#' @param return_data If TRUE, the result is returned as a list of tables.
#' @param ... Not used at thispoint (more info to be added).
#'
#' @return A list containting tibbles for all scores.
#' @export

verify_spatial <- function(start_date,
                           end_date,
                           parameter,
                           model                = harpSpatial_conf$model,
                           lead_time            = harpSpatial_conf$lead_time, # seq(0,36,3)
                           lt_unit              = harpSpatial_conf$lt_unit, #"h",
                           by                   = harpSpatial_conf$by, # "12h",
                           scores               = NULL,
                           members              = harpSpatial_conf$members, #NULL,
#                           members_out          = members,
#                           lags                 = harpSpatial_conf$lags, #NULL,
                           fc_file_path         = harpSpatial_conf$fc_file_path, # "",
                           fc_file_template     = harpSpatial_conf$fc_file_template, #"",
                           fc_file_format       = harpSpatial_conf$fc_file_format, #"fa",
                           fc_file_options      = harpSpatial_conf$fc_file_options, #list(),
                           fc_domain            = harpSpatial_conf$fc_domain, #NULL,
                           fc_interp_method     = harpSpatial_conf$fc_interp_method, #"closest",
                           fc_accumulation      = harpSpatial_conf$fc_accumulation, #NULL,
                           ob_file_path         = harpSpatial_conf$ob_file_path, #"",
                           ob_file_template     = harpSpatial_conf$ob_file_template, #"",
                           ob_file_format       = harpSpatial_conf$ob_file_format, #"hdf5",
                           ob_file_options      = harpSpatial_conf$ob_file_options, #list(),
                           ob_domain            = harpSpatial_conf$ob_domain, #NULL,
                           ob_interp_method     = harpSpatial_conf$ob_interp_method, #"closest",
                           ob_accumulation      = harpSpatial_conf$ob_accumulation, #"15m",
                           verif_domain         = harpSpatial_conf$verif_domain, #NULL,
                           use_mask             = harpSpatial_conf$use_mask, #FALSE,
                           window_sizes         = harpSpatial_conf$window_sizes, #c(1, 3, 5, 11, 21),
                           thresholds           = harpSpatial_conf$thresholds, #c(0.1, 1, 5, 10),
                           sqlite_path          = harpSpatial_conf$sqlite_path, #NULL,
                           sqlite_file          = harpSpatial_conf$sqlite_file, #"harp_spatial_scores.sqlite",
                           return_data          = FALSE) {

  # TODO: we may need more options! masked interpolation, options by score, 
  prm <- harpIO::parse_harp_parameter(parameter)
  by_secs <- harpIO:::units_multiplier(by) * readr::parse_number(by)

  # For efficiency, we use a slightly counter-intuitive loop order
  # we don't loop over forecast date and then lead time,
  # because that would cause excessive re-reading (or caching) of observations.
  # Rather, we loop over all observation times 
  # and then over all forecasts valid for those times.
  # Close to start_date and end_date you must make sure not to read beyond the time window.
  # TODO: to be even more efficient, we could try to
  #       - open (&parse) FC fields only once
  #       - read accumulated fields only once (e.g. "acc3h = 6h - 3h" also re-uses 
  #                                            the 3h field) 
  #       But that would require extensive "caching", which may end up even slower.
  # Alternative strategy: loop by fcdate, and cache all obs in a list by leadtime
  #     next fcdate -> "ldt -= by" ; drop negative ldt ; read missing obs
  # For an accumulated variable (precip), the minimum lead time is 
  #     the accumulation time. Otherwise zero.

  # some date handling first : create vectors of date_time class
  sdate <- lubridate::ymd_hm(start_date, tz = "UTC", truncated = 4)
  if (missing(end_date)) {
    edate <- sdate
  } else {
    edate <- lubridate::ymd_hm(end_date, tz = "UTC", truncated = 4)
    # add the last fc hour if the date is just YYYYMMDD
    if (nchar(end_date) < 10) edate <- edate + (24 * 3600 / by_secs - 1) * by_secs
  }
  # convert lead_time to seconds and remove lead_times smaller than accum
  # we don't have 3h precip at 0h forecast.
  # also, we probably want lead_time in steps of the accumulation
  lt_scale <- harpIO:::units_multiplier(lt_unit)
  lead_time <- lead_time * lt_scale
  if (prm$accum > 0) {
    lead_time <- lead_time[which(lead_time >= prm$accum & lead_time %% prm$accum == 0)]
  }
  all_fc_dates <- 
    seq(as.numeric(sdate), as.numeric(edate), by_secs) %>% 
      lubridate::as_datetime()
  all_ob_dates <- (rep(all_fc_dates, each=length(lead_time)) + lead_time ) %>%
    unique() %>%
    sort()

  message("Running spatial verification.")
  message("Forecast dates: ", paste(all_fc_dates, collapse = " "))
  message("Lead times: ", paste(lead_time / lt_scale, collapse = " "))
  message("Observation dates: ", paste(all_ob_dates, collapse = " ; "))
  # the re-gridding weights will come here:
  init <- list()

  # The read functions
  # Most read functions can't deal with accumulated parameters like AccPcp1h
  # We will need special "accumulator functions"

  # FIXME: also, we must MODIFY the parameter!
  # - correct accumulation
  # - maybe even a different field -> need a "modifier"???
  # if (!is.null(ob_param$accum)
  ob_param <- prm
  ob_param$accum <- readr::parse_number(ob_accumulation) * 
                    harpIO:::units_multiplier(ob_accumulation)
  # FIXME: avoid reading domain information for every file (obs and fc)
  #        BUT: we need it once to initialise the regridding. Use "get_domain(file)".
  # FIXME: should we do the regridding within the read_grid call?

  get_ob <- function(obdate) {
    obfile <- get_filenames(
      file_date     = format(obdate, "%Y%m%d%H"),
      file_path     = ob_file_path,
      file_template = ob_file_template,
      parameter     = ob_param
    )
    do.call(harpIO::read_grid, 
      c(list(file_name=obfile, file_format=ob_file_format,
                   parameter = ob_param, file_format_opts = ob_file_options)))
  }

  # FIXME: if (!is.null(members) && length(members) > 1)
  # NOTE: lead_time in original units
  # FIXME: det_model vs eps_model...
  # FIXME: for eps, either read_grid on single netcdf/grib2 or on multiple FA/GRIB
  # if "members" is defined (and length > 1) but {MBRx} is not in the template -> single file
  if (is.null(members)) {
    get_fc <- function(fcdate, lead_time) {
      fcfile <- get_filenames(
        file_date = fcdate,
        lead_time = lead_time,
        parameter = parameter,
        det_model = model,
        file_path = fc_file_path,
        file_template = fc_file_template)
      do.call(harpIO::read_grid,
        c(list(file_name = fcfile, file_format = fc_file_format, 
                   parameter = parameter, lead_time = lead_time,
                       file_format_opts = fc_file_options)))
    }
  } else {
    # for EPS models, we try to get the members into a 3D geogrid array.
    # very fast for passing to Rccp functions.
    get_fc <- function(fcdate, lead_time) {
      fcfile <- get_filenames(
        file_date = fcdate,
        lead_time = lead_time,
        parameter = parameter,
        eps_model = model,
        members   = members,
        file_path = fc_file_path,
        file_template = fc_file_template)
      if (length(fcfile) == 1) {
        do.call(harpIO::read_grid,
                c(list(file_name = fcfile, file_format = fc_file_format, 
                  parameter = parameter, lead_time = lead_time,
		  file_format_opts = fc_file_options)))
      } else {
        lapply(fcfile, harpIO::read_grid, file_format = fc_file_format, 
                parameter = parameter, lead_time = lead_time, members=members,
                unlist(fc_file_options))
      }
    }
  }
  # We will write to SQL only at the end (more efficient),
  ncases <- length(all_fc_dates) * length(lead_time)
  message("expected ncases= ", ncases)

  if (is.null(scores)) {
    score_list <- spatial_scores()
  } else {
    score_list <- spatial_scores()[scores]
  }
#  score_templates <- lapply(names(score_list), function(sc) spatial_scores(score = sc))
#  names(score_templates) <- score_list

  # Define the list of score tables.
  score_tables <- vector("list", length(score_list))
  names(score_tables) <- names(score_list)
  # some score funtions calculate several scores together
  # we don't want to call them twice...
  score_function_list <- unique(sapply(score_list, function(x) x$func))
  # And conversely, for every such "multiscore", we need the list of scores that depend on it
  score_function_subset <- sapply(score_function_list, function(msc)
                                  names(which(sapply(score_list, function(x) x$func == msc))))


  message("score functions: ", paste(score_function_list, collapse=" "))
  # MAIN LOOP
  case <- 1
  for (ob in seq_along(all_ob_dates)) {  # (obdate in all_ob_dates) looses POSIXct class
    obdate <- all_ob_dates[ob]
    message("=====\nobdate: ", format(obdate, "%Y%m%d-%H%M"))
    obfield <- get_ob(obdate)
    if (inherits(obfield, "try-error")) { # e.g. missing observation
      if (harpenv$verbose) message("Observation not found. Skipping.\n")
      next
    }
    if (prm$accum > 0) { # an accumulated field like precipitation
      if (is.null(ob_accumulation) || ob_accumulation < 0) {
        # RARE: observation is an accumulated field (e.g. reference run)
        # TODO: This does not look very useful, unless get_ob() can somehow deal with it.
        warning("Accumulated observation fields not yet validated.", immediate.=TRUE)
        obfield <- obfield - get_ob(obdate - prm$accum)
      } else {
        ostep <- readr::parse_number(ob_accumulation) * harpIO:::units_multiplier(ob_accumulation)
        if (ostep == prm$accum) { # this is easy !
          # nothing to do
        } else if  (ostep > prm$accum) { # this is easy !
          stop("The chosen accumulation time is smaller than that of the observations!")
        } else {
          nstep <- prm$accum / ostep
          for (i in 1:(nstep-1)) {
            obfield <- obfield + get_ob(obdate - i*ostep)
          }
        }
      }
    }

    # convert to common verification grid
    if (!is.null(ob_interp_method)) {
      if (is.null(init$regrid_ob)) {
        message("Initialising observation regridding.")
        if (is.null(ob_domain)) {
          if (!is.null(ob_file_options$domain)) {
            ob_domain <- ob_file_options$domain
          } else {
            ob_domain <- obfield
          }
        }
        init$regrid_ob <- meteogrid::regrid.init(
          olddomain = ob_domain,
          newdomain = verif_domain,
          method    = ob_interp_method)
      }
      obfield <- meteogrid::regrid(obfield, weights = init$regrid_ob)
    }

    # find forecasts valid for this date/time
    # intersect drops the POSIXct class
    # valid_fc_dates <- intersect(obdate - lead_time, all_fc_dates)
    valid_fc_dates <- (obdate - lead_time)[which((obdate - lead_time) %in% all_fc_dates)]
    message("valid FC dates: ", paste(format(valid_fc_dates, "%Y%m%d-%H%M"), collapse = " "))
    # inner loop
    for (fc in seq_along(valid_fc_dates)) {
      fcdate <- valid_fc_dates[fc]
      ldt <- (as.numeric(obdate) - as.numeric(fcdate)) # in seconds !
      message(
        "   +++ fcdate = ", format(fcdate,"%Y%m%d-%H%M"),
        " +++ ldt = ", ldt / lt_scale, lt_unit
      )

      fcfield <- get_fc(fcdate, ldt/lt_scale)
      if (inherits(fcfield, "try-error")) { # e.g. missing forecast run
        if (harpenv$verbose) message("..... Forecast not found. Skipping.",
                                     .immediate = TRUE)
        next
      }
      if (prm$accum > 0) {
        if (is.null(fc_accumulation) || fc_accumulation < 0) {
          if (ldt > prm$accum) { # if ldt==accum, you don't need to decumulate
            fcfield <- fcfield - get_fc(fcdate, (ldt - prm$accum) / lt_scale)
          }
        } else {
          # In rare cases the forecast model needs "accumulating" rather than "decumulating"
          #       e.g. when verifying INCA against radar
          fstep <- readr::parse_number(fc_accumulation) * 
            harpIO:::units_multiplier(fc_accumulation)
          if (fstep == prm$accum) { # this is easy !
          # nothing to do
          } else if  (fstep > prm$accum) { # this is easy !
            stop("The chosen accumulation time is smaller than what is available in the forecasts!")
          } else {
            nstep <- prm$accum / fstep
            for (i in 1:(nstep - 1)) {
              fcfield <- fcfield + get_fc(fcdate, (ldt - i * fstep) / lt_scale)
            }
          }
        }
      }
      # convert forecast to common verification grid
      if (!is.null(fc_interp_method)) {
        if (is.null(init$regrid_fc)) {
          message("Initialising fc regridding.")
          if (is.null(fc_domain)) {
            if (!is.null(fc_file_options$domain)) {
              fc_domain <- fc_file_options$domain
            } else {
              fc_domain <- fcfield
            }
          }
          init$regrid_fc <- 
            meteogrid::regrid.init(
              olddomain = fc_domain,
              newdomain = verif_domain,
              method = fc_interp_method
            )
        }
        fcfield <- meteogrid::regrid(fcfield, weights = init$regrid_fc)
      }

      #############################
      ### NOW WE COMPUTE SCORES ###
      #############################

      # FIXME: Some scores (e.g. SAL) have various other parameters that we can't pass yet...
      #        While others don't need any
      for (sf in score_function_list) {
        # get the required arguments for this function
        # NOTE: args() only works if the function is found
        #       so either exported, or only internal use
#        arglist <- names(as.list(args(sf)))
        myargs <- list(obfield=obfield, fcfield=fcfield,
                         thresholds = thresholds,
                         scales = window_sizes)
        message("--> Calling ", sf)
        multiscore <- do.call(sf, myargs)

        if (!is.null(multiscore)) {
          # nrows per case for this score
          message("output dim : ", paste(dim(multiscore), collapse="x"))
          nrow <- dim(multiscore)[1]
          # interval of rows for this case in full score table
          intv <- seq_len(nrow) + (case - 1) * nrow
          for (sn in score_function_subset[[sf]]) {
            message("-----> Calling score ", sn)
            sc <- multiscore[,c(score_list[[sn]]$primary, score_list[[sn]]$fields)]
            if (is.null(score_tables[[sn]])) {
              template <- spatial_score_table(score_list[[sn]])
              tbl_struct <- lapply(template$fields, 
                                   function(x) switch(x, 
                                             "CHARACTER" = NA_character_,
                                             "INTEGER"   = NA_integer_, 
                                             "REAL"      = NA_real_,
                                             NA_real_))
              score_tables[[sn]] <- do.call(tibble::tibble, c(tbl_struct, .rows = ncases * nrow))
              # we can already fill some constant columns
              if ("model" %in% names(score_tables[[sn]])) score_tables[[sn]]$model <- model
              if ("prm" %in% names(score_tables[[sn]]))   score_tables[[sn]]$prm   <- parameter
            }
            # which interval of the score table is to be filled (may be only 1 row -> score[case, ...])
            # NOTE: save fcdate as unix date, leadtime in seconds !
            score_tables[[sn]]$fcdate[intv] <- as.numeric(fcdate)
            score_tables[[sn]]$leadtime[intv] <- ldt
            score_tables[[sn]][intv, names(sc)] <- sc
          }
        }
      }

      case <- case + 1
    } # fcdate
  } #obdate
  if (case < ncases + 1) {
    message("There were ", ncases + 1 - case, "missing cases out of ", ncases, ".")
    ncases <- case - 1
  }

  ## write to SQLite
  if (!is.null(sqlite_file)) {
    save_spatial_verif(score_tables, sqlite_path, sqlite_file)
  }

  if (return_data) invisible(score_tables)
  else invisible(NULL)
}

#' Save spatial scores to SQLite
#' @param scores A list of spatial score tables
#' @param sqlite_path The path for the sqlite file
#' @param sqlite_file The file to which the tables are written or added
save_spatial_verif <- function(score_tables, sqlite_path, sqlite_file) {
  if (is.null(sqlite_path)) db_file <- sqlite_file
  else db_file <- paste(sqlite_path, sqlite_file, sep = "/")
  message("Writing to SQLite file ", db_file)
  db <- harpIO:::dbopen(db_file)
  for (sc in names(score_tables)) {
    # check for score table and create if necessary
    tab <- spatial_score_table(spatial_scores(score = sc))
    # drop empty rows (missing cases)
    harpIO:::create_table(db, sc, tab$fields, tab$primary)
    # TODO: should we drop all cases were any field is missing?
    ok <- !is.na(score_tables[[sc]][, "fcdate"])
    message(sc, ":", dim(score_tables[[sc]])[1], " rows, ", sum(ok), " non-NA.")
    harpIO:::dbwrite(db, sc, score_tables[[sc]][ok, ])
  }

  harpIO:::dbclose(db)
}

