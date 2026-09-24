#' Pharmacokinetic Dosing Models for Ribavirin and Favipiravir
#'
#' Functions to compute dosing events and return concentration-time profiles
#' for the two antivirals. All time is in **days**.
#'
#' @name pk_models
NULL

#' Ribavirin dosing schedule (HFRS trial regimen)
#'
#' IV loading dose 33 mg/kg, then 16 mg/kg q6h x 4 days,
#' then 8 mg/kg q8h x 6 days (total 10 days).
#'
#' @param dose_mg Loading dose in mg (default 33 mg/kg * body_weight_kg)
#' @param body_weight_kg Patient body weight in kg
#' @param start_day Day of treatment initiation (relative to symptom onset)
#' @param duration_days Total treatment duration (default 10)
#' @return Data frame with columns: time (days), dose_mg, route
#' @export
ribavirin_dosing_schedule <- function(dose_mg = NULL,
                                       body_weight_kg = 75,
                                       start_day = 3,
                                       duration_days = 10) {
  if (is.null(dose_mg)) {
    dose_mg <- 33 * body_weight_kg  # loading dose
  }

  events <- data.frame(
    time    = numeric(),
    dose_mg = numeric(),
    route   = character(),
    stringsAsFactors = FALSE
  )

  # Loading dose at start
  events <- rbind(events, data.frame(
    time = start_day, dose_mg = dose_mg, route = "IV"
  ))

  # Phase 1: 16 mg/kg q6h for 4 days, beginning 6 h after loading.
  dose_phase1 <- 16 * body_weight_kg
  for (day in 0:3) {
    for (hour in c(6, 12, 18, 24)) {
      t <- start_day + day + hour / 24
      if (t <= start_day + duration_days) {
        events <- rbind(events, data.frame(
          time = t, dose_mg = dose_phase1, route = "IV"
        ))
      }
    }
  }

  # Phase 2: 8 mg/kg q8h (3 times/day) from day 5 to the end of the course.
  # The loop bound is derived from duration_days so a longer course actually
  # emits more doses; at the default duration_days = 10 this is 4:9 exactly as
  # before, so existing results are unchanged.
  dose_phase2 <- 8 * body_weight_kg
  for (day in seq.int(4, max(4, duration_days - 1))) {
    for (hour in c(8, 16, 24)) {
      t <- start_day + day + hour / 24
      if (t <= start_day + duration_days) {
        events <- rbind(events, data.frame(
          time = t, dose_mg = dose_phase2, route = "IV"
        ))
      }
    }
  }

  events <- events[order(events$time), ]
  rownames(events) <- NULL
  events
}

#' Favipiravir dosing schedule
#'
#' Standard: 1600 mg BID day 1, then 600 mg BID x 14 days.
#' High dose: 2400 mg BID day 1, then 1200 mg BID x 14 days.
#'
#' @param regimen "standard" or "high"
#' @param start_day Day of treatment initiation
#' @param duration_days Total course length in days (default 15)
#' @return Data frame with columns: time (days), dose_mg, route
#' @export
favipiravir_dosing_schedule <- function(regimen = "standard",
                                         start_day = 3,
                                         duration_days = 15) {
  if (regimen == "standard") {
    dose_day1 <- 1600
    dose_maint <- 600
  } else if (regimen == "high") {
    dose_day1 <- 2400
    dose_maint <- 1200
  } else {
    stop("regimen must be 'standard' or 'high'")
  }

  events <- data.frame(
    time    = numeric(),
    dose_mg = numeric(),
    route   = character(),
    stringsAsFactors = FALSE
  )

  # Day 1: BID (every 12h)
  for (hour in c(0, 12)) {
    events <- rbind(events, data.frame(
      time = start_day + hour / 24, dose_mg = dose_day1, route = "PO"
    ))
  }

  # Days 2 onward: maintenance BID, to the end of the course. At the default
  # duration_days = 15 this is 1:14 exactly as before.
  for (day in seq_len(max(1, duration_days - 1))) {
    for (hour in c(0, 12)) {
      events <- rbind(events, data.frame(
        time = start_day + day + hour / 24, dose_mg = dose_maint, route = "PO"
      ))
    }
  }

  events <- events[order(events$time), ]
  rownames(events) <- NULL
  events
}

#' Build IV bolus events for deSolve
#'
#' Converts a dosing schedule into an events data frame compatible with
#' `deSolve::lsoda()` via the `events` argument.
#'
#' @param schedule Data frame from [ribavirin_dosing_schedule()] or
#'   [favipiravir_dosing_schedule()]
#' @param comp_name Compartment name in the state vector
#' @param body_weight_kg Patient body weight (for dose scaling)
#' @param Vd Volume of distribution in L (to convert mg -> ug/mL)
#' @return Data frame with columns: time, var, value, method
#' @export
build_deSolve_events <- function(schedule, comp_name,
                                  body_weight_kg = 75,
                                  Vd = 45) {
  if (nrow(schedule) == 0) {
    return(data.frame(
      time = numeric(), var = character(),
      value = numeric(), method = character(),
      stringsAsFactors = FALSE
    ))
  }

  # Convert IV doses to concentration increments. The favipiravir gut
  # compartment is an amount compartment (mg), not a concentration.
  delta_C <- if (comp_name == "C_FAV_gut") {
    schedule$dose_mg
  } else {
    # dose_mg / Vd_L = mg/L = ug/mL
    schedule$dose_mg / Vd
  }

  data.frame(
    time   = schedule$time,
    var    = comp_name,
    value  = delta_C,
    method = "add",
    stringsAsFactors = FALSE
  )
}

#' Combine multiple event schedules
#'
#' @param ... Named event data frames from [build_deSolve_events()]
#' @return Combined events data frame sorted by time
#' @export
combine_events <- function(...) {
  evts <- do.call(rbind, list(...))
  evts <- evts[order(evts$time), ]
  rownames(evts) <- NULL
  evts
}

#' Renal function-scaled Ribavirin clearance
#'
#' @param CL_RBV_std Standard clearance at eGFR = 90 mL/min (L/day)
#' @param eGFR Patient eGFR in mL/min
#' @param eGFR_ref Reference eGFR (default 90)
#' @return Scaled clearance in L/day
#' @export
scale_ribavirin_clearance <- function(CL_RBV_std, eGFR, eGFR_ref = 90) {
  CL_RBV_std * (eGFR / eGFR_ref)
}
