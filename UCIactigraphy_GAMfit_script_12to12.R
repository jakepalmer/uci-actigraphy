##########
# Actigraphy GAM fitting
# Created by: Jake Palmer
# Email: jake.palmer@sydney.edu.au
# Last edit: 23.05.2020
##########

#----- INSTRUCTIONS
# See associated documentation file - UCIactigraphy_GAMfit_documentation.html

# INPUT:
# - Trimmed Actiware output (files ending with _trimmed.csv)
# OUTPUT:
# - Plot of smoothed model fit per file
# - Diagnostic plots per file
# - File with collated output variables for each file included in the analysis batch

#----- FUNCTIONS
# The functions are seperated into 3 broad steps:
# 1. data prep
# 2. model fitting
# 3. plotting
# There are a number of low-level functions that contribute to each step. These low-level functions
# are combined in a function to wrap each of the above steps. Therefore the structure of this script is:
# - low-level data prep functions
# - wrapper of low-level data prep functions
# - low-level model fitting functions
# - wrapper of low-level model fitting functions
# - low-level plotting functions
# - wrapper of low-level model fitting functions
# - setup for analysis
# - loop over each file running each step with wrapper functions

#----- SETUP

# Check if required packages are installed, load if so, install if not
check.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, library, character.only = TRUE)
}


# Create empty dataframe for output
create_output <- function(fname, cols) {
  df_out <- data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(df_out) <- cols
  suppressWarnings(
    if (!file.exists(fname)) {
      write.table(df_out, fname,
        sep = ",",
        col.names = TRUE, append = TRUE, row.names = TRUE
      )
    }
  )
  assign("df_out", df_out, envir = .GlobalEnv)
}


# Save plots to jpg easily
save_plot <- function(fname, plot) {
  ggsave(fname, plot = plot, device = "jpeg", width = 8, height = 6, units = "in")
}


#----- 1. Data prep

# Drop rows with missing data
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


# Drop missing and calculate % missing for output
handle_missing <- function(df) {
  perc_miss <- mean(is.na(df$Activity)) * 100
  assign("perc_miss", perc_miss, envir = .GlobalEnv)
  df <- completeFun(df, "Activity")
  return(df)
}


# Check the format of time variables, reformat if needed
timeCheck <- function(time_val, date.format = "%I:%M:%S %p") {
  tryCatch(!is.na(as.Date(time_val, date.format)),
    error = function(err) {
      FALSE
    }
  )
}


# Convert from H:M:S to numeric representation of time and
# shift time axis for model fitting. Value of min_time represents
# the earliest time on x-axis (e.g. if min_time = 3.00 then earliest
# time on x-axis will be 03:00).
time_HMS_to_numeric <- function(df, min_time, col) {
  df$hrs <- as.numeric(format(strptime(df[[col]], "%H:%M:%S"), format = "%H"))
  df$min <- as.numeric(format(strptime(df[[col]], "%H:%M:%S"), format = "%M"))
  df$ss <- as.numeric(format(strptime(df[[col]], "%H:%M:%S"), format = "%S"))
  df$timehour <- df$hrs + (df$min / 60)
  df$timehour <- round(df$timehour, 2)
  df$timehour_shifted <- ifelse(df$timehour < min_time, df$timehour + 24, df$timehour)
  ### NEW: Added to bring minimum time back to midday (0.00 numeric)
  df$timehour_shifted <- df$timehour_shifted - min_time
  return(df)
}


# Calculate log(Activity) and aggregate by mean
aggLog <- function(df) {
  df$LogActivity <- log10(df$Activity + 1)
  df <- df %>%
    group_by(timehour_shifted) %>%
    summarise(agglogActivity = mean(LogActivity))
  return(df)
}


# WRAPPER for PREP STEPS
prep_data <- function(df) {
  df <- handle_missing(df)
  if (timeCheck(df$Time[1])) {
    df$Time <- format(strptime(df$Time, "%I:%M:%S %p"), format = "%H:%M:%S")
  }
  df <- time_HMS_to_numeric(df, min_time = 12.0, "Time")
  df <- aggLog(df)
  return(df)
}


#----- 2. Model fitting


# Fit GAM
fit_model <- function(formula, df) {
  mod <- gam(formula, data = df, method = "REML")
  rsq <- summary(mod)
  rsq <- rsq$r.sq
  assign("rsq", rsq, envir = .GlobalEnv)
  return(mod)
}


# Calculate UP/DOWN slope values and times
calc_slope <- function(df, mod) {
  # n_deriv <- nrow(df) / 2
  n_deriv <- nrow(df)
  first_deriv <- gratia::derivatives(mod, type = "central", n = n_deriv, order = 1)
  ### NEW - split at 15.00
  UPdf <- subset(first_deriv, data > 15.0)
  DOWNdf <- subset(first_deriv, data <= 15.0)
  UP <- UPdf[which.max(UPdf$derivative), ]
  UPslope <- UP$derivative
  UPtime <- UP$data
  DOWN <- DOWNdf[which.min(DOWNdf$derivative), ]
  DOWNslope <- DOWN$derivative
  DOWNtime <- DOWN$data
  assign("UPslope", UPslope, envir = .GlobalEnv)
  assign("UPtime", UPtime, envir = .GlobalEnv)
  assign("DOWNslope", DOWNslope, envir = .GlobalEnv)
  assign("DOWNtime", DOWNtime, envir = .GlobalEnv)
}


# Calulate turning points relative to UP/DOWN slope times
calc_turn_points <- function(df, mod, UPtime, DOWNtime) {
  df$smoothed_vals <- predict(mod, newdata = df)
  # Calculate first peak
  UPdf <- subset(df, timehour_shifted > UPtime)
  first_turns <- pastecs::turnpoints(as.numeric(UPdf$smoothed_vals))
  first_peaks_time <- UPdf$timehour_shifted[extract(first_turns, no.tp = FALSE, peak = TRUE, pit = FALSE)]
  first_peak <- min(first_peaks_time)
  # Calculate last peak
  DOWNdf <- subset(df, timehour_shifted <= DOWNtime)
  last_turns <- pastecs::turnpoints(as.numeric(DOWNdf$smoothed_vals))
  last_peaks_time <- DOWNdf$timehour_shifted[extract(last_turns, no.tp = FALSE, peak = TRUE, pit = FALSE)]
  last_peak <- max(last_peaks_time)
  # Assign
  assign("first_peak", first_peak, envir = .GlobalEnv)
  assign("last_peak", last_peak, envir = .GlobalEnv)
}


# Write output
write_output <- function(outfile) {
  df_out <- rbind(df_out, c(
    subj, round(UPtime, 2), round(UPslope, 2), round(DOWNtime, 2),
    round(DOWNslope, 2), round(first_peak, 2), round(last_peak, 2),
    round(rsq, 2), round(perc_miss, 2)
  ))
  write.table(df_out, outfile, sep = ",", col.names = FALSE, append = TRUE, row.names = FALSE)
}


# WRAPPER for MODEL FITTING
modelling <- function(df) {
  mod <- fit_model(formula = agglogActivity ~ s(timehour_shifted), df = df_agg)
  calc_slope(df_agg, mod)
  calc_turn_points(df_agg, mod, UPtime, DOWNtime)
  write_output(outfile)
  p_check <- appraise(mod)
  p_check_fname <- paste0(subj, "_GAMcheck.jpg")
  save_plot(p_check_fname, p_check)
  assign("gam_mod", mod, envir = .GlobalEnv)
}


#----- 3. Plotting


# Function to return basic features of the plot
base_plot <- function(df, x, y, smoothed_vals) {
  p <- df %>%
    ggplot(aes(x = x, group = 1)) +
    geom_point(aes(y = y), size = 1, alpha = 0.5, na.rm = TRUE) +
    geom_line(aes(y = smoothed_vals), colour = "red") +
    theme_classic() +
    # xlab(paste0("Time (", round(min(x), 0), " to ", round(max(x), 0), ")")) +
    xlab("Time") +
    ylab("log(Activity)") +
    scale_x_continuous(breaks = c(0, 12, 24),
                         labels = paste0(c("12:00pm", "12:00am", "12:00pm")))
  return(p)
}


# Plot with shaded habitual sleep periods, builds on base plot
plot_w_sleeptimes <- function(df, x, y, smoothed_vals) {
  sleep_data <- read.csv("SPRiNT_actigraphy SLEEP.csv") %>% filter(ID == subj)
  sleep_data <- time_HMS_to_numeric(sleep_data, min_time = 12.0, "Sleep_Onset_Time")
  sleep_data <- sleep_data %>% rename(Sleep_Onset_Time_shifted = timehour_shifted)
  sleep_data <- time_HMS_to_numeric(sleep_data, 12.0, "Sleep_Offset_Time")
  sleep_data <- sleep_data %>% rename(Sleep_Offset_Time_shifted = timehour_shifted)
  sleep_onset <- sleep_data$Sleep_Onset_Time_shifted[1]
  sleep_offset <- sleep_data$Sleep_Offset_Time_shifted[1]
  p_base <- base_plot(df, x, y, smoothed_vals)
  p_sleep <- p_base +
    geom_vline(xintercept = sleep_onset, linetype = "dashed") +
    geom_vline(xintercept = sleep_offset, linetype = "dashed") +
    annotate(geom = "rect", xmin = sleep_onset, xmax = sleep_offset, ymin = -Inf, ymax = Inf, fill = "slategray1", alpha = 0.5)
  return(p_sleep)
}


# Plot with dashed lines representing first and last peaks
plot_w_peaks <- function(df, x, y, smoothed_vals) {
  p_base <- base_plot(df, x, y, smoothed_vals)
  p_peaks <- p_base +
    geom_vline(xintercept = first_peak, linetype = "dashed") +
    geom_vline(xintercept = last_peak, linetype = "dashed")
}


# Plot with dashed lines representing UP and DOWN slopes
plot_w_slopes <- function(df, x, y, smoothed_vals) {
  p_base <- base_plot(df, x, y, smoothed_vals)
  p_slopes <- p_base +
    geom_vline(xintercept = UPtime, linetype = "dashed") +
    geom_vline(xintercept = DOWNtime, linetype = "dashed")
}


# TO EDIT
# plot_w_dlmo <- function(df, x, y, smoothed_vals) {
#
# }


# WRAPPER for PLOTTING
plotting <- function(df) {
  x <- df$timehour_shifted
  y <- df$agglogActivity
  smoothed_vals <- predict(gam_mod, newdata = df)
  #-- Base plot
  p_base <- base_plot(df, x, y, smoothed_vals)
  p_base_name <- paste0(subj, "_GAMplot.jpg")
  save_plot(p_base_name, p_base)
  #-- Plot with UP/DOWN slopes
  p_slope <- plot_w_slopes(df, x, y, smoothed_vals)
  p_slope_name <- paste0(subj, "_GAMplot_slopes.jpg")
  save_plot(p_slope_name, p_slope)
  #-- Plot with turning points
  p_tp <- plot_w_peaks(df, x, y, smoothed_vals)
  p_tp_name <- paste0(subj, "_GAMplot_peaks.jpg")
  save_plot(p_tp_name, p_tp)
  #-- Plot with sleep times
  p_sleep <- plot_w_sleeptimes(df, x, y, smoothed_vals)
  p_sleep_name <- paste0(subj, "_GAMplot_sleeptimes.jpg")
  save_plot(p_sleep_name, p_sleep)
  #-- Plot with DLMO - TO EDIT
  # p_dlmo <- plot_w_sleeptimes(df, x, y, smoothed_vals)
  # p_dlmo_name <- paste0(subj, "_GAMplot_DLMO.jpg")
  # save_plot(p_dlmo_name, p_dlmo)
}


#----- RUN COMPLETE PROCESS
# Setup

setwd("/home/UCIactig/data")
packages <- c("dplyr", "lubridate", "ggplot2", "gratia", "mgcv", "pastecs")
check.packages(packages)

# Prepare output file

outfile <- "GAMfit_output.csv"
output_cols <- c(
  "Subject", "UPslope_time", "UPslope", "DOWNslope_time",
  "DOWNslope", "FIRSTpeak_time", "LASTpeak_time", "rsq",
  "percent_missing"
)
create_output(outfile, output_cols)

# Get list of files to be processed

file_list <- Sys.glob("*_trimmed.csv")

### Run process for all files ###

for (file in file_list) {
  print(paste0("Working on: ", file))
  subj <- strsplit(file, "_trimmed.csv")[[1]]
  df_subj <- read.csv(file)
  df_agg <- prep_data(df_subj)
  modelling(df_agg)
  plotting(df_agg)
}

print("SCRIPT RUN COMPLETE")
