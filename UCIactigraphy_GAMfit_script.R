##########
# Actigraphy GAM fitting
# Created by: Jake Palmer
# Email: jake.palmer@sydney.edu.au
# Last edit: 13.10.2019
##########

#####
# INSTRUCTIONS
#####
# See associated documentation file - UCIactigraphy_GAMfit_documentation.html

# INPUT:
# - Trimmed Actiware output (files ending with _trimmed.csv)
# OUTPUT:
# - Plot of smoothed model fit per file
# - Diagnostic plots per file
# - File with collated output variables for each file included in the analysis batch

#####
# GAM FITTING
#####

setwd('/Users/mq44848301/Desktop')

check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, library, character.only = TRUE)
}
packages <- c("dplyr", "lubridate", "ggplot2", "gratia", "mgcv")
check.packages(packages)

# Create empty dataframe for output
df_out <- data.frame(matrix(ncol = 7, nrow = 0))
cols <- c("Subject", "UPslope_time", "UPslope", "DOWNslope_time",
          "DOWNslope", "rsq", "percent_missing")
colnames(df_out) <- cols
suppressWarnings(write.table(df_out, "GAMfit_output.csv", sep = ",",
            col.names = TRUE, append = TRUE,
            row.names = FALSE))

# Function to be used to drop missing data
completeFun <- function(data, desiredCols) {
    completeVec <- complete.cases(data[, desiredCols])
    return(data[completeVec, ])
}

timeCheck <- function(time_val, date.format = "%I:%M:%S %p") {
    tryCatch(!is.na(as.Date(time_val, date.format)),  
             error = function(err) {FALSE})  
}

# Get list of files to be processed
file_list <- Sys.glob("*_trimmed.csv")

# Main routine for fitting GAM
print("Aggregating on MEAN:")
fit <- lapply(file_list, function (file) {
    print(paste0("    Working on: ", file))
    subj <- strsplit(file, '_trimmed.csv')[[1]]
    df_subj <- read.csv(file)
    perc_miss <- mean(is.na(df_subj$Activity)) * 100
    df_subj <- completeFun(df_subj, 'Activity')
    time_val <- df_subj$Time[1]
    if (timeCheck(time_val)) {
        df_subj$Time <- format(strptime(df_subj$Time, "%I:%M:%S %p"), format = "%H:%M:%S")
    }
    df_subj$hrs <- as.numeric(format(strptime(df_subj$Time, "%H:%M:%S"), format = "%H"))
    df_subj$min <- as.numeric(format(strptime(df_subj$Time, "%H:%M:%S"), format = "%M"))
    df_subj$ss <- as.numeric(format(strptime(df_subj$Time, "%H:%M:%S"), format = "%S"))
    df_subj$timehour <- df_subj$hrs + (df_subj$min / 60)
    df_subj$timehour <- round(df_subj$timehour, 2)
    df_subj$timehour_shifted <- ifelse(df_subj$timehour < 3.0, df_subj$timehour + 24, df_subj$timehour)
    df_subj$LogActivity <- log10(df_subj$Activity + 1)
    df_agg <- df_subj %>% 
        group_by(timehour_shifted) %>%
        summarise(agglogActivity = mean(LogActivity))
    gam_mod <- gam(agglogActivity ~ s(timehour_shifted), data = df_agg, method = "REML")
    n_deriv <- nrow(df_agg) / 2
    deriv <- derivatives(gam_mod, type = "central", n = n_deriv)
    max_slope <- deriv[which.max(deriv$derivative),]
    max_slope_slope <- max_slope$derivative
    max_slope_time <- max_slope$data
    min_slope <- deriv[which.min(deriv$derivative),]
    min_slope_slope <- min_slope$derivative
    min_slope_time <- min_slope$data
    rsq <- summary(gam_mod)
    rsq <- rsq$r.sq
    
    # Write output
    df_out <- df_out %>% add_row("Subject" = subj, "UPslope_time" = round(max_slope_time, 2),
                                 "UPslope" = round(max_slope_slope, 2), "DOWNslope_time" = round(min_slope_time, 2),
                                 "DOWNslope" = round(min_slope_slope, 2), "rsq" = round(rsq, 2), 
                                 "percent_missing" = round(perc_miss, 2))
    write.table(df_out, "GAMfit_output.csv", sep = ",",
                col.names = FALSE, append = TRUE, row.names = FALSE)
    
    # Checking and plotting
    rsqstr <- toString(round(rsq, digits = 2))
    perc_missstr <- toString(round(perc_miss, digits = 2))
    p_check <- appraise(gam_mod)
    ggsave(paste0(subj, '_GAMcheck.jpg'), plot = p_check, device = "jpeg", width = 8,
           height = 6, units = 'in')
    pred_data <- data.frame(time = df_agg$timehour_shifted,
                            activity = df_agg$agglogActivity,
                            predicted_values = predict(gam_mod, newdata = df_agg))
    p <- ggplot(pred_data, aes(x = time)) + 
        geom_point(aes(y = activity), size = 1, alpha = 0.5, na.rm=TRUE) +
        geom_line(aes(y = predicted_values), colour = "red") +
        geom_vline(xintercept = max_slope_time, linetype="dashed") +
        geom_vline(xintercept = min_slope_time, linetype="dashed") +
        xlab('24hr Time (03:00 to 27:00 for model fitting)') +
        ylab('log(Activity)') +
        theme_classic() +
        labs(title = paste0(subj),
             subtitle = paste0('R2 = ',rsqstr,'; Percent activity missing = ',perc_missstr,'; Aggregated on mean'),
             caption = "Note: Dashed verticle lines depict calculated time of steepest UP and DOWN slope") +
        theme(plot.caption = element_text(hjust = 0))
    ggsave(paste0(subj, '_GAMplot.jpg'), plot = p, device = "jpeg", width = 8,
           height = 6, units = 'in')
})
print('SCRIPT RUN COMPLETE')