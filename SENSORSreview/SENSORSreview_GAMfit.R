library(dplyr)
library(lubridate)
library(ggplot2)
library(gratia)
library(mgcv)
library(sjlabelled)

setwd('/Users/mq44848301/Desktop/SENSORSreview/data/processed')

fname <- 'EbyE post_include+group 2019.07.13.sav'
timepoint <- 'post'

# Create empty dataframe for output
df_out <- data.frame(matrix(ncol = 9, nrow = 0))
cols <- c("Subject", "group", "UPslope_time", "UPslope", "DOWNslope_time",
          "DOWNslope", "rsq", "n_days", "percent_missing")
colnames(df_out) <- cols
write.table(df_out, "Activity_Slope_GAMoutput.csv", sep = ",",
            col.names = !file.exists("Activity_Slope_GAM_output.csv"), append = T,
            row.names = FALSE)

# Function to be used to drop missing data
completeFun <- function(data, desiredCols) {
    completeVec <- complete.cases(data[, desiredCols])
    return(data[completeVec, ])
}

# Read in complete data set.
df_all <- read_spss(fname)

# Get subject list from file
subj_list <- pull(df_all, ID)
subj_list <- unique(subj_list)

# Main routine for fitting GAM
print("Aggregating on MEAN:")
fit <- lapply(subj_list, function (subj) {
    print(paste0("Working on: ", subj))
    df_subj <- df_all %>% filter(ID == subj)
    colidx <- grep("Nights", colnames(df_all))
    n_days <- df_subj[1, colidx]
    group <- df_subj$Group[1]
    perc_miss <- mean(is.na(df_subj$LogActivity)) * 100
    df_subj <- completeFun(df_subj, 'LogActivity')
    df_subj$timehour_shifted <- ifelse(df_subj$timehour < 3.0, df_subj$timehour + 24, df_subj$timehour)
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
    df_out <- df_out %>% add_row("Subject" = subj, "group" = group, "UPslope_time" = max_slope_time,
                                 "UPslope" = max_slope_slope, "DOWNslope_time" = min_slope_time,
                                 "DOWNslope" = min_slope_slope, "rsq" = rsq, "n_days" = n_days,
                                 "percent_missing" = perc_miss)
    write.table(df_out, "Activity_Slope_GAMoutput.csv", sep = ",",
                col.names = !file.exists("Activity_Slope_GAMoutput.csv"),
                append = T, row.names = FALSE)
    
    # Checking and plotting
    rsqstr <- toString(round(rsq, digits = 2))
    perc_missstr <- toString(round(perc_miss, digits = 2))
    p_check <- appraise(gam_mod)
    ggsave(paste0(subj,'_',timepoint,'_GAMcheck.jpg'), plot = p_check, device = "jpeg", width = 8,
          height = 6, units = 'in')
    pred_data <- data.frame(time = df_agg$timehour_shifted,
                       activity = df_agg$agglogActivity,
                       predicted_values = predict(gam_mod, newdata = df_agg))
    p <- ggplot(pred_data, aes(x = time)) + 
        geom_point(aes(y = activity), size = 1, alpha = 0.5, , na.rm=TRUE) +
        geom_line(aes(y = predicted_values), colour = "red") +
        geom_vline(xintercept = max_slope_time, linetype="dashed") +
        geom_vline(xintercept = min_slope_time, linetype="dashed") +
        xlab('24hr Time (03:00 to 27:00 for model fitting)') +
        ylab('log(Activity)') +
        theme_classic() +
        labs(title = paste0(subj,' - ',timepoint),
             subtitle = paste0('R2 = ',rsqstr,'; Percent activity missing = ',perc_missstr,'; Aggregated on mean'),
             caption = "Note: Dashed verticle lines depict calculated time of steepest UP and DOWN slope") +
        theme(plot.caption = element_text(hjust = 0))
    ggsave(paste0(subj,'_',timepoint,'_GAMplot.jpg'), plot = p, device = "jpeg", width = 8,
          height = 6, units = 'in')
})