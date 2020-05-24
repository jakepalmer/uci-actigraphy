##########
# Actigraphy GAM fitting
# Created by: Jake Palmer
# Email: jake.palmer@sydney.edu.au
# Last edit: 09.12.2019
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

setwd('/Users/cheniy1/Desktop/test1/')

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
    df_subj$timehour_shifted <- ifelse(df_subj$timehour < 12.0, df_subj$timehour + 24, df_subj$timehour)
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
    
    #----- NEW
    # time_convert below is a function to convert time from H:M:S to numeric (eg. 23:30 = 23.50).
    # This format is easier to handle in modeling and plotting.You'll see code similar to this used above from row 67 as same process was used for Time values in fitting GAM.
    # I've created a function because it is best practice to move code you use more than once into a function to make it easier to maintain. Here we need it for onset and offset.
    # If you copy and paste some code, it should probably be in a function.
    time_convert <- function(x) {
        if (timeCheck(x)) {
            x <- format(strptime(x, "%I:%M:%S %p"), format = "%H:%M:%S") # If format is 11:30:00 pm, change to 24-hr 23:30:00
        }
        hrs <- as.numeric(format(strptime(x, "%H:%M:%S"), format = "%H")) # select hours
        min <- as.numeric(format(strptime(x, "%H:%M:%S"), format = "%M")) # select mins
        timehour <- round(hrs + (min / 60), 2) # convert to continuous value and round to 2 decimals
        timehour_shifted <- ifelse(timehour < 12.0, timehour + 24, timehour) ###Ivy-Time shifted to 12:00 to 36:00
        timehour <- timehour_shifted
        return(timehour)
    }

    # Read data and filter to include the row for the subject currently being processed (subj defined on row 62)
    sleep_data <- read.csv('SPRiNT_actigraphy SLEEP.csv') %>% filter(ID == subj)
    print("Sleep data in dataframe format, filtered for current subject:")
    print(sleep_data) # can remove this, only here for example
    # Select just the onset/offset time and save it in a variable (sleep_onset) the $ is used to select the column and the [1] selects the row by index.
    # Because we filter to only include the column names and one row for the subj, this will always be 1, referring to the first row of data.
    sleep_onset <- sleep_data$Sleep_Onset_Time[1]
    print(paste0("Sleep onset original format = ", sleep_onset)) # can remove this, only here for example
    sleep_offset <- sleep_data$Sleep_Offset_Time[1]
    # ONSET #
    # Pass sleep onset to time_convert function
    sleep_onset <- time_convert(sleep_onset)
    print(paste0("Sleep onset modified to numeric = ", sleep_onset)) # can remove this, only here for example
    # OFFSET #
    sleep_offset <- time_convert(sleep_offset)
    
    ### Ivy edited for DLMO ###
    ## Read data and filter to include the row for the subject currently being processed (subj defined on row 62)
    DLMO_data <- read.csv('DLMO_all.csv') %>% filter(ID == subj)
    print("DLMO data in dataframe format, filtered for current subject:")
    print(DLMO_data) # can remove this, only here for example
    ## Select just the DLMO time and save it in a variable (DLMO_time) the $ is used to select the column and the [1] selects the row by index.
    ## Because we filter to only include the column names and one row for the subj, this will always be 1, referring to the first row of data.
    DLMO_time <- DLMO_data$DLMO_time[1]
    print(paste0("DLMO time original format = ", DLMO_time)) ## can remove this, only here for example
    ## Pass DLMO time to time_convert function
    DLMO_time <- time_convert(DLMO_time)
    print(paste0("DLMO modified to numeric = ", DLMO_time)) ## can remove this, only here for example
    
    ### Ivy edited for melatonin concentration ###
    ## Read data and filter to include the row for the subject currently being processed (subj defined on row 62)
    mel_data <- read.csv('Melatonin_FINAL.csv') %>% filter(ID == subj)
    print("Melatonin data in dataframe format, filtered for current subject:")
    print(mel_data) # can remove this, only here for example
    ## Select just the collection time and save it in a variable (collection_time) the $ is used to select the column and the [1] selects the row by index.
    ## Because we filter to only include the column names and one row for the subj, this will always be 1, referring to the first row of data.
    collection_time <- mel_data$Time[1:7]
    print(paste0("Collection time original format = ", collection_time)) # can remove this, only here for example
    ## Pass collection time to time_convert function
    collection_time <- time_convert(collection_time)
    print(paste0("Collection time modified to numeric = ", collection_time)) # can remove this, only here for example
    
    #-----
    
    rsqstr <- toString(round(rsq, digits = 2))
    perc_missstr <- toString(round(perc_miss, digits = 2))
    p_check <- appraise(gam_mod)
    ggsave(paste0(subj, '_GAMcheck.jpg'), plot = p_check, device = "jpeg", width = 8,
           height = 6, units = 'in')
    pred_data <- data.frame(time = df_agg$timehour_shifted,
                            activity = df_agg$agglogActivity, 
                            predicted_values = predict(gam_mod, newdata = df_agg))
    
    ###Ivy- Melatonin concentration plot ###
    melatonin_data <- data.frame(time = collection_time,
                                 melatonin = mel_data$Concentration)
    ###Ivy- Summary Table###
    library(gridExtra)
    sum_table <- data.frame("Sleep_onset" = sleep_onset,
                            "Sleep_offset" = sleep_offset - 24,
                            "DOWN_slope" = round(min_slope_time, 2),
                            "UP_slope" = round(max_slope_time, 2) - 24,
                            "DLMO" = DLMO_time)
    sum_table = tableGrob(t(sum_table), theme = ttheme_default(base_size = 8))
    
    ###Ivy- Melatonin curve as main plot, log(Activity) as secondary
        
    p <- ggplot(data = melatonin_data, aes(x = time)) + 
        geom_point(aes(y = melatonin), color = "turquoise4") +
        geom_line(aes(y = melatonin), color = "turquoise4", size = 0.8) +
        geom_segment(aes(x = -Inf, y = 3, xend = 23, yend = 3), linetype = "dotted", color = "turquoise4") + ##Ivy-Melatonin thershold line at 3pg/mL
        annotate(geom = "text", x = 13, y = 4.5, label = "Threshold: 3 pg/mL", color = "turquoise4", size = 3) +
    
        geom_point(data = pred_data, aes(x = time, y = activity*20), size = 1, alpha = 0.5, na.rm=TRUE) +
        geom_line(data = pred_data, aes(x = time, y = predicted_values*20), colour = "red") +
        geom_vline(xintercept = max_slope_time, linetype="dashed") +
        geom_vline(xintercept = min_slope_time, linetype="dashed") +
        scale_y_continuous(sec.axis = sec_axis(~./20, name = "log(Activity)")) +
        
        geom_vline(xintercept = DLMO_time, col = "turquoise4") + ###Ivy-Add a vertical line to indicate DLMO
        annotate(geom = "label", x = DLMO_time, y = 30, label = "DLMO", col = "turquoise4") + ###Ivy-Add a rectangle label for DLMO
        
        #----- NEW
        # Using 3 and 27 as min and max because time axis was shifted by 3 hrs for model fitting
        # Can change colour ("gray" might be better than the green example here) of fill and the opacity (alpha value where 1 = solid colour)
        # Seperate into two annotations if want min and max regions shaded
        #annotate(geom = "rect", xmin = 3.00, xmax = sleep_offset, ymin = -Inf, ymax = Inf,
         #        fill = "slategray1", alpha = 0.5) +
        #annotate(geom = "rect", xmin = sleep_onset, xmax = 27.00, ymin = -Inf, ymax = Inf, 
         #        fill = "slategray1", alpha = 0.5) +
        ## Only need one annotation if want middle part shaded
         annotate(geom = "rect", xmin = sleep_onset, xmax = sleep_offset, ymin = -Inf, ymax = Inf, 
                 fill = "slategray1", alpha = 0.5) +
        
        ###Ivy-Add text to shaded rect area
        #annotate(geom = "text", x = (3.00 + sleep_offset) / 2, y = 2.5, label = "Habitual sleep period") +  
        #annotate(geom = "text", x = (sleep_onset + 27) / 2, y = 2.5, label = "Habitual sleep period") +  
        annotate(geom = "text", x = (sleep_onset + sleep_offset) / 2, y = 50, label = "Habitual sleep period") +
        
        ###Ivy-Add summary table to plot###
        annotation_custom(sum_table, xmin = 12, xmax = 15, ymin = 10, ymax = 25) +
    
        xlab('24hr Time (12:00 to 36:00 for model fitting)') +
        ylab('Salivary Melatonin (pg/mL)') +
        
        ###Ivy- plot melatonin curve###
        #geom_point(data = melatonin_data, aes(x = time, y = melatonin/20), color = "turquoise4") +
        #geom_line(data = melatonin_data, aes(x = time, y = melatonin/20), color = "turquoise4", size = 0.8) +
        #scale_y_continuous(sec.axis = sec_axis(~.*20, name = "Salivary Melatonin (pg/mL)")) + ##Ivy-Secondary y axis
        #geom_segment(aes(x = 15, y = 3/20, xend = Inf, yend = 3/20), linetype = "dashed", color = "turquoise4") + ##Ivy-Melatonin thershold line at 3pg/mL
        
        theme_classic() +
        ###Ivy-change y axis (melatonin) color###
        theme(            
            axis.title.y.left = element_text(color = "turquoise4"),
            axis.text.y.left = element_text(color = "turquoise4"),
            axis.line.y.left = element_line(color = "turquoise4"), 
            axis.ticks.y.left = element_line(color = "turquoise4")
            ) +
        labs(title = paste0(subj),
             subtitle = paste0('R2 = ',rsqstr,'; Percent activity missing = ',perc_missstr,'; Aggregated on mean'),
             caption = "Note: Dashed verticle lines depict calculated time of steepest UP and DOWN slope. ") +
        theme(plot.caption = element_text(hjust = 0))
    ggsave(paste0(subj, '_GAMplot_mela_12to12.jpg'), plot = p, device = "jpeg", width = 8,
           height = 6, units = 'in')
})
print('SCRIPT RUN COMPLETE')