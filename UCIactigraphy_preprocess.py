##########
# Actigraphy Pre-Processing
# Created by: Jake Palmer
# Email: jake.palmer@sydney.edu.au
# Last edit: 12.12.2019
##########

#####
# INSTRUCTIONS
#####
# This notebook has two sections:
# 1. Code to extract the standard actigraphy measures (such as Sleep Onset and Offset Time).
# 2. Code to trim the files to only include only epoch-by-epoch activity data (complete days) to be used for further processing.

# Before running...
# 1. Check the variable 'base' has the correct path to the data that you intend to process
# 2. Choose which sections of the script to run (i.e. extract standard variables/trim file or both). 
#    To do this, just fill the 'extract' and 'trim' variables with 'yes' or 'no' (make sure to keep the quotes)

#####
# EDIT THESE VARIABLES
#####
base = '/Users/mq44848301/Desktop/test'
extract = 'yes'
trim = 'yes'
compile_trimmed_files = 'yes'
# If compiling trimmed files, assign filename for SPSS output file (must have .sav extension)
compiled_file = 'compiled_output_test.sav'

#####
# SECTION 1 - Extract standard variables
#####

import glob
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pyreadstat
pd.options.mode.chained_assignment = None

base_dir = Path(base)
os.chdir(base_dir)

# FUNCTIONS


def get_subj_ID(f):
    """Get subject related variables for file naming"""
    subj = f.split(".")
    subj_ID = subj[0]
    subj_trim_fname = subj_ID + '.csv'
    subj_orig_fname = subj_ID + '_orig' + '.csv'
    return subj_ID, subj_trim_fname, subj_orig_fname


def read_file(f):
    """Read in csv files"""
    # Am providing dummy column names as irregular files can occasionally be problematic
    cols = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
    df = pd.read_csv(f, header=None, low_memory=False, names=cols)
    df = df.dropna(axis=0, how='all')
    df.reset_index(drop=True, inplace=True)
    return df


def mean_time(df, col):
    df[col] = pd.to_datetime(df[col])
    df[col] = df[col].dt.strftime('%H:%M:%S')
    times = df[col].tolist()
    times = [x for x in times if x != 'NaT']
    t = (time.split(':') for time in times)
    seconds = ((int(s) + int(m) * 60 + int(h) * 3600)
               for h, m, s in t)
    to_seconds = [s for s in seconds]
    to_seconds = [s + 86400 if s <= 43200 else s for s in to_seconds]
    mean_seconds = sum(to_seconds) / len(to_seconds)
    mean_h, mean_m = divmod(mean_seconds, 3600)
    mean_m, mean_s = divmod(mean_m, 60)
    if mean_h >= 24:
        mean_h = mean_h - 24
    sd_seconds = np.std(to_seconds)
    sd_h, sd_m = divmod(sd_seconds, 3600)
    sd_m, sd_s = divmod(sd_m, 60)
    mean = '%02i:%02i:%02i' % (mean_h, mean_m, mean_s)
    sd = '%02i.%02i' % (sd_m, sd_s)
    return mean, sd


def get_value(col, df):
    if col in df.columns:
        mean = df.loc[df.index[3], col]
        var = df.loc[df.index[4], col]
    else:
        mean = 'NA'
        var = 'NA'
    return mean, var


# MAIN ROUTINE


def extract_standard_values(file, interval_type, interval_type_sum, out_df, period):
    """Function that calls main routine to extract standard values.
    Moved to a function so can be called for each of REST, ACTIVITY and SLEEP."""
    subj_ID, subj_trim_fname, subj_orig_fname = get_subj_ID(file)
    df = read_file(file)
    # Get summary data
    mark_start = ['------------------------ Statistics ------------------------']
    mark_end = ['--------------------- Marker/Score List --------------------']
    mark_start_index = df[df['a'].isin(mark_start)].index[0]
    mark_start_index = mark_start_index + 1
    mark_end_index = df[df['a'].isin(mark_end)].index[0]
    df = df.iloc[mark_start_index:mark_end_index]
    new_header = df.iloc[0]  # Get first row for the header (instead of integers as header)
    new_header = [x.strip().replace(' ', '_') for x in new_header]  # Remove white space from new header
    new_header = [x.strip().replace('-', '_') for x in new_header]
    df = df[1:]  # Get df minus the header row
    df.columns = new_header  # Set the header row as the df header
    df = df.dropna(subset=['Interval_Type'])
    # Get data for each interval type
    df_interval = df[df.Interval_Type == interval_type]  # Select the rows that are for given interval type only
    df_interval_sum = df[df.Interval_Type == interval_type_sum]  # Select the rows that are for rest_summary only
    # Get specific summary values
    if period in ('REST', 'ACTIVE'):
        onset_time, onset_time_var = mean_time(df_interval, 'Start_Time')  # Get M and SD for onset times
        offset_time, offset_time_var = mean_time(df_interval, 'End_Time')  # Get M and SD for offset times
        TiB, TiB_var = get_value('Duration', df_interval_sum)
        TST, TST_var = get_value('Sleep_Time', df_interval_sum)
        Time_Awake, Time_Awake_var = get_value('Wake_Time', df_interval_sum)
        perc_Time_Awake, perc_Time_Awake_var = get_value('%Wake', df_interval_sum)
        Wake_bouts, Wake_bouts_var = get_value('#Wake_Bouts', df_interval_sum)
        total_AC, total_AC_var = get_value('Total_AC', df_interval_sum)
        # Put summary values in df for output
        out_df = out_df.append({'ID': subj_ID, 'Onset_Time': onset_time, 'Onset_Time_var': onset_time_var,
                                'Offset_Time': offset_time, 'Offset_Time_var': offset_time_var, 'TiB': TiB,
                                'TiB_var': TiB_var, 'TST': TST, 'TST_var': TST_var, 'Time_Awake': Time_Awake,
                                'Time_Awake_var': Time_Awake_var, 'Time_Awake_perc': perc_Time_Awake,
                                'Time_Awake_perc_var': perc_Time_Awake_var, 'Wake_bouts': Wake_bouts,
                                'Wake_bouts_var': Wake_bouts_var, 'Total_AC': total_AC}, ignore_index=True)
    elif period == 'SLEEP':
        onset_time, onset_time_var = mean_time(df_interval, 'Start_Time')  # Get M and SD for onset times
        offset_time, offset_time_var = mean_time(df_interval, 'End_Time')  # Get M and SD for offset times
        TiB, TiB_var = get_value('Duration', df_interval_sum)
        OffWrist_mins, OffWrist_mins_var = get_value('Off_Wrist', df_interval_sum)
        OffWrist_perc, OffWrist_perc_var = get_value('%Off_Wrist', df_interval_sum)
        SOL, SOL_var = get_value('Onset_Latency', df_interval_sum)
        SE, SE_var = get_value('Efficiency', df_interval_sum)
        WASO, WASO_var = get_value('WASO', df_interval_sum)
        Time_Awake, Time_Awake_var = get_value('Wake_Time', df_interval_sum)
        perc_Time_Awake, perc_Time_Awake_var = get_value('%Wake', df_interval_sum)
        TST, TST_var = get_value('Sleep_Time', df_interval_sum)
        sleep_perc, sleep_perc_var = get_value('%Sleep', df_interval_sum)
        # Put summary values in df for output
        out_df = out_df.append({'ID': subj_ID, 'Onset_Time': onset_time, 'Onset_Time_var': onset_time_var,
                                'Offset_Time': offset_time, 'Offset_Time_var': offset_time_var, 'TiB': TiB,
                                'TiB_var': TiB_var, 'OffWrist_mins': OffWrist_mins, 'OffWrist_mins_var': OffWrist_mins_var,
                                'OffWrist_percent': OffWrist_perc, 'OffWrist_percecnt_var': OffWrist_perc_var,
                                'SOL': SOL, 'SOL_var': SOL_var, 'SE': SE, 'SE_var': SE_var, 'WASO': WASO, 'WASO_var': WASO_var,
                                'Time_Awake': Time_Awake, 'Time_Awake_var': Time_Awake_var, 'Time_Awake_percent': perc_Time_Awake,
                                'Time_Awake_percent_var': perc_Time_Awake_var, 'TST': TST, 'TST_var': TST_var,
                                'Sleep_percent': sleep_perc, 'Sleep_percent_var': sleep_perc_var}, ignore_index=True)
    return out_df


# Initiate empty data frames to append values to
REST_out_df = pd.DataFrame(
    columns=['ID', 'Onset_Time', 'Onset_Time_var', 'Offset_Time', 'Offset_Time_var', 'TiB', 'TiB_var', 'TST', 'TST_var', 'Time_Awake',
             'Time_Awake_var', 'Time_Awake_perc', 'Time_Awake_perc_var', 'Wake_bouts', 'Wake_bouts_var', 'Total_AC'])
ACTIVE_out_df = pd.DataFrame(
    columns=['ID', 'Onset_Time', 'Onset_Time_var', 'Offset_Time', 'Offset_Time_var', 'TiB', 'TiB_var', 'TST', 'TST_var', 'Time_Awake',
             'Time_Awake_var', 'Time_Awake_perc', 'Time_Awake_perc_var', 'Wake_bouts', 'Wake_bouts_var', 'Total_AC'])
SLEEP_out_df = pd.DataFrame(
    columns=['ID', 'Onset_Time', 'Onset_Time_var', 'Offset_Time', 'Offset_Time_var', 'TiB', 'TiB_var', 'OffWrist_mins', 'OffWrist_mins_var',
             'OffWrist_percent', 'OffWrist_percecnt_var', 'SOL', 'SOL_var', 'SE', 'SE_var', 'WASO', 'WASO_var', 'Time_Awake', 'Time_Awake_var',
             'Time_Awake_percent', 'Time_Awake_percent_var', 'TST', 'TST_var', 'Sleep_percent', 'Sleep_percent_var'])

### Add if then here for runnning extract
if extract == 'yes':
    print('Extracting standard variables:')
    for file in glob.glob('*.csv'):
        if "standard_variables" in file:
            continue
        elif "trimmed" in file:
            continue
        else:
            print('    Working on', file,'...')
            # Extract standard variables
            REST_out_df = extract_standard_values(file, 'REST', 'Rest Summary', REST_out_df, 'REST')
            ACTIVE_out_df = extract_standard_values(file, 'ACTIVE', 'Active Summary', ACTIVE_out_df, 'ACTIVE')
            SLEEP_out_df = extract_standard_values(file, 'SLEEP*', 'Sleep Summary', SLEEP_out_df, 'SLEEP')
            # Write output data frames to file
            rest_file = base_dir / 'REST_standard_variables.csv'
            active_file = base_dir / 'ACTIVE_standard_variables.csv'
            sleep_file = base_dir / 'SLEEP_standard_variables.csv'
            REST_out_df.to_csv(rest_file, header=True, index=False)
            ACTIVE_out_df.to_csv(active_file, header=True, index=False)
            SLEEP_out_df.to_csv(sleep_file, header=True, index=False)
    print("Done extracting standard variables!")
else:
    print("Skipping standard variable extraction...")


#####
# SECTION 2 - Trim files
#####

# FUNCTIONS


def get_subj_ID(f):
    """Get subject related variables for file naming"""
    subj = f.split(".")
    subj_ID = subj[0]
    subj_trim_fname = subj_ID + '_trimmed.csv'
    return subj_ID, subj_trim_fname


def read_file(f):
    """Read in csv files"""
    cols = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
    df = pd.read_csv(f, header=None, low_memory=False, names=cols)
    df = df.dropna(axis=0, how='all')  # Drop blank rows
    df.reset_index(drop=True, inplace=True)
    return df


# MAIN ROUTINE
if trim == 'yes':
    print('Trimming files:')
    for file in glob.glob('*.csv'):
        if "standard_variables" in file:
            continue
        elif "trimmed" in file:
            continue
        else:
            print('    Trimming file', file, '...')
            subj_ID, subj_trim_fname = get_subj_ID(file)
            df = read_file(file)
            mark_start = ['S/W Status']
            mark_start_index = df[df['m'].isin(mark_start)].index[0]
            df = df.iloc[mark_start_index:]
            # Get first row for the header (instead of integers as header)
            new_header = df.iloc[0]
            df = df[1:]  # Get df minus the header row
            df.columns = new_header  # Set the header row as the df header
            df = df.dropna(axis=1, how='all')
            df.to_csv(subj_trim_fname, index=False, na_rep='NaN')
    print("Done trimming files!")
else:
    print("Skipping file trimming...")

#####
# SECTION 3 - Compile trimmed files
#####

# MAIN ROUTINE
if compile_trimmed_files == 'yes':
    print('Compiling files:')
    dfs = []
    files = [f for f in glob.glob("*_trimmed.csv")]
    files.sort()
    if len(files) == 0:
        print('ERROR: No trimmed files to compile...')
    else:
        for file in files:
            print('    Adding file', file, '...')
            ID = file.rstrip('_trimmed.csv')
            df_name = ID + '_df'
            df_name = pd.read_csv(file)
            df_name['ID'] = ID
            dfs.append(df_name)
        compiled_df = pd.concat(dfs, axis=0, ignore_index=True, sort=False)
        cols = compiled_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        compiled_df = compiled_df[cols]
        compiled_df.columns = compiled_df.columns.str.replace(' ', '_')
        compiled_df.columns = compiled_df.columns.str.replace('-', '_')
        compiled_df.columns = compiled_df.columns.str.replace('/', '_')
        print('Writing SPSS output...')
        compiled_df = compiled_df.astype(str)
        pyreadstat.write_sav(compiled_df, compiled_file)

print('SCRIPT RUN COMPLETE')
