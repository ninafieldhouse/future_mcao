#!/bin/bash
source "$(dirname "$0")/../config/config.sh"

module load dyn_tools

MODEL=cesm; export MODEL
NETCDF_FORMAT=CF; export NETCDF_FORMAT
LAGRANTO=${DYN_TOOLS}/lagranto.cesm; export LAGRANTO

## Check if event_index is provided
if [ -z "$1" ]; then
    echo "Error: No event_index provided. Usage: $0 <event_index>"
    exit 1
fi

event_index=$1  # Get the event index from the first argument
member_id=$3  # Get the member ID from the third argument

## Go to run directory
run_dir="${CAO_RUN_DIR}/NA_present"  # Change to your run directory
mkdir -p $run_dir  # Create the directory if it doesn't exist
cd $run_dir  # Change to your run directory

rm folder_traj  
rm folder_trace
rm folder_input

## Initialization parameters
hours_backward_tracing=-168
radius=250	
## Set up directories
folder_traj=${CAO_DATA_ROOT}/${member_id}/lagranto/traj/NA  # Change path to your data folder (here trajectory files are saved)
folder_input=${CAO_DATA_ROOT}/${member_id}/lagranto/cao_input/NA # Change path to your data folder (text files as input for Lagranto)
rm -rf folder_traj folder_input  # Remove existing links to avoid conflicts

mkdir -p ${folder_traj}  # Create the trajectory folder if it doesn't exist
mkdir -p ${folder_input}  # Create the input folder if it doesn't exist
ln -sf ${folder_traj} folder_traj
ln -sf ${folder_input} folder_input

# Input file 
file_events=$2 # List of events, passed as the second argument

# Read the specific event date
date_old=$(awk -v line=$event_index 'NR==line{print $1}' $file_events) # Extract the date from the specified line
year_old=$(echo -e ${date_old} | cut -c1-4)
month_old=$(echo -e ${date_old} | cut -c6-7)
day_old=$(echo -e ${date_old} | cut -c9-10)
hour_old=$(echo -e ${date_old} | cut -c12-13)

date="${year_old}${month_old}${day_old}_${hour_old}"

if [ "${month_old}" -eq 01 ]; then
    if [ "${day_old}" -lt 07 ]; then
        prevyear=$(echo -e $(newtime ${date} -8760) | cut -c1-4) 
    fi
fi	 

echo "Processing event for date: $date"

# startf as starting positions of trajectories
startf=folder_input/${year_old}/CAO_${date_old}
echo "Start file: $startf"
if [ ! -f "$startf" ]; then
    echo "Error: startf file $startf does not exist! Exiting."
    exit 1
fi

# Create links to files so that they are available for Lagranto
data_path=${CAO_CESM_ARCHIVE}/b.e212.BHISTcmip6.f09_g17.${member_id}/archive/atm/hist
ln -s ${data_path}/b.e212.BHISTcmip6.f09_g17.${member_id}.cam.h3.${year_old}-01-01-00000.nc P${year_old}
if [ "${month_old}" -eq 01 ]; then
    if [ "${day_old}" -lt 07 ]; then
        ln -s ${data_path}/b.e212.BHISTcmip6.f09_g17.${member_id}.cam.h3.${prevyear}-01-01-00000.nc P${prevyear}
    fi
fi

# outfile
outfile_back=B_${date}

# set up tracevars file for tracing trajectories
echo "Q 1000. P" >& tracevars
echo "T 1. P" >> tracevars
echo "PTTEND 3600. P" >> tracevars #scale here (check from ncview) [K/6h] --> scale with 10^4
echo "SHFLX 1 P" >> tracevars
echo "LHFLX 1 P" >> tracevars
echo "SST 1. P" >> tracevars
echo "PSL 0.01 P" >> tracevars

# calculate 7-day backward trajectories from all positions stored in startf
caltra ${startf} -168 ${outfile_back} -o 360 -i 360 -ref ${date} -ts 30 -cesm2

# move traje and traj file to the right folder
mv $outfile_back folder_traj/$outfile_back
            
# remove links
rm P${year_old}
if [ "${month_old}" -eq 01 ]; then
    if [ "${day_old}" -lt 07 ]; then
        rm  P${prevyear}
    fi
fi

\rm -f caltra.param trace.param startf

echo "Processing complete for event: $date"