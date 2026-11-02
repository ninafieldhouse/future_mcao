#!/bin/bash
source "$(dirname "$0")/../config/config.sh"

module load dyn_tools

MODEL=cesm; export MODEL
NETCDF_FORMAT=CF; export NETCDF_FORMAT
LAGRANTO=${DYN_TOOLS}/lagranto.cesm; export LAGRANTO

member_id=$1
timeframe=$2  # present or future
region=$3
date_old=$4


if [ "$timeframe" == "present" ]; then
    cesm_type="HIST"
elif [ "$timeframe" == "future" ]; then
    cesm_type="SSP370"
else
    echo "Invalid timeframe. Use 'present' or 'future'."
    exit 1
fi

##  go to run directory
cd "$CAO_RUN_DIR/"  # change to your run directory

rm folder_traj  
rm folder_trace
rm folder_input

## Initialization parameters
hours_backward_tracing=-168
radius=250	
F
## set up directories
folder_traj=${CAO_DATA_ROOT}/${member_id}/lagranto/traj/  # change path to your data folder (here trajectory files are saved)
# folder_trace=/net/thermo/atmosdyn/fschnyde/data/lagranto/trace # change path to your data folder (here traced trajectory files are saved)
folder_input=${CAO_DATA_ROOT}/${member_id}/lagranto/cao_input/${region} # change path to your data folder (text files as input for Lagranto)
echo $folder_input
rm -rf folder_traj folder_input  # Remove existing links to avoid conflicts

mkdir -p ${folder_traj}  # Create the trajectory folder if it doesn't exist
mkdir -p ${folder_input}  # Create the input folder if it doesn't exist

ln -sf ${folder_traj} folder_traj		# Pay attention to the links
# ln -sf ${folder_trace} folder_trace
ln -sf ${folder_input} folder_input

# --------------------
## For one event
# --------------------
# Specify the event index (line number) to process
event_index=$4  # Change this to the desired event line number (e.g., 1 for the first event)

# Read the specific event date
# date_old=$(awk -v line=$event_index 'NR==line{print $1}' $file_events) # Extract the date from the specified line
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

# Create links to files so that they are available for Lagranto
data_path=${CAO_CESM_ARCHIVE}/b.e212.B${cesm_type}cmip6.f09_g17.${member_id}/archive/atm/hist
ln -s ${data_path}/b.e212.B${cesm_type}cmip6.f09_g17.${member_id}.cam.h3.${year_old}-01-01-00000.nc P${year_old}
if [ "${month_old}" -eq 01 ]; then
    if [ "${day_old}" -lt 07 ]; then
        ln -s ${data_path}/b.e212.B${cesm_type}cmip6.f09_g17.${member_id}.cam.h3.${prevyear}-01-01-00000.nc P${prevyear}
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

echo "startf: ${folder_input}/${year_old}/CAO_${date_old}"
echo "outfile_back: ${outfile_back}"
echo "date: ${date}"

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

# --------------------
## For all events
# --------------------
 
# nlines=$(wc -l $file_events | awk '{print $1}') # reads how many events you have in the list
# echo $nlines

# for ((i=1; i<=$nlines; i++)); do   # loop through all events
# 	# read date from datelist
# 	date_old=$(awk -v line=$i 'NR==line{print $1}' $file_events) # dates have wrong format --> need to be refomated such that have fromat %Y%m%d_%H

# 	# date_old="1999-01-10-00"
# 	year_old=$(echo -e ${date_old} | cut -c1-4)
# 	month_old=$(echo -e ${date_old} | cut -c6-7)
# 	day_old=$(echo -e ${date_old} | cut -c9-10)
# 	hour_old=$(echo -e ${date_old} | cut -c12-13)

# 	date="${year_old}${month_old}${day_old}_${hour_old}"

# 	if [ "${month_old}" -eq 01 ]; then
# 		if [ "${day_old}" -lt 07 ]; then
# 			prevyear=$(echo -e $(newtime ${date} -8760) | cut -c1-4) 
# 		fi
# 	fi	 

# 	echo $year_old 

# 	# startf as starting positions of trajectories
# 	startf=folder_input/${year_old}/CAO_${date_old}
# 	echo $startf

# 	# Create links to files so that they are available for Lagranto
# 	data_path=/net/meso/climphys/cesm212/b.e212.BHISTcmip6.f09_g17.1000/archive/atm/hist 
# 	ln -s ${data_path}/b.e212.BHISTcmip6.f09_g17.1000.cam.h3.${year_old}-01-01-00000.nc P${year_old}
# 	if [ "${month_old}" -eq 01 ]; then
# 		if [ "${day_old}" -lt 07 ]; then
# 			ln -s ${data_path}/b.e212.BHISTcmip6.f09_g17.1000.cam.h3.${prevyear}-01-01-00000.nc P${prevyear}
# 		fi
# 	fi

# 	# outfile
# 	outfile_back=B_${date}
# 	# outfile_trace_back=BT_${date}

# 	# set up tracevars file for tracing trajectories
# 	echo "Q 1000. P" >& tracevars
# 	echo "T 1. P" >> tracevars
# 	echo "PTTEND 1. P" >> tracevars
# 	echo "SHFLX 1 P" >> tracevars
# 	echo "LHFLX 1 P" >> tracevars
# 	echo "SST 1. P" >> tracevars
# 	echo "PSL 0.01 P" >> tracevars

# 	# calculate 7day backward trajectories from all positions stored in startf. (caltra startfile time outfile -o 360 -i 360 -ref reference_date -rs 30 -cesm2)
# 	caltra ${startf} -168 ${outfile_back} -o 360 -i 360 -ref ${date} -ts 30 -cesm2

# 	# move traje and traj file to right folder
# 	mv $outfile_back folder_traj/$outfile_back
			
# 	# remove links
# 	rm P${year_old}
# 	if [ "${month_old}" -eq 01 ]; then
# 		if [ "${day_old}" -lt 07 ]; then
# 			rm  P${prevyear}
# 		fi
# 	fi

# 	\rm -f caltra.param trace.param startf


# done 	# end loop


