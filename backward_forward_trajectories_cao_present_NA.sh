#!/bin/bash
module load dyn_tools

MODEL=cesm; export MODEL
NETCDF_FORMAT=CF; export NETCDF_FORMAT
LAGRANTO=${DYN_TOOLS}/lagranto.cesm; export LAGRANTO

##  go to run directory
cd /home/nfieldhouse/ma/run/  # change to your run directory

rm folder_traj  
rm folder_trace
rm folder_input

## Initialization parameters
hours_backward_tracing=-168
radius=250	

## set up directories
folder_traj=/net/litho/atmosdyn2/nfieldhouse/data/lagranto/traj/NA  # change path to your data folder (here trajectory files are saved)
# folder_trace=/net/thermo/atmosdyn/fschnyde/data/lagranto/trace # change path to your data folder (here traced trajectory files are saved)
folder_input=/net/litho/atmosdyn2/nfieldhouse/data/lagranto/cao_input/NA # change path to your data folder (text files as input for Lagranto)
folder_dates=/net/litho/atmosdyn2/nfieldhouse/data/cao/events_d  # folder with list of events 
ln -sf ${folder_traj} folder_traj		# Pay attention to the links
ln -sf ${folder_trace} folder_trace
ln -sf ${folder_input} folder_input

echo "folder_input points to: $(readlink -f folder_input)"

# Input file 
file_events=$folder_dates/summary_events_low_CAOI_DJF10_overstep_5_present_NA # list of events
cat $file_events

# Read the number of lines in the file_events
nlines=$(wc -l $file_events | awk '{print $1}')
echo "Number of events: $nlines"

# Loop through all events
for ((i=1; i<=$nlines; i++)); do
    # Read the specific event date
    date_old=$(awk -v line=$i 'NR==line{print $1}' $file_events) # Extract the date from the specified line
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
    echo "Checking startf file: $startf"
    if [ ! -f "$startf" ]; then
        echo "Error: startf file $startf does not exist!"
        continue
    fi
    echo "Contents of startf file:"
    cat $startf | head -n 5  # Print the first 5 lines of the file

    # Create links to files so that they are available for Lagranto
    data_path=/net/meso/climphys/cesm212/b.e212.BHISTcmip6.f09_g17.1000/archive/atm/hist 
    ln -s ${data_path}/b.e212.BHISTcmip6.f09_g17.1000.cam.h3.${year_old}-01-01-00000.nc P${year_old}
    if [ "${month_old}" -eq 01 ]; then
        if [ "${day_old}" -lt 07 ]; then
            ln -s ${data_path}/b.e212.BHISTcmip6.f09_g17.1000.cam.h3.${prevyear}-01-01-00000.nc P${prevyear}
        fi
    fi

    echo "Checking CESM data file: P${year_old}"
    if [ ! -f "P${year_old}" ]; then
        echo "Error: CESM data file P${year_old} does not exist!"
        continue
    fi
    if [ "${month_old}" -eq 01 ] && [ "${day_old}" -lt 07 ]; then
        echo "Checking CESM data file for previous year: P${prevyear}"
        if [ ! -f "P${prevyear}" ]; then
            echo "Error: CESM data file P${prevyear} does not exist!"
            continue
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
done

echo "All events processed."