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
folder_traj=/net/litho/atmosdyn2/nfieldhouse/data/lagranto/traj/JS  # change path to your data folder (here trajectory files are saved)
# folder_trace=/net/thermo/atmosdyn/fschnyde/data/lagranto/trace # change path to your data folder (here traced trajectory files are saved)
folder_input=/net/litho/atmosdyn2/nfieldhouse/data/lagranto/cao_input/JS # change path to your data folder (text files as input for Lagranto)
folder_dates=/net/litho/atmosdyn2/nfieldhouse/data/cao/events_d  # folder with list of events 
ln -sf ${folder_traj} folder_traj		# Pay attention to the links
ln -sf ${folder_trace} folder_trace
ln -sf ${folder_input} folder_input

# Input file 
file_events=$folder_dates/dates_CAO_max_intens_1990:2100_JS # list of events
 
nlines=$(wc -l $file_events | awk '{print $1}') # reads how many events you have in the list
echo $nlines

for ((i=1; i<=$nlines; i++)); do   # loop through all events
	# read date from datelist
	date_old=$(awk -v line=$i 'NR==line{print $1}' $file_events) # dates have wrong format --> need to be refomated such that have fromat %Y%m%d_%H

	# date_old="1999-01-10-00"
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

	echo $year_old
	
	echo "Processing event for date: $date"

	# startf as starting positions of trajectories
	startf=folder_input/${year_old}/CAO_${date_old}
	echo "Start file: $startf"

	# Determine the cesm_type based on the year
	if [ "${year_old}" -ge 1990 ] && [ "${year_old}" -le 2000 ]; then
		cesm_type="HIST"
	elif [ "${year_old}" -ge 2090 ] && [ "${year_old}" -le 2100 ]; then
		cesm_type="SSP370"
	else
		echo "Year ${year_old} is not within the supported range (1990-2000 or 2090-2100)."
		exit 1
	fi

	# Set the data path based on the cesm_type
	data_path="/net/meso/climphys/cesm212/b.e212.B${cesm_type}cmip6.f09_g17.1000/archive/atm/hist"

	# Create links to files so that they are available for Lagranto
	ln -s ${data_path}/b.e212.B${cesm_type}cmip6.f09_g17.1000.cam.h3.${year_old}-01-01-00000.nc P${year_old}
	if [ "${month_old}" -eq 01 ]; then
		if [ "${day_old}" -lt 07 ]; then
			ln -s ${data_path}/b.e212.B${cesm_type}cmip6.f09_g17.1000.cam.h3.${prevyear}-01-01-00000.nc P${prevyear}
		fi
	fi

	# outfile
	outfile_back=B_${date}
	# outfile_trace_back=BT_${date}

	# set up tracevars file for tracing trajectories
	echo "Q 1000. P" >& tracevars
	echo "T 1. P" >> tracevars
	echo "PTTEND 3600. P" >> tracevars
	echo "SHFLX 1 P" >> tracevars
	echo "LHFLX 1 P" >> tracevars
	echo "SST 1. P" >> tracevars
	echo "PSL 0.01 P" >> tracevars

	# calculate 7day backward trajectories from all positions stored in startf. (caltra startfile time outfile -o 360 -i 360 -ref reference_date -rs 30 -cesm2)
	caltra ${startf} -168 ${outfile_back} -o 360 -i 360 -ref ${date} -ts 30 -cesm2

	# move traje and traj file to right folder
	mv $outfile_back folder_traj/$outfile_back
			
	# remove links
	rm P${year_old}
	if [ "${month_old}" -eq 01 ]; then
		if [ "${day_old}" -lt 07 ]; then
			rm  P${prevyear}
		fi
	fi

	\rm -f caltra.param trace.param startf


done 	# end loop


