# 1) Do a ranked time list of all the time steps in a given SEASON and remove all the time steps below a certain threshold (90th)
# 2) Then the events must be defined iteratively, beginning with the most intense heat flux. If consecutive time steps are in the list, add them to this event. 
# 3) If the time step is not found, keep looking for 12/24 hours and see if it resumes. If yes, attribute it to that event.
# 4) The first and last time step of each event are the beginning and the end. Now eliminate from the list the time steps from the beginning to the end.
# 5) Start again until there are no more time steps in the 90th percentile. 

###### NEW VERSION: Now in the file sumflux there are all the three variables. To avoid changing too much the code, create help files with simply date and value.
cd /home/nfieldhouse/ma/run
percentile=90			#Select 90 or 95
overstep_hours_max=24		#Free selection (must be positive, preferably 12/24/48h)
variable="CAOI"			#Select sshf or slhf
season="DJF"			#DJF, MAM, JJA, SON or ALL for all years
timestep=6

member_id=$1  # Get from command line
outdir="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/events_d"
mkdir -p ${outdir}

# Define regions and timeframes
regions=("JS" "NA")
timeframes=("present" "future") # "present" for 1990-2000, "future" for 2090-2100

# Loop through regions and timeframes
for region in "${regions[@]}"
do
	for timeframe in "${timeframes[@]}"
	do
		base_path="/net/litho/atmosdyn2/nfieldhouse/data/${member_id}/cao/clim/${region}"
		flux_list_file_original="${base_path}/CAOI_${timeframe}.csv"
		flux_list_file_timeframe="${outdir}/clim_${variable}_${timeframe}"

		# Check input file
		if [[ ! -s "${flux_list_file_original}" ]]; then
			echo "WARNING: Input file ${flux_list_file_original} does not exist or is empty. Skipping."
			continue
		fi

		awk '{print $1,$2*(-1.0)}' ${flux_list_file_original} > ${flux_list_file_timeframe} #multiply by -1.0 to make highest 10th percentile correspond to lowest 10th percentile
		if [[ ! -s "${flux_list_file_timeframe}" ]]; then
			echo "ERROR: Output file ${flux_list_file_timeframe} was not created or is empty after awk."
			continue
		fi

		echo "Processing region=${region}, timeframe=${timeframe}, variable=${variable}"
		echo "Input: ${flux_list_file_original}"
		echo "Output: ${flux_list_file_timeframe}"

		#Let's start with 1979_2022
		flux_list_file=${flux_list_file_timeframe}
		flux_list_sorted=${flux_list_file}_sorted

		events_list_file="${outdir}/events_${season}${percentile}_${timeframe}_overstep_${overstep_hours_max}_${variable}_${region}"
		events_list_file_summary="${outdir}/summary_events_${season}${percentile}_${timeframe}_overstep_${overstep_hours_max}_${variable}_${region}"

		#events_list_file="TEST"
		rm -f ${events_list_file} ${events_list_file_summary}
		touch -f ${events_list_file} ${events_list_file_summary}

		#Creation of file must be handled with care. Creation of three different files.


		#Sort and create backup of data

		LC_ALL=C sort -g -k 2,2 ${flux_list_file} > ${flux_list_sorted}
		cat ${flux_list_sorted} > ${flux_list_sorted}_backup


		#Retain only sshf/slhf values for each season if you select a certain season (starting from winter 1979/1980, excluding Jan and Feb 2000). 

		awk '{print $1,$2*(1.0)}' ${flux_list_sorted} > shorter_file_${season}

		if [[ "${season}" == "DJF" ]] ; then
			if [[ "${timeframe}" == "present" ]]; then
				# Handle present timeframe (e.g., 1990–2000)
				awk 'substr($1,6,2)=="12" && substr($1,1,4)<2000 || substr($1,6,2)=="01" && substr($1,1,4)>1990 || substr($1,6,2)=="02" && substr($1,1,4)>1990 {print $1,$2*(1.0)}' ${flux_list_sorted} > shorter_file_${season}
			elif [[ "${timeframe}" == "future" ]]; then
				# Handle future timeframe (e.g., 2090–2100)
				awk 'substr($1,6,2)=="12" && substr($1,1,4)<2100 || substr($1,6,2)=="01" && substr($1,1,4)>2090 || substr($1,6,2)=="02" && substr($1,1,4)>2090 {print $1,$2*(1.0)}' ${flux_list_sorted} > shorter_file_${season}
			fi
		fi

		mv shorter_file_${season} ${flux_list_sorted}

		#Extract upper 90 percentile and 95 percentile

		perc90=$(awk '{all[NR] = $2} END{print all[int(NR*0.9 - 0.5)]}' ${flux_list_sorted})
		echo -e ${perc90}
		perc10=$(awk '{all[NR] = $2} END{print all[int(NR*0.1 + 0.5)]}' ${flux_list_sorted})
		echo -e ${perc10}

		perc95=$(awk '{all[NR] = $2} END{print all[int(NR*0.95 - 0.5)]}' ${flux_list_sorted})
		echo -e ${perc95}
		perc05=$(awk '{all[NR] = $2} END{print all[int(NR*0.05 + 0.5)]}' ${flux_list_sorted})
		echo -e ${perc05}

		if (("$percentile" == "90")); then

			perc_up=${perc90}
			perc_down=${perc10}
		fi

		if (("$percentile" == "95")); then

			perc_up=${perc95}
			perc_down=${perc05}
		fi

		#Re-sort after extraction of data

		LC_ALL=C sort -g -k 2,2 ${flux_list_sorted} > ${flux_list_sorted}_${season}
		mv ${flux_list_sorted}_${season} ${flux_list_sorted}

		#Now remove the ones below the 90th percentile

		awk -v perc=${perc_up} '$2>perc{print $0}' ${flux_list_sorted} > ${flux_list_sorted}_${percentile}
		cp ${flux_list_sorted}_${percentile} ${flux_list_sorted}_timesteplist_${region}
		mv ${flux_list_sorted}_${percentile} ${flux_list_sorted}
		# Check sorted file
		if [[ ! -s "${flux_list_sorted}" ]]; then
			echo "ERROR: Sorted file ${flux_list_sorted} was not created or is empty."
			continue
		fi
		#mv ${flux_list_sorted}_${percentile} TEST

		#Now I have to write a loop to do it

		fluxmax_val=$(tail -1 ${flux_list_sorted} | awk '{print $2}')
		nevent=0

		#First loop, to consider all the values down to the 90th percentile
		while [ "$(bc <<< "${fluxmax_val} > ${perc_up}")" == "1" ]; do 

			#Extract line at the end of the file
			year_max=$(tail -1 ${flux_list_sorted} | awk '{print substr($1,1,4)}')
			month_max=$(tail -1 ${flux_list_sorted} | awk '{print substr($1,6,2)}')
			day_max=$(tail -1 ${flux_list_sorted} | awk '{print substr($1,9,2)}')
			time_max=$(tail -1 ${flux_list_sorted} | awk '{print substr($1,12,2)}')
			fluxmax_val=$(tail -1 ${flux_list_sorted} | awk '{print $2}')
			d_max="${year_max}-${month_max}-${day_max} ${time_max}:00:00"				# This remains fixed
			
			nevent=$((${nevent}+1))
			echo -e ""
			echo -e "Event" ${nevent} ${d_max}
			echo -e ""

			\rm -f overstep_file		#If previously existent after abrupt interruption
			event_file="event_${nevent}"
			touch -f ${event_file}
			event_file_forward="event_forward_${nevent}"
			touch -f ${event_file_forward}
			event_file_backward="event_backward_${nevent}"
			touch -f ${event_file_backward}

			#Attribute line to the event and write it in ${event_file_forward}
			line_summary=$(awk -v year=${year_max} -v month=${month_max} -v day=${day_max} -v hour=${time_max} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_file})		#I extract if from this file to include overstepping time steps
			echo -e ${line_summary} ${nevent} >> ${event_file_forward}

			#Eliminate line with maximum
			awk -v year=${year_max} -v month=${month_max} -v day=${day_max} -v hour=${time_max} 'substr($1,1,4)!=year || substr($1,6,2)!=month || substr($1,9,2)!=day || substr($1,12,2)!=hour' ${flux_list_sorted} > shorter_file
			mv shorter_file ${flux_list_sorted}

			#Begin FORWARD loop to see if also other dates belong to the event
			overstep_counter=0
			d="${year_max}-${month_max}-${day_max} ${time_max}:00:00"				# This is modified
			while [ "$(bc <<< "${overstep_counter} <= ${overstep_hours_max}")" == "1" ]; do 

				touch -f overstep_file		#In case overstepping is required

				d=$(date -d "${d} ${timestep}hour" +"%Y-%m-%d %T")
				echo -e "FORWARD" ${d}

					year_n=$(echo -e ${d} | cut -c1-4)
					month_n=$(echo -e ${d} | cut -c6-7)
					day_n=$(echo -e ${d} | cut -c9-10)
					time_n=$(echo -e ${d} | cut -c12-13)

				#Check for overstep
				awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_sorted} > check_file
				check_flag=$(wc -l check_file  | awk '{print $1}')	#It should be one if the date is present, 0 otherwise

				echo -e ${check_flag} ${overstep_counter}

				#If time step is missing, increase by one the overstep_counter and add the line to overstep_file
				if [ "${check_flag}" -eq "0" ]; then
					overstep_counter=$((${overstep_counter}+${timestep}))
					line=$(awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_file})			#I extract if from this file to include overstepping time steps
					touch -f overstep_file
					echo -e ${line} ${nevent} >> overstep_file
				fi
				#If time step is present, no overstep is required. If previously overstepped time steps were recorded, then add them to event_file_forward
				if [ "${check_flag}" -eq "1" ]; then
					touch -f overstep_file
					cat overstep_file >> ${event_file_forward}
					\rm -f overstep_file
					#Attribute line to the event and write it in ${event_file_forward}
					line=$(awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_file})			#I extract if from this file to include overstepping time steps
					echo -e ${line} ${nevent} >> ${event_file_forward}
					overstep_counter=0
				fi

				#Eliminate line from file (if present)

				awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)!=year || substr($1,6,2)!=month || substr($1,9,2)!=day || substr($1,12,2)!=hour' ${flux_list_sorted} > shorter_file
				mv shorter_file ${flux_list_sorted}
			done	#End FORWARD loop

			\rm -f overstep_file	#This file must be removed, so that it does not appear on the backward file

			#Begin BACKWARD loop to see if also other dates belong to the event
			overstep_counter=0
			d="${year_max}-${month_max}-${day_max} ${time_max}:00:00"	#Re-initialize date
			while [ "$(bc <<< "${overstep_counter} <= ${overstep_hours_max}")" == "1" ]; do 

				d=$(date -d "${d} ${timestep}hour ago" +"%Y-%m-%d %T")
				echo -e "BACKWARD" ${d}

					year_n=$(echo -e ${d} | cut -c1-4)
					month_n=$(echo -e ${d} | cut -c6-7)
					day_n=$(echo -e ${d} | cut -c9-10)
					time_n=$(echo -e ${d} | cut -c12-13)

				#Check for overstep
				awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_sorted} > check_file
				check_flag=$(wc -l check_file  | awk '{print $1}')	#It should be one if the date is present, 0 otherwise

				echo -e ${check_flag} ${overstep_counter}

				#If time step is missing, increase by one the overstep_counter. Write line in overstep_file
				if [ "${check_flag}" -eq "0" ]; then
					overstep_counter=$((${overstep_counter}+${timestep}))
					line=$(awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_file})			#I extract if from this file to include overstepping time steps
					touch -f overstep_file
					echo -e ${line} ${nevent} >> overstep_file
				fi
				#If time step is present, no overstep is required
				if [ "${check_flag}" -eq "1" ]; then
					touch -f overstep_file
					cat overstep_file >> ${event_file_backward}
					\rm -f overstep_file
					#Attribute line to the event and write it in ${event_file_backward}
					line=$(awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)==year && substr($1,6,2)==month && substr($1,9,2)==day && substr($1,12,2)==hour' ${flux_list_file})			#I extract if from this file to include overstepping time steps
					echo -e ${line} ${nevent} >> ${event_file_backward}
					overstep_counter=0
				fi

				#Eliminate line from file (if present)
				awk -v year=${year_n} -v month=${month_n} -v day=${day_n} -v hour=${time_n} 'substr($1,1,4)!=year || substr($1,6,2)!=month || substr($1,9,2)!=day || substr($1,12,2)!=hour'  ${flux_list_sorted} > shorter_file
				mv shorter_file ${flux_list_sorted}
			done	#End BACKWARD loop

			#Merge results in a unique file and attach it to file with all events
			tac ${event_file_backward} >> ${event_file}
			cat ${event_file_forward} >> ${event_file}
			cat ${event_file} >> ${events_list_file}

			#Compute length of the event and add event to summary file with maxima
			length_event=$(wc -l ${event_file}  | awk '{print $1*6}') # multiply by 6 to convert from lines to hours
			echo -e ${line_summary} ${nevent} ${length_event} >> ${events_list_file_summary}

			#Remove files
			\rm -f ${event_file_backward} ${event_file} ${event_file_forward} overstep_file check_file

			#Check value of heat flux for while loop
			fluxmax_val=$(tail -1 ${flux_list_sorted} | awk '{print $2}')

		done	#End first loop
		\rm -f ${event_file_backward} ${event_file} ${event_file_forward} overstep_file check_file
		\rm -f ${flux_list_file} ${flux_list_sorted}_backup ${flux_list_sorted}

		awk '{OFS=" "; $2 = $2 * -1.0; print $0}' ${events_list_file} > ${events_list_file}.tmp && mv ${events_list_file}.tmp ${events_list_file}
		awk '{OFS=" "; $2 = $2 * -1.0; print $0}' ${events_list_file_summary} > ${events_list_file_summary}.tmp && mv ${events_list_file_summary}.tmp ${events_list_file_summary}
		awk '$4>24' ${events_list_file_summary} > ${events_list_file_summary}_filtered # Filter out events longer than 24 hours
	done	#End loop through timeframes
done 	#End loop through regions
exit

