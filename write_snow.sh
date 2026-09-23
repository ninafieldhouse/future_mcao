#!/bin/bash

set -euo pipefail

output_dir="/net/litho/atmosdyn2/nfieldhouse/data/cao/clim"
member_ids=(1000 1100 1200 1300 1400)

mkdir -p "$output_dir"

select_djf_files() {
	local input_dir="$1"
	local file_prefix="$2"
	local first_december_year="$3"
	local last_february_year="$4"
	local year

	for ((year=first_december_year; year<last_february_year; year++)); do
		printf '%s/%s.%04d-12.nc\n' "$input_dir" "$file_prefix" "$year"
		next_year=$((year + 1))
		printf '%s/%s.%04d-01.nc\n' "$input_dir" "$file_prefix" "$next_year"
		printf '%s/%s.%04d-02.nc\n' "$input_dir" "$file_prefix" "$next_year"
	done
}

write_member_mean() {
	local member_id="$1"
	local timeframe="$2"
	local output_file="$3"
	local cesm_type
	local first_december_year
	local last_february_year
	local input_dir
	local file_prefix
	local temp_dir
	local -a source_files=()
	local -a selected_files=()
	local source_file
	local selected_file
	local index=0

	if [[ "$timeframe" == "present" ]]; then
		cesm_type="HIST"
		first_december_year=1990
		last_february_year=2000
	else
		cesm_type="SSP370"
		first_december_year=2090
		last_february_year=2100
	fi

	input_dir="/net/meso/climphys/cesm212/b.e212.B${cesm_type}cmip6.f09_g17.${member_id}/archive/lnd/hist"
	file_prefix="b.e212.B${cesm_type}cmip6.f09_g17.${member_id}.clm2.h0"
	mapfile -t source_files < <(select_djf_files "$input_dir" "$file_prefix" "$first_december_year" "$last_february_year")

	for source_file in "${source_files[@]}"; do
		if [[ ! -f "$source_file" ]]; then
			echo "Missing input file: $source_file" >&2
			return 1
		fi
	done

	temp_dir=$(mktemp -d)
	trap 'rm -rf "$temp_dir"' RETURN

	for source_file in "${source_files[@]}"; do
		selected_file=$(printf '%s/selected_%02d.nc' "$temp_dir" "$index")
		cdo -L -selname,FSNO,SNODP "$source_file" "$selected_file"
		selected_files+=("$selected_file")
		index=$((index + 1))
	done

	echo "Writing ${timeframe} DJF mean for member ${member_id} to ${output_file}"
	cdo -L -timmean -mergetime "${selected_files[@]}" "$output_file"

	trap - RETURN
	rm -rf "$temp_dir"
}

write_timeframe_mean() {
	local timeframe="$1"
	local first_december_year
	local last_february_year
	local final_output
	local temp_dir
	local -a member_outputs=()
	local member_id

	if [[ "$timeframe" == "present" ]]; then
		first_december_year=1990
		last_february_year=2000
	else
		first_december_year=2090
		last_february_year=2100
	fi

	final_output="${output_dir}/snow_djf_mean_${timeframe}_${first_december_year}_${last_february_year}.nc"
	temp_dir=$(mktemp -d)
	trap 'rm -rf "$temp_dir"' RETURN

	for member_id in "${member_ids[@]}"; do
		member_outputs+=("${temp_dir}/${timeframe}_${member_id}.nc")
		write_member_mean "$member_id" "$timeframe" "${temp_dir}/${timeframe}_${member_id}.nc"
	done

	echo "Writing ensemble-mean ${timeframe} DJF mean to ${final_output}"
	cdo -L ensmean "${member_outputs[@]}" "$final_output"

	trap - RETURN
	rm -rf "$temp_dir"
}

write_timeframe_mean present
write_timeframe_mean future

echo "Done."
