#!/bin/bash
. /etc/profile.d/modules.sh

# -------------------------------------------------------
# configuration
# -------------------------------------------------------

. /home/svoigt/Master_Thesis/config_traj_barents_cesm2.sh

# -------------------------------------------------------
# set up configuration
# -------------------------------------------------------

# load modules and set environmental variables
module unload dyn_tools
module load dyn_tools
cesm
NETCDF_FORMAT=CF;export $NETCDF_FORMAT
LAGRANTO=/usr/local/dyn_tools//lagranto.cesm;export LAGRANTO
MODEL=cesm;export MODEL

# create output directory
#mkdir -pv ${destdir} 

# go to runnning directory
#if [ ${cre_new_rundir} -eq 1 ]; then
#    rm -r ${rundir} &> /dev/null
#fi
#mkdir -pv ${rundir} &> /dev/null
cd ${rundir} &> /dev/null

# -------------------------------------------------------
# define cesm time period an member setup
# -------------------------------------------------------
mems=( 1100 ) #002 003 004 )
styrs=( 2100 )
enyrs=( 2100 )
#styrs=( 1992 1996 2000 )
#enyrs=( 1992 1996 2000 )
stday=( 0101 ) #0201 0301 0401 0501 0601 0701 0801 0901 1001 1101 1201 )
enday=( 0131 ) #0228 0331 0430 0531 0630 0731 0831 0930 1031 1130 1231 )
flag="cesm2"
string1='-01-01-00000'
string1l='-01-02-00000'
subdir='archive/atm/hist/'
subdirl='archive/lnd/hist/'

# -------------------------------------------------------
# define procedures
# -------------------------------------------------------
# function to link CESM netcdf files into running directory for current process
# (${rundir}/proc${iproc})
crefiles() {
    
    . /home/svoigt/Master_Thesis/config_traj_barents_cesm2.sh

    # the function receives the date and process ids as arguments
    datum=$1
    memnr=$2
    iproc=$3
    dt_end=$4
    
    echo "  --> Prepare files for "${datum}
    
    # go to running directory
    cd ${rundir}/proc${iproc}
    
    # determine time range of all files needed for trajectory computation
    # first option: forward
    # second option: backward
    if [ ${traj_len} -gt 0 ]; then
        tstart=`echo "${dt_start} - ${dt_input}" | bc`
        tend=`echo "${dt_end} + ${traj_len} + ${dt_input}" | bc`
    else
        tstart=`echo "${dt_start} + ${traj_len} - ${dt_input}" | bc`
        tend=`echo "${dt_end} + ${dt_input}" | bc`        
    fi

    # link all files into run directory
    for j in $(seq ${tstart} ${dt_input} ${tend})
    do

        datum2=`newtime ${datum} ${j}`

        yr2=`echo ${datum2} | cut -c1-4`
        mn2=`echo ${datum2} | cut -c5-6`
        day2=`echo ${datum2} | cut -c7-8`
        hr2=`echo ${datum2} | cut -c10-11`

        # --------- #
        #  link/prepare files...
        # --------- #
        string0='.cam.h3.'
        string0s='.cam.h4.'        
        string0l='.clm2.h1.'
        #let prevyear=$year-1 #change sven
        year=`echo ${datum} | cut -c1-4`
        let prevyear=$year-1
        if [ ${year} -ge 2100 ]; then
            yr2=$prevyear
        fi
        
        if [ ${yr2} -lt 2015 ]; then
            premem='b.e212.BHISTcmip6.f09_g17.'
        else
            premem='b.e212.BSSP370cmip6.f09_g17.'
        fi

        PFILE=$premem$memnr$string0$yr2$string1'.nc'
        SFILE=$premem$memnr$string0s$yr2$string1'.nc'
        #SFILE='S'$premem$memnr$string0$yr2$string1'.nc'
        CFILE='mincl'$yr2$'.nc'
        LFILE=$premem$memnr$string0l$yr2$string1l'.nc'
        #link P-files
        if [ ! -f P${yr2} ]; then
            ln -s $cesmdir'/'$premem$memnr'/'$subdir'/'$PFILE './P'$yr2 &> /dev/null
        fi
        #link S-files
        if [ ! -f S${yr2} ]; then
        	ncks -A -v datesec,date,date_written,lev $cesmdir'/'$premem$memnr'/'$subdir'/'$PFILE './S'$yr2
        	ncks -A -v ICEFRAC,QREFHT,TREFHT,U10 $cesmdir'/'$premem$memnr'/'$subdir'/'$SFILE './S'$yr2 
        fi
        
        if [ ! -f P${yr2} ]; then
            echo "  --> error preparing files: "${datum2}
            echo "exit"
            exit
        fi

    done
    
    return 0

}
export -f crefiles


# function to calculate trajectories for date ${datum}
# requires date and process id as input, e.g., 
# caltraj 20180101_00 1
caltraj() {
    
    . /home/svoigt/Master_Thesis/config_traj_barents_cesm2.sh
    
    # setup of environment for trajectory calculations
    module unload dyn_tools
    module load dyn_tools
    cesm
    NETCDF_FORMAT=CF;export $NETCDF_FORMAT
    LAGRANTO=/usr/local/dyn_tools//lagranto.cesm;export LAGRANTO
    MODEL=cesm;export MODEL
 
    # the function receives the date and process ids as arguments
    datum=$1
    memnr=$2
    iproc=$3
    dt_end=$4
        
    # create run directory for current process
    mkdir -pv ${rundir}/proc${iproc} &> /dev/null
    cd ${rundir}/proc${iproc} &> /dev/null
    
    # output
    echo "(job-id = ${iproc}) Compute trajectories for trackid ${trackid} at ${datum}"
    
    # call function to link all input files into run directory
    crefiles ${datum} ${memnr} ${iproc} ${dt_end}
    
    # ------------------------------------------------------------------
    # 1. compute trajectories
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # 1.1 create start files
    # ------------------------------------------------------------------ 

    # variables to be traced
    # for lagranto fortran, the file must have the name tracevars
    echo "PRECL 1000000. P" >& tracevars
    echo "PRECC 1000000. P" >> tracevars    
        
    # loop over all trajectory initialization times
    for tdiff in $(seq ${dt_start} ${dt_ini} ${dt_end})
    do
    
        # compute current date
        curdate=`newtime ${datum} ${tdiff}`
        
        year=`echo ${datm} | cut -c1-4` # sven
        year2=`echo ${curdate} | cut -c1-4`
        let prevyear=$year-1
        if [ ${year2} -gt 2150 ]; then
            mn2=`echo ${datum2} | cut -c5-6`
            day2=`echo ${datum2} | cut -c7-8`
            hr2=`echo ${datum2} | cut -c10-11`
            curdate=$prevyear$mn2$day2_$hr2
        fi 
    
        # create start file
        #ln -s $mskdir/$mskfile start_mask.nc
        if [ ! -f mask_bk ]; then
        	ln -s $mskdir/$mskfile mask_bk
        fi
        echo 'This is the date: '${curdate}
        echo ${plev_list}
        #startf ${curdate} startf 'box.grid(20,60,70,80)@list('${plev_list}')@'${plev_unit}
        startf ${curdate} startf 'file(mask_bk)@list('${plev_list}')@'${plev_unit}
        #startf ${curdate} startf 'mask.grid(start_mask.nc,MSK)@list('${plev_list}')@'${plev_unit}  # -criterion GT:0:PRECL | GT:0:PRECC -t tracefile
                
        caltra startf 6 outtrace.1 -ref ${curdate}
        
        select.cesm outtrace.1 startf.1 'GT:PRECC:0.01:FIRST | GT:PRECL:0.01:FIRST' -startf

        
        #startf ${curdate} startf 'box.grid(20,60,70,80)@list('${plev_list}')@'${plev_unit}'@GT:0.000000001:PRECC|GT:0.000000001:PRECC'


                
        tail -n+5 startf.1 | awk -v CURTDIFF="${tdiff}" '{if (NF > 1) {printf "%9.2f %9.2f %9.2f %9.2f\n", CURTDIFF, $2, $3, $4}}' >> startf_all
       
        # clean up
        rm startf # startf_selected.1 
        rm startf.1 #
        rm outtrace.1
                     
    done
    
    # ------------------------------------------------------------------
    # 1.2 compute trajectories
    # ------------------------------------------------------------------
    
    # Here the faster fortran verions for lagranto is used. This Version has the 
    # advantage that the I/O is reduced because input netcdf files are read only 
    # once. The disadvantage is that trajectories with different starting times 
    # are all in the same file.

    # add variables to be traced
    echo "Q 1000. P" >> tracevars
    echo "T 1. P" >> tracevars
    echo "SST 1. P" >> tracevars    
    echo "PS .01 P" >> tracevars    
    echo "PSL 0.01 P" >> tracevars    
    echo "LHFLX 1. P" >> tracevars
    echo "QRAIN 1000. P" >> tracevars
    echo "QSNOW 1000. P" >> tracevars
    echo "ICEFRAC 1. S" >> tracevars
    echo "QREFHT 1000. S" >> tracevars
    echo "TREFHT 1. S" >> tracevars
    echo "U10 1. S" >> tracevars  

        
    # input timestep in minutes
    dt_input_min=`echo "${dt_input} * 60" | bc -l`
    
    # run trajectories    
    caltra startf_all ${traj_len} traj${datum}.1 -i ${dt_input_min} -o ${dt_input_min} -ts 30 -ref ${datum} -$flag
        
    # clean up
    rm tracevars 
    
    # save data
    mv startf_all ${destdir}/startf_all_${datum}.${memnr}
    mv traj${datum}.1 ${destdir}/traj${datum}.${memnr}
    
    # clean up
    if [ ${clean_rundir} -eq 1 ]; then
        cd ${rundir}
        rm -r proc${iproc}
    fi

    return 0
    
}
export -f caltraj


# ------------------------------------------------------------------
# Execution of program
# ------------------------------------------------------------------

for y in $(seq 0 1 $(( ${#styrs[@]} - 1 )) ); do

styr=${styrs[$y]}
enyr=${enyrs[$y]}

for memnr in ${mems[@]}; do

for yr in $(seq $styr 1 $enyr);do

for std in $(seq 0 1 $(( ${#stday[@]} - 1 )) ); do

stdate=${yr}${stday[$std]}_00
endate=${yr}${enday[$std]}_18
deltat=`gettidiff $endate $stdate`
# for h in $(seq 0 24 ${deltat}); do
# 
# nowdate=`newtime $stdate $h`

caltraj $stdate $memnr 4 $deltat

# done

done

done

done

done

# for parallel calculation for all dates in startdates array
#parallel -j ${npar} --lb 'caltraj {1} {#}' ::: ${startdates[@]}


exit 0

