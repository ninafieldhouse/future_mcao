# program directory
scriptdir=$PWD

#destination directory
destdir=/net/tropo/atmosdyn/svoigt/lagranto

# cesm data directory 
cesmdir=/net/meso/climphys/cesm212

# mask dir
mskdir=/net/tropo/atmosdyn/svoigt/lagranto
mskfile=mask_barent_kara.csv
mskvar=MSK

# rundir
# for each process, a separate run directory is created within rundir
# with the name proc${iproc} where iproc is the process id
rundir=${destdir}/run_cesm2

# parallel computing units
npar=1

# new run directory
cre_new_rundir=0

# clean up after each task (should be on by default)
clean_rundir=0		#1=default

# ---------------------------------------------------------------------
# options for trajectory calculations
# ---------------------------------------------------------------------

# length of trajectories [h]
# positive values: forward trajectories
# negative values: backward trajectories
traj_len=-240

# time increment of input data [h]
dt_input=6

# start dates 
# trajectories are computed relative to start dates in time range specified by dt_start and dt_end
dt_ini=6                    # time interval between initialization times
dt_start=0                  # first date relative to start date
#dt_end=738           	    # last date relative to start date (90 days: 2160h)

# spacing of regular starting grid [km]
dx_traj=80   

# trajectory starting region
# KAN_L 67.0955 49.9513 670 2008-09-01 2020
# KAN_M 67.0670 48.8355 1270 2008-09-02 2020
# KAN_U 67.0003 47.0253 1840 2009-04-04 2020
# BOX: -50.15, -48, 67.00, 67.25
# lon1=-50.15
# lon2=-48
# lat1=67
# lat2=67.25

# vertical levels of regular starting grid    
dp=50 
pmax=550
pmin=10
plev=$(seq ${pmin} ${dp} ${pmax})
plev_unit='hPa,agl'


# ---------------------------------------------------------------------
# some definitions
# ---------------------------------------------------------------------

plev_list=""

for curplev in ${plev[@]}
do
    plev_list=${plev_list}${curplev}","
done
plev_list=${plev_list::-1}


# mass of trajectory
trajmass=`echo "${dp} * 10^2 * ${dx_traj} * ${dx_traj} * 10^6 / 9.81" | bc -l`




