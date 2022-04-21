#!/usr/bin/ksh

# -------------------------------------------------------------------------------
# Important notes
# -------------------------------------------------------------------------------

# - Digital elevation model (DEM)
#   GLOBE data suffers from serious artefacts in some regions. E.g. in High 
#   Mountain Asia, deviations exceeding 1000 m were found relative to more 
#   accurate and recent DEM data. Regions in the tropics seem to be affected 
#   as well. It is generally recommended to use the more accurate and recent 
#   MERIT data set.
# - Surface albedo
#   Note the correct pairing of albedo settings in EXTPAR ('ialb_type') and 
#   COSMO ('itype_albedo'). A valid combination that revealed good results is 
#   'ialb_type = 1' and 'itype_albedo = 3'.
# - Sequential and parallel run of EXTPAR
#   For most cases, EXTPAR can be run in parallel. However, in case MERIT
#   data is used ('itopo type = 3') for large spatial domains, memory issues
#   can occur. These issues can be avoided with a sequential run
#   ('run_parallel = false').

# -------------------------------------------------------------------------------

# Import functions to launch Extpar executables
. ../test/testsuite/bin/runcontrol_functions.sh

ulimit -s unlimited
ulimit -c unlimited

# Get hostname
hostname="`echo $HOSTNAME`"
logfile="extpar_runscript.log"

# -------------------------------------------------------------------------------
# Options to define by user
# -------------------------------------------------------------------------------

# Output file format and name
netcdf_output_filename='extpar_EAS_ext_12km_merit.nc'

# Output directory (-> ensure that this location provides enough disk space)
sandboxdir=/scratch/snx3000/csteger/EXTPAR/output/EAS_12km_1118x670

# Parallel or sequential EXTPAR run
run_parallel=false

# -------------------------------------------------------------------------------
# Auto-set paths
# -------------------------------------------------------------------------------

# Directory of runscripts -> need to be as in original repository
scriptdir=`pwd`
src_python=${scriptdir}/../python/lib

# Change dir to src_python to get absolute path
cd $src_python
export PYTHONPATH=$PYTHONPATH:$(pwd)
cd - > /dev/null 2>&1

# Directory of compiled extpar executables
exedir=$scriptdir/../bin

# -------------------------------------------------------------------------------
# Host-dependent paths and variables
# -------------------------------------------------------------------------------

# CSCS-machines
if [[ $hostname == tsa* || $hostname == arolla* || $hostname == nid* \
	|| $hostname == daint* ]]; then

    # NetCDF raw data for external parameter
    #data_dir=/store/c2sm/extpar_raw_data/linked_data
    data_dir=/scratch/snx3000/csteger/EXTPAR/input_linked

else

    # Exit script in case of unknown host
    echo ERROR: Unknown host: $hostname >> ${logfile}
    exit 1

fi

# -------------------------------------------------------------------------------
# Names of Python/Fortran executables
# -------------------------------------------------------------------------------

# Python
binary_alb=extpar_alb_to_buffer.py
binary_ndvi=extpar_ndvi_to_buffer.py
binary_tclim=extpar_cru_to_buffer.py

# Fortran
binary_lu=extpar_landuse_to_buffer.exe
binary_topo=extpar_topo_to_buffer.exe
binary_aot=extpar_aot_to_buffer.exe
binary_soil=extpar_soil_to_buffer.exe
binary_flake=extpar_flake_to_buffer.exe
binary_consistency_check=extpar_consistency_check.exe

# -------------------------------------------------------------------------------
# Names of input fields 
# -------------------------------------------------------------------------------

# Albedo
raw_data_alb='alb_new.nc'
raw_data_alnid='alnid_new.nc'
raw_data_aluvd='aluvd_new.nc'
buffer_alb='month_alb_buffer.nc'

# Aerosols
raw_data_aot='aod_AeroCom1.nc'
buffer_aot='extpar_buffer_aot.nc'

# 2 m air temperature climatology
raw_data_tclim_fine='CRU_T_SOIL_clim.nc'
buffer_tclim='crutemp_clim_extpar_buffer.nc'

# GLCC land use
raw_data_glcc='GLCC_usgs_class_byte.nc'
buffer_glcc='glcc_landuse_buffer.nc'

# Globcover 2009 land use
raw_data_globcover_0='GLOBCOVER_0_16bit.nc'
raw_data_globcover_1='GLOBCOVER_1_16bit.nc'
raw_data_globcover_2='GLOBCOVER_2_16bit.nc'
raw_data_globcover_3='GLOBCOVER_3_16bit.nc'
raw_data_globcover_4='GLOBCOVER_4_16bit.nc'
raw_data_globcover_5='GLOBCOVER_5_16bit.nc'

buffer_lu='extpar_landuse_buffer.nc'

# MERIT digital elevation model
raw_data_merit_N90_N60_W180_W150='MERIT_N90-N60_W180-W150.nc'
raw_data_merit_N90_N60_W150_W120='MERIT_N90-N60_W150-W120.nc'
raw_data_merit_N90_N60_W120_W090='MERIT_N90-N60_W120-W090.nc'
raw_data_merit_N90_N60_W090_W060='MERIT_N90-N60_W090-W060.nc'
raw_data_merit_N90_N60_W060_W030='MERIT_N90-N60_W060-W030.nc'
raw_data_merit_N90_N60_W030_E000='MERIT_N90-N60_W030-E000.nc'
raw_data_merit_N90_N60_E000_E030='MERIT_N90-N60_E000-E030.nc'
raw_data_merit_N90_N60_E030_E060='MERIT_N90-N60_E030-E060.nc'
raw_data_merit_N90_N60_E060_E090='MERIT_N90-N60_E060-E090.nc'
raw_data_merit_N90_N60_E090_E120='MERIT_N90-N60_E090-E120.nc'
raw_data_merit_N90_N60_E120_E150='MERIT_N90-N60_E120-E150.nc'
raw_data_merit_N90_N60_E150_E180='MERIT_N90-N60_E150-E180.nc'
raw_data_merit_N60_N30_W180_W150='MERIT_N60-N30_W180-W150.nc'
raw_data_merit_N60_N30_W150_W120='MERIT_N60-N30_W150-W120.nc'
raw_data_merit_N60_N30_W120_W090='MERIT_N60-N30_W120-W090.nc'
raw_data_merit_N60_N30_W090_W060='MERIT_N60-N30_W090-W060.nc'
raw_data_merit_N60_N30_W060_W030='MERIT_N60-N30_W060-W030.nc'
raw_data_merit_N60_N30_W030_E000='MERIT_N60-N30_W030-E000.nc'
raw_data_merit_N60_N30_E000_E030='MERIT_N60-N30_E000-E030.nc'
raw_data_merit_N60_N30_E030_E060='MERIT_N60-N30_E030-E060.nc'
raw_data_merit_N60_N30_E060_E090='MERIT_N60-N30_E060-E090.nc'
raw_data_merit_N60_N30_E090_E120='MERIT_N60-N30_E090-E120.nc'
raw_data_merit_N60_N30_E120_E150='MERIT_N60-N30_E120-E150.nc'
raw_data_merit_N60_N30_E150_E180='MERIT_N60-N30_E150-E180.nc'
raw_data_merit_N30_N00_W180_W150='MERIT_N30-N00_W180-W150.nc'
raw_data_merit_N30_N00_W150_W120='MERIT_N30-N00_W150-W120.nc'
raw_data_merit_N30_N00_W120_W090='MERIT_N30-N00_W120-W090.nc'
raw_data_merit_N30_N00_W090_W060='MERIT_N30-N00_W090-W060.nc'
raw_data_merit_N30_N00_W060_W030='MERIT_N30-N00_W060-W030.nc'
raw_data_merit_N30_N00_W030_E000='MERIT_N30-N00_W030-E000.nc'
raw_data_merit_N30_N00_E000_E030='MERIT_N30-N00_E000-E030.nc'
raw_data_merit_N30_N00_E030_E060='MERIT_N30-N00_E030-E060.nc'
raw_data_merit_N30_N00_E060_E090='MERIT_N30-N00_E060-E090.nc'
raw_data_merit_N30_N00_E090_E120='MERIT_N30-N00_E090-E120.nc'
raw_data_merit_N30_N00_E120_E150='MERIT_N30-N00_E120-E150.nc'
raw_data_merit_N30_N00_E150_E180='MERIT_N30-N00_E150-E180.nc'
raw_data_merit_N00_S30_W180_W150='MERIT_N00-S30_W180-W150.nc'
raw_data_merit_N00_S30_W150_W120='MERIT_N00-S30_W150-W120.nc'
raw_data_merit_N00_S30_W120_W090='MERIT_N00-S30_W120-W090.nc'
raw_data_merit_N00_S30_W090_W060='MERIT_N00-S30_W090-W060.nc'
raw_data_merit_N00_S30_W060_W030='MERIT_N00-S30_W060-W030.nc'
raw_data_merit_N00_S30_W030_E000='MERIT_N00-S30_W030-E000.nc'
raw_data_merit_N00_S30_E000_E030='MERIT_N00-S30_E000-E030.nc'
raw_data_merit_N00_S30_E030_E060='MERIT_N00-S30_E030-E060.nc'
raw_data_merit_N00_S30_E060_E090='MERIT_N00-S30_E060-E090.nc'
raw_data_merit_N00_S30_E090_E120='MERIT_N00-S30_E090-E120.nc'
raw_data_merit_N00_S30_E120_E150='MERIT_N00-S30_E120-E150.nc'
raw_data_merit_N00_S30_E150_E180='MERIT_N00-S30_E150-E180.nc'

buffer_topo='topography_buffer.nc'
output_topo='topography_COSMO.nc'

# Names of generated SGSL files
raw_data_sgsl_N90_N60_W180_W150='S_ORO_N90-N60_W180-W150.nc'
raw_data_sgsl_N90_N60_W150_W120='S_ORO_N90-N60_W150-W120.nc'
raw_data_sgsl_N90_N60_W120_W090='S_ORO_N90-N60_W120-W090.nc'
raw_data_sgsl_N90_N60_W090_W060='S_ORO_N90-N60_W090-W060.nc'
raw_data_sgsl_N90_N60_W060_W030='S_ORO_N90-N60_W060-W030.nc'
raw_data_sgsl_N90_N60_W030_E000='S_ORO_N90-N60_W030-E000.nc'
raw_data_sgsl_N90_N60_E000_E030='S_ORO_N90-N60_E000-E030.nc'
raw_data_sgsl_N90_N60_E030_E060='S_ORO_N90-N60_E030-E060.nc'
raw_data_sgsl_N90_N60_E060_E090='S_ORO_N90-N60_E060-E090.nc'
raw_data_sgsl_N90_N60_E090_E120='S_ORO_N90-N60_E090-E120.nc'
raw_data_sgsl_N90_N60_E120_E150='S_ORO_N90-N60_E120-E150.nc'
raw_data_sgsl_N90_N60_E150_E180='S_ORO_N90-N60_E150-E180.nc'
raw_data_sgsl_N60_N30_W180_W150='S_ORO_N60-N30_W180-W150.nc'
raw_data_sgsl_N60_N30_W150_W120='S_ORO_N60-N30_W150-W120.nc'
raw_data_sgsl_N60_N30_W120_W090='S_ORO_N60-N30_W120-W090.nc'
raw_data_sgsl_N60_N30_W090_W060='S_ORO_N60-N30_W090-W060.nc'
raw_data_sgsl_N60_N30_W060_W030='S_ORO_N60-N30_W060-W030.nc'
raw_data_sgsl_N60_N30_W030_E000='S_ORO_N60-N30_W030-E000.nc'
raw_data_sgsl_N60_N30_E000_E030='S_ORO_N60-N30_E000-E030.nc'
raw_data_sgsl_N60_N30_E030_E060='S_ORO_N60-N30_E030-E060.nc'
raw_data_sgsl_N60_N30_E060_E090='S_ORO_N60-N30_E060-E090.nc'
raw_data_sgsl_N60_N30_E090_E120='S_ORO_N60-N30_E090-E120.nc'
raw_data_sgsl_N60_N30_E120_E150='S_ORO_N60-N30_E120-E150.nc'
raw_data_sgsl_N60_N30_E150_E180='S_ORO_N60-N30_E150-E180.nc'
raw_data_sgsl_N30_N00_W180_W150='S_ORO_N30-N00_W180-W150.nc'
raw_data_sgsl_N30_N00_W150_W120='S_ORO_N30-N00_W150-W120.nc'
raw_data_sgsl_N30_N00_W120_W090='S_ORO_N30-N00_W120-W090.nc'
raw_data_sgsl_N30_N00_W090_W060='S_ORO_N30-N00_W090-W060.nc'
raw_data_sgsl_N30_N00_W060_W030='S_ORO_N30-N00_W060-W030.nc'
raw_data_sgsl_N30_N00_W030_E000='S_ORO_N30-N00_W030-E000.nc'
raw_data_sgsl_N30_N00_E000_E030='S_ORO_N30-N00_E000-E030.nc'
raw_data_sgsl_N30_N00_E030_E060='S_ORO_N30-N00_E030-E060.nc'
raw_data_sgsl_N30_N00_E060_E090='S_ORO_N30-N00_E060-E090.nc'
raw_data_sgsl_N30_N00_E090_E120='S_ORO_N30-N00_E090-E120.nc'
raw_data_sgsl_N30_N00_E120_E150='S_ORO_N30-N00_E120-E150.nc'
raw_data_sgsl_N30_N00_E150_E180='S_ORO_N30-N00_E150-E180.nc'
raw_data_sgsl_N00_S30_W180_W150='S_ORO_N00-S30_W180-W150.nc'
raw_data_sgsl_N00_S30_W150_W120='S_ORO_N00-S30_W150-W120.nc'
raw_data_sgsl_N00_S30_W120_W090='S_ORO_N00-S30_W120-W090.nc'
raw_data_sgsl_N00_S30_W090_W060='S_ORO_N00-S30_W090-W060.nc'
raw_data_sgsl_N00_S30_W060_W030='S_ORO_N00-S30_W060-W030.nc'
raw_data_sgsl_N00_S30_W030_E000='S_ORO_N00-S30_W030-E000.nc'
raw_data_sgsl_N00_S30_E000_E030='S_ORO_N00-S30_E000-E030.nc'
raw_data_sgsl_N00_S30_E030_E060='S_ORO_N00-S30_E030-E060.nc'
raw_data_sgsl_N00_S30_E060_E090='S_ORO_N00-S30_E060-E090.nc'
raw_data_sgsl_N00_S30_E090_E120='S_ORO_N00-S30_E090-E120.nc'
raw_data_sgsl_N00_S30_E120_E150='S_ORO_N00-S30_E120-E150.nc'
raw_data_sgsl_N00_S30_E150_E180='S_ORO_N00-S30_E150-E180.nc'

# NDVI
raw_data_ndvi='NDVI_1998_2003.nc'
buffer_ndvi='ndvi_buffer.nc'

# Soil
raw_data_soil_HWSD='HWSD0_30_topsoil.nc'
raw_data_deep_soil='HWSD30_100_subsoil.nc'
buffer_soil='soil_buffer.nc'
raw_lookup_table_HWSD='LU_TAB_HWSD_UF.data'
raw_HWSD_data='HWSD_DATA_COSMO.data'
raw_HWSD_data_deep='HWSD_DATA_COSMO_S.data'
raw_HWSD_data_extpar='HWSD_DATA_COSMO_EXTPAR.asc'

# Lakes
raw_data_flake='GLDB_lakedepth.nc'
buffer_flake='flake_buffer.nc'

# -------------------------------------------------------------------------------
# Prepare working directory
# -------------------------------------------------------------------------------

if [[ ! -d ${sandboxdir} ]] ; then
  mkdir -p ${sandboxdir}
fi

# Link raw data files to local workdir
ln -s -f ${data_dir}/*.nc ${sandboxdir}

cp $scriptdir/* ${sandboxdir}/.
cp $exedir/* ${sandboxdir}/.
cd ${sandboxdir}
if [[ -e ${logfile} ]] ; then
  rm -f ${logfile}
fi

cd ${sandboxdir}

echo "\n>>>> Data will be processed and produced in `pwd` <<<<\n"

echo PYTHONPATH: ${PYTHONPATH} >> ${logfile}

# -------------------------------------------------------------------------------
# Create input namelists 
# -------------------------------------------------------------------------------

cat > namelist.py << EOF_namelist_python
input_alb = {
        'ialb_type': 1,
        'raw_data_alb_path': '',
        'raw_data_alb_filename': '${raw_data_alb}',
        'raw_data_alnid_filename': '${raw_data_alnid}',
        'raw_data_aluvd_filename': '${raw_data_aluvd}',
        'alb_buffer_file': '${buffer_alb}',
        }

input_tclim = {
        'raw_data_t_clim_path': '',
        'raw_data_tclim_coarse': '',
        'raw_data_tclim_fine': '${raw_data_tclim_fine}',
        't_clim_buffer_file': '${buffer_tclim}',
        'it_cl_type': 1
        }

input_ndvi = {
        'raw_data_ndvi_path': '',
        'raw_data_ndvi_filename': '${raw_data_ndvi}',
        'ndvi_buffer_file': '${buffer_ndvi}',
        }
EOF_namelist_python

# Set target grid definition 
cat > INPUT_grid_org << EOF_go
&GRID_DEF 
 igrid_type = 2,
 domain_def_namelist='INPUT_COSMO_GRID'
/ 
EOF_go
cat > INPUT_COSMO_GRID << EOF_grid
&lmgrid
 pollon=-63.70, 
 pollat=61.00, 
 startlon_tot=-58.05,
 startlat_tot=-27.69,
 dlon=0.11,
 dlat=0.11,
 ie_tot=1118,
 je_tot=670,
/
EOF_grid

cat > INPUT_AOT << EOF_aot
&aerosol_raw_data
  raw_data_aot_path='',
  raw_data_aot_filename='${raw_data_aot}'
  iaot_type=2
/  
&aerosol_io_extpar
  aot_buffer_file='${buffer_aot}',
/
EOF_aot

cat > INPUT_LU << EOF_lu
&lu_raw_data
   raw_data_lu_path='',
   raw_data_lu_filename = \
'${raw_data_globcover_0}' '${raw_data_globcover_1}' \
'${raw_data_globcover_2}' '${raw_data_globcover_3}' \
'${raw_data_globcover_4}' '${raw_data_globcover_5}',
   i_landuse_data=1,
   ilookup_table_lu=1
/
&lu_io_extpar
   lu_buffer_file='${buffer_lu}',
/
&glcc_raw_data
   raw_data_glcc_path='',
   raw_data_glcc_filename='${raw_data_glcc}'
/
&glcc_io_extpar
   glcc_buffer_file='${buffer_glcc}',
/
EOF_lu

cat > INPUT_ORO << EOF_oro
&oro_runcontrol
  lcompute_sgsl=.TRUE. ,
  /
&orography_io_extpar
  orography_buffer_file='${buffer_topo}',
  orography_output_file='${output_topo}'
/
&orography_raw_data
 itopo_type = 3,
 lsso_param = .TRUE.,
 raw_data_orography_path='',
 ntiles_column = 12,
 ntiles_row = 4,
 topo_files = \
'${raw_data_merit_N90_N60_W180_W150}'
'${raw_data_merit_N90_N60_W150_W120}'
'${raw_data_merit_N90_N60_W120_W090}'
'${raw_data_merit_N90_N60_W090_W060}'
'${raw_data_merit_N90_N60_W060_W030}'
'${raw_data_merit_N90_N60_W030_E000}'
'${raw_data_merit_N90_N60_E000_E030}'
'${raw_data_merit_N90_N60_E030_E060}'
'${raw_data_merit_N90_N60_E060_E090}'
'${raw_data_merit_N90_N60_E090_E120}'
'${raw_data_merit_N90_N60_E120_E150}'
'${raw_data_merit_N90_N60_E150_E180}'
'${raw_data_merit_N60_N30_W180_W150}'
'${raw_data_merit_N60_N30_W150_W120}'
'${raw_data_merit_N60_N30_W120_W090}'
'${raw_data_merit_N60_N30_W090_W060}'
'${raw_data_merit_N60_N30_W060_W030}'
'${raw_data_merit_N60_N30_W030_E000}'
'${raw_data_merit_N60_N30_E000_E030}'
'${raw_data_merit_N60_N30_E030_E060}'
'${raw_data_merit_N60_N30_E060_E090}'
'${raw_data_merit_N60_N30_E090_E120}'
'${raw_data_merit_N60_N30_E120_E150}'
'${raw_data_merit_N60_N30_E150_E180}'
'${raw_data_merit_N30_N00_W180_W150}'
'${raw_data_merit_N30_N00_W150_W120}'
'${raw_data_merit_N30_N00_W120_W090}'
'${raw_data_merit_N30_N00_W090_W060}'
'${raw_data_merit_N30_N00_W060_W030}'
'${raw_data_merit_N30_N00_W030_E000}'
'${raw_data_merit_N30_N00_E000_E030}'
'${raw_data_merit_N30_N00_E030_E060}'
'${raw_data_merit_N30_N00_E060_E090}'
'${raw_data_merit_N30_N00_E090_E120}'
'${raw_data_merit_N30_N00_E120_E150}'
'${raw_data_merit_N30_N00_E150_E180}'
'${raw_data_merit_N00_S30_W180_W150}'
'${raw_data_merit_N00_S30_W150_W120}'
'${raw_data_merit_N00_S30_W120_W090}'
'${raw_data_merit_N00_S30_W090_W060}'
'${raw_data_merit_N00_S30_W060_W030}'
'${raw_data_merit_N00_S30_W030_E000}'
'${raw_data_merit_N00_S30_E000_E030}'
'${raw_data_merit_N00_S30_E030_E060}'
'${raw_data_merit_N00_S30_E060_E090}'
'${raw_data_merit_N00_S30_E090_E120}'
'${raw_data_merit_N00_S30_E120_E150}'
'${raw_data_merit_N00_S30_E150_E180}'
/
&sgsl_io_extpar
 lpreproc_oro=.TRUE.,
 sgsl_buffer_file='sgsl_buffer.nc',
 sgsl_files = \
'${raw_data_sgsl_N90_N60_W180_W150}'
'${raw_data_sgsl_N90_N60_W150_W120}'
'${raw_data_sgsl_N90_N60_W120_W090}'
'${raw_data_sgsl_N90_N60_W090_W060}'
'${raw_data_sgsl_N90_N60_W060_W030}'
'${raw_data_sgsl_N90_N60_W030_E000}'
'${raw_data_sgsl_N90_N60_E000_E030}'
'${raw_data_sgsl_N90_N60_E030_E060}'
'${raw_data_sgsl_N90_N60_E060_E090}'
'${raw_data_sgsl_N90_N60_E090_E120}'
'${raw_data_sgsl_N90_N60_E120_E150}'
'${raw_data_sgsl_N90_N60_E150_E180}'
'${raw_data_sgsl_N60_N30_W180_W150}'
'${raw_data_sgsl_N60_N30_W150_W120}'
'${raw_data_sgsl_N60_N30_W120_W090}'
'${raw_data_sgsl_N60_N30_W090_W060}'
'${raw_data_sgsl_N60_N30_W060_W030}'
'${raw_data_sgsl_N60_N30_W030_E000}'
'${raw_data_sgsl_N60_N30_E000_E030}'
'${raw_data_sgsl_N60_N30_E030_E060}'
'${raw_data_sgsl_N60_N30_E060_E090}'
'${raw_data_sgsl_N60_N30_E090_E120}'
'${raw_data_sgsl_N60_N30_E120_E150}'
'${raw_data_sgsl_N60_N30_E150_E180}'
'${raw_data_sgsl_N30_N00_W180_W150}'
'${raw_data_sgsl_N30_N00_W150_W120}'
'${raw_data_sgsl_N30_N00_W120_W090}'
'${raw_data_sgsl_N30_N00_W090_W060}'
'${raw_data_sgsl_N30_N00_W060_W030}'
'${raw_data_sgsl_N30_N00_W030_E000}'
'${raw_data_sgsl_N30_N00_E000_E030}'
'${raw_data_sgsl_N30_N00_E030_E060}'
'${raw_data_sgsl_N30_N00_E060_E090}'
'${raw_data_sgsl_N30_N00_E090_E120}'
'${raw_data_sgsl_N30_N00_E120_E150}'
'${raw_data_sgsl_N30_N00_E150_E180}'
'${raw_data_sgsl_N00_S30_W180_W150}'
'${raw_data_sgsl_N00_S30_W150_W120}'
'${raw_data_sgsl_N00_S30_W120_W090}'
'${raw_data_sgsl_N00_S30_W090_W060}'
'${raw_data_sgsl_N00_S30_W060_W030}'
'${raw_data_sgsl_N00_S30_W030_E000}'
'${raw_data_sgsl_N00_S30_E000_E030}'
'${raw_data_sgsl_N00_S30_E030_E060}'
'${raw_data_sgsl_N00_S30_E060_E090}'
'${raw_data_sgsl_N00_S30_E090_E120}'
'${raw_data_sgsl_N00_S30_E120_E150}'
'${raw_data_sgsl_N00_S30_E150_E180}'
/
EOF_oro

cat > INPUT_OROSMOOTH << EOF_orosm
&orography_smoothing
  lfilter_oro=.TRUE.,
  ilow_pass_oro=1,
  numfilt_oro=1,
  ilow_pass_xso=0,
  lxso_first=.FALSE.,
  numfilt_xso=1,
  rxso_mask=0.0,
  eps_filter=0.1,
  rfill_valley=0.0,
  ifill_valley=0
/
EOF_orosm

cat > INPUT_RADTOPO << EOF_radtopo
&radtopo
  lradtopo=.FALSE.,
  nhori=24
/
EOF_radtopo

cat > INPUT_SCALE_SEP << EOF_scale_sep
&scale_separated_raw_data
  lscale_separation = .FALSE.,
  raw_data_scale_sep_path = '',
  scale_sep_files = ''
/
EOF_scale_sep

cat > INPUT_SOIL << EOF_soil
&soil_raw_data
 isoil_data = 3,
 ldeep_soil = .FALSE.,
 raw_data_soil_path='',
 raw_data_soil_filename='${raw_data_soil_HWSD}',
 raw_data_deep_soil_filename='${raw_data_deep_soil}'
/
&soil_io_extpar
  soil_buffer_file='${buffer_soil}',
/
&HWSD_index_files
 path_HWSD_index_files='',
 lookup_table_HWSD='${raw_lookup_table_HWSD}', 
 HWSD_data='${raw_HWSD_data}',
 HWSD_data_deep='${raw_HWSD_data_deep}',
 HWSD_data_extpar='${raw_HWSD_data_extpar}'
/
EOF_soil

cat > INPUT_FLAKE << EOF_flake
&flake_raw_data
   raw_data_flake_path='',
   raw_data_flake_filename='${raw_data_flake}'
/
&flake_io_extpar
   flake_buffer_file='${buffer_flake}'
/
EOF_flake

# Consistency check
cat > INPUT_CHECK << EOF_check
&extpar_consistency_check_io
  grib_output_filename="${grib_output_filename}",
  grib_sample="${grib_sample}",
  netcdf_output_filename="${netcdf_output_filename}",
  i_lsm_data=1,
  land_sea_mask_file="",
  number_special_points=0,
  lwrite_grib=.FALSE.,
/  
EOF_check

# -------------------------------------------------------------------------------
# Launch EXTPAR executables in parallel
# -------------------------------------------------------------------------------
if $run_parallel ; then

    echo ">>>> EXTPAR is run in parallel <<<<"
        
    run_parallel ${binary_alb}
    run_parallel ${binary_aot}
    run_parallel ${binary_tclim}
    run_parallel ${binary_lu}
    run_parallel ${binary_soil} 
    run_parallel ${binary_flake}
    run_parallel ${binary_ndvi} 
    run_parallel ${binary_topo} 

    # Wait for all parallel executables to end
    wait

    # Count non-zero exit status
    error_count=0
    check_exit_status ${binary_alb} error_count
    check_exit_status ${binary_aot} error_count
    check_exit_status ${binary_tclim} error_count
    check_exit_status ${binary_lu} error_count
    check_exit_status ${binary_topo}  error_count
    check_exit_status ${binary_ndvi}  error_count
    check_exit_status ${binary_soil}  error_count
    check_exit_status ${binary_flake} error_count

    # Exit script in case some EXTPAR executables failed
    if [[ $error_count > 0 ]]; then

        echo "*****************************************"
        echo ""
        echo "Some Extpar executables did not terminate correctly!"
        echo "See ${logfile} for more information"
        echo ""
        echo "*****************************************"
        exit 

    fi

    run_sequential ${binary_consistency_check}
    # -> the consistency check requires the output of ${binary_aot},
    #    ${binary_tclim}, ${binary_lu}, ${binary_globe}, ${binary_ndvi},
    #    ${binary_soil} and ${binary_flake}

    # Clean-up
    rm exit_status_*
    rm time_*

# -------------------------------------------------------------------------------
# Launch EXTPAR executables sequentially
# -------------------------------------------------------------------------------
else 

    echo ">>>> EXTPAR is run sequentially <<<<"

    run_sequential ${binary_alb}
    run_sequential ${binary_aot}
    run_sequential ${binary_tclim}
    run_sequential ${binary_lu}
    run_sequential ${binary_soil}
    run_sequential ${binary_flake}
    run_sequential ${binary_ndvi}
    run_sequential ${binary_topo}
    
    # -> error checking is performed within 'run_sequential' functions

    run_sequential ${binary_consistency_check}
    # -> the consistency check requires the output of ${binary_aot},
    #    ${binary_tclim}, ${binary_lu}, ${binary_globe}, ${binary_ndvi},
    #    ${binary_soil} and ${binary_flake}

fi
# -------------------------------------------------------------------------------

echo ">>>> External parameters for COSMO model generated <<<<"
