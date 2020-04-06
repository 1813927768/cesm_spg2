# $1: the relative path of root case directory
cd $1

# setup env_mach_pes.xml
./xmlchange -file ./env_mach_pes.xml -id NTASKS_ATM  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_LND  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_ICE  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_OCN  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_CPL  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_GLC  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_ROF  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_WAV  -val 48

# setup env_run.xml
./xmlchange -file env_run.xml -id GET_REFCASE  -val FALSE
./xmlchange -file env_run.xml -id RUN_TYPE  -val branch
./xmlchange -file env_run.xml -id RUN_STARTDATE  -val 0001-01-01
./xmlchange -file env_run.xml -id START_TOD  -val 0
./xmlchange -file env_run.xml -id RUN_REFCASE  -val mycase
./xmlchange -file env_run.xml -id RUN_REFDATE  -val 0052-11-01
./xmlchange -file env_run.xml -id RUN_REFTOD  -val 00000
./xmlchange -file env_run.xml -id BRNCH_RETAIN_CASENAME  -val TRUE
./xmlchange -file env_run.xml -id STOP_OPTION  -val ndays
./xmlchange -file env_run.xml -id STOP_N  -val 15

# call cesm_setup
./cesm_setup -clean
./cesm_setup

# prepare input data
cp -r /BIGDATA1/iocas_mmu_3/support/input/. ./run
