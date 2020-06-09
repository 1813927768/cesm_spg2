# create a case
create_newcase -mach Tianhe2 -case $1 -compset F -res f09_g16

cd $1

# modify env_mach_pes.xml
./xmlchange -file ./env_mach_pes.xml -id NTASKS_ATM  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_LND  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_ICE  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_OCN  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_CPL  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_GLC  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_ROF  -val 48
./xmlchange -file ./env_mach_pes.xml -id NTASKS_WAV  -val 48

# modify env_run.xml
./xmlchange -file env_run.xml -id GET_REFCASE  -val FALSE
./xmlchange -file env_run.xml -id RUN_TYPE  -val branch
./xmlchange -file env_run.xml -id RUN_STARTDATE  -val 0001-01-01
./xmlchange -file env_run.xml -id START_TOD  -val 0
./xmlchange -file env_run.xml -id RUN_REFCASE  -val mycase
./xmlchange -file env_run.xml -id RUN_REFDATE  -val 0052-11-15
./xmlchange -file env_run.xml -id RUN_REFTOD  -val 00000
./xmlchange -file env_run.xml -id BRNCH_RETAIN_CASENAME  -val TRUE
./xmlchange -file env_run.xml -id STOP_OPTION  -val ndays
./xmlchange -file env_run.xml -id STOP_N  -val 15

# cesm_setup
./cesm_setup -clean
./cesm_setup

# prepare input_file
cp -r /BIGDATA1/iocas_mmu_3/support/input_15/. ./run

# edit name_list and run file
echo "nhtfrq = -24" >> user_nl_cam
sed -i '/cesm\.exe/a\mkdir \$CASEROOT\/finished' mycase.run

# build
./mycase.clean_build
./mycase.build

# run (submit to the job system)
yhbatch -n 48 mycase.run
