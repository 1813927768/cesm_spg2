## Operation Steps

1. create a case
```
create_newcase -mach Tianhe2 -case /BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/1/mycase -compset F -res f09_g16
```

2. setup the case

    2.1 setup xml files
    ```
    ./setup_xml.sh ./1/mycase
    ```

    2.2 edit user_nl_cam
    ```
    nhtfrq = -24
    ```

    2.3 run cesm_setup
    ```
    ./cesm_setup -clean
    ./cesm_setup
    ```

    2.4 edit mycase.run
    ```
    ## yhrun -N ${SLURM_NNODES} -n ${SLURM_NTASKS}  ${EXEROOT}/cesm.exe 
    cd $CASEROOT
    mkdir finished
    ## echo "`date` -- CSM EXECUTION HAS FINISHED" 
    ```

3. build the case 

```
./mycase.clean_build
./mycase.bulid
```

4. run 
```
yhbatch -n 48 ./mycase.run
```