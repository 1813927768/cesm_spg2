## Algorithm
 
1. set init parameters: delta, x0
2. calculate func(x0), gradient(x0)
3. While the stopping criterion is not satisfied 

    3.1 xnew = lineSearch(x,g)
    
    3.2 gnew = gradient(xnew)


## Operation Steps

> Updates: Now you can run ./setup_xml.sh ./1/mycase to complete the whole process before step 4 running

1. create a case
```
create_newcase -mach Tianhe2 -case /BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/3/mycase -compset F -res f09_g16
```

2. setup the case

    2.1 setup xml files
    ```
    ./setup_xml.sh ./pre/mycase
    ```

    2.2 edit user_nl_cam
    ```
    nhtfrq = -24
    ```

    2.3 edit mycase.run
    
    ```
    ## yhrun -N ${SLURM_NNODES} -n ${SLURM_NTASKS}  ${EXEROOT}/cesm.exe 
    # cd $CASEROOT
    # mkdir finished
    ## echo "`date` -- CSM EXECUTION HAS FINISHED" 
    sed -i '/cesm\.exe/a\mkdir \$CASEROOT\/finished' mycase.run
    ```

3. build the case 

```
./mycase.clean_build
./mycase.build
```

4. run 
```
yhbatch -n 48 ./mycase.run
```