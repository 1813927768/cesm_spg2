# useful bash gists

# kill all running tasks
yhcancel `yhq | grep -Eo [0-9]{7}`

# find unsuccessful case by log
 grep -c "CSM EXECUTION HAS FINISHED" ./*.out | grep "0$"
 
# get linesearch nums for iter:1-12
 for i in $(seq 1 12); do echo $(egrep "func start \($i\)" process.txt | wc -l); done

# count words in tex
 texcount contents/Chapter3.tex contents/Chapter2.tex contents/Chapter1.tex main.tex