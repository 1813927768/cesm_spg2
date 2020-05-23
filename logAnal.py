import os,re,json
import numpy as np

def get_LS_num(fileName):
    countList = [0 for i in range(50)]
    with open("log/"+fileName,'r') as fr:
        for line in fr.readlines():
            matchObj = re.search(r"func start \((\d+)\)",line)
            if not matchObj is None:
                tmp = matchObj.group(1)
                countList[int(tmp)] += 1
    return countList

def LS_analy(fileName):
    '''
    analyze details of line search 
    '''
    switch = False
    count = 0
    alpha = 1
    with open("log/"+fileName,'r') as fr:
        for line in fr.readlines():
            matchObj = re.search(r"line search (start|done) \((\d+)",line)
            if not matchObj is None:
                LSnum = int(matchObj.group(2))
                switch = not switch
                if switch:
                    count = 0
                else:
                    # print("iter %d: line search nums %d"%(LSnum,count))
                    pass
            if switch:
                matchObj = re.search(r"old f = (-\d+\.\d*)",line)
                if matchObj:
                    oldf = float(matchObj.group(1))
                    print("old f = %f"%(oldf))
                
                matchObj = re.search(r"^d = (.*)",line.strip())
                if matchObj:
                    d = json.loads(matchObj.group(1))
                    dmax = np.linalg.norm(np.array(d),np.inf)
                    print("dmax = %f"%(dmax))
                
                matchObj = re.search(r"alpha = (.*)",line.strip())
                if matchObj:
                    alpha = json.loads(matchObj.group(1))
                
                matchObj = re.search(r"value = (-\d+\.\d*)",line)
                if matchObj:
                    count += 1
                    val = float(matchObj.group(1))
                    print("dmax %.15f: %.10f"%(alpha*dmax,val-oldf))


def change_Anal():
    change = np.zeros((50,))
    with open("log/change.out",'r') as fr:
        for line in fr.readlines():
            parts = line.split("=")
            body = json.loads(parts[1])
            change += np.array(body)
    print(list(change))



if __name__ == "__main__":
    get_LS_num("process_0.01_42.txt")
    # LS_analy("process_0.01_18.txt")
    # change_Anal()