import numpy as np

def logistic(x, a, b, c): # Logistic function for getting FBS
    return a / (1 + np.exp(-b * (x - c)))

def FBStoPPV(FBS,logPPV): # Translating FBS to PFC using formula previously picked up from the database
    coefficients = logPPV
    if 'Nan' not in coefficients:
        coefficients = [float(i) for i in coefficients.split(',')]
        a, b, c = coefficients
        return np.round(np.clip(logistic(np.clip(FBS, 0, np.max(FBS)),a, b, c),0,1),3)
    return 0
    
def getNewFbs(row,indexesToConsider,ppvFunctions,pfc_constant=0.001): 
    # Getting all llr information and excluding evidences according to query. 
    # indexesToConsider includes all evidence to be considered in the calculations
    newFBS=""
    newLLR_evidence=""
    newPPV=""
    newMaxPPV=0
    for gsIndex,llrPerGs in enumerate(row['llr_evidence'].split(";")): # Iterating over gold standards
        evLLrs=llrPerGs.split("|")[0].split(",")
        fbs=0
        for i,e in enumerate(evLLrs): # Iterating over evidence LLRs
            if i in indexesToConsider:
                fbs+=float(e)
                newLLR_evidence+=str(e)+","
            else:
                newLLR_evidence+="0," # Setting excluded evidences LLR to 0
        newLLR_evidence=newLLR_evidence[:-1]+"|"+llrPerGs.split("|")[1]+";" # Appending all species information untouched
        newFBS=newFBS+str(round(fbs,3))+","
        # pfc=round(1/(1+np.exp(-np.log(pfc_constant)-fbs)),3)
        if ppvFunctions[gsIndex]!="NONE": # If a function could be extracted from the database
            ppv = FBStoPPV(fbs,ppvFunctions[gsIndex])
        else:
            ppv=0 # If no function, ppv i set to 0
        newPPV+=str(ppv)+","
        if ppv>newMaxPPV: # Setting the maximum of the new ppv
            newMaxPPV=ppv
    row["fbs_goldstandard"]=newFBS[:-1] # just removing last comma
    row["llr_evidence"]=newLLR_evidence[:-1]
    row["ppv"]=newPPV[:-1]
    row["max_ppv"]=newMaxPPV
    return row # Returning row with updated values

def getNewFbsSpecies(row,indexesToConsider, ppvFunctions, pfc_constant=0.001):
    # Getting all llr information and excluding species according to query. 
    # indexesToConsider includes all evidence to be considered in the calculations
    newFBS=""
    newLLR_evidence=""
    newPPV=""
    newMaxPPV=0
    for gsIndex,llrPerGs in enumerate(row['llr_evidence'].split(";")): # Iterating over gold standards
        spLLrs=llrPerGs.split("|")[1].split(",")
        newLLR_evidence+=llrPerGs.split("|")[0]+"|" # Appending evidence information untouched
        fbs=0
        for i,e in enumerate(spLLrs): # Iterate over species
            if i in indexesToConsider:
                fbs+=float(e)
                newLLR_evidence+=str(e)+"," 
            else:
                newLLR_evidence+="0," # Setting excluded species LLR to 0
        newLLR_evidence=newLLR_evidence[:-1]+";"
        newFBS=newFBS+str(round(fbs,3))+","
        # pfc=round(1/(1+np.exp(-np.log(pfc_constant)-fbs)),3)
        if ppvFunctions[gsIndex]!="NONE": # If a function could be extracted from the database
            ppv = FBStoPPV(fbs,ppvFunctions[gsIndex])
        else:
            ppv=0 # If no function, ppv i set to 0
        newPPV+=str(ppv)+","
        if ppv>newMaxPPV: # Setting the maximum of the new ppv
            newMaxPPV=ppv
    row["fbs_goldstandard"]=newFBS[:-1]
    row["llr_evidence"]=newLLR_evidence[:-1]
    row["ppv"]=newPPV[:-1]
    row["max_ppv"]=newMaxPPV
    return row # Returning row with updated values

# def getDirection(GRG_dir, GS_dir, grgLLR, globalDirection):
#     dir=0
#     if GRG_dir>0 and grgLLR>0: # SET A REASONABLE CUTOFF FOR WHEN TO SHOW AN ARROW HERE!
#         dir=GRG_dir
#     if GS_dir>1 and GS_dir!=dir and dir!=5:
#         dir+=GS_dir
#     if globalDirection==0:
#         globalDirection=dir
#     elif dir!=5 and globalDirection!=dir:
#         globalDirection+=dir
#     return globalDirection, dir

def getDirection(GRG_dir, GS_dir, grgLLR, grgLLRThreshold, ppv, ppvThreshold, globalDirection):
    localDirection = 0
    if GS_dir == 5:
        globalDirection,localDirection = 5,5
    else:
        if GS_dir in [2,3]:
            localDirection = GS_dir
        if ppv>=ppvThreshold:
            if GRG_dir != localDirection and grgLLR > grgLLRThreshold: # SET A REASONABLE CUTOFF FOR WHEN TO SHOW AN ARROW HERE!
                localDirection += GRG_dir
            if globalDirection == 5: next
            elif localDirection == 5:
                globalDirection = 5
            elif localDirection in [2,3] and globalDirection in [0,2,3] and globalDirection!=localDirection:
                globalDirection+=localDirection

    return globalDirection,localDirection

