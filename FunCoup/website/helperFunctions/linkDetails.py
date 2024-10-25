import numpy as np
from sklearn.preprocessing import MinMaxScaler
from data.models import *
from django.db.models import Q
from django.apps import apps
from .colors import *


def getMetadataText(metadata, type):
    # Constructing strings of detailed information. 
    # Information needs to be presented differently depending on evidence type

    textString=""
    if metadata!=None and metadata!="Null" and metadata!="null" and metadata!="":
        if type=="PHP":
            textString="Co-conserved in "+str(metadata.split(",")[0])+" out of "+str(metadata.split(",")[1])+" species"
        elif type=="MEX":
            mexExperiment,score=metadata.split("|")
            textString="Correlation of "+str(score)+" in "+mexExperiment
        elif type=="SCL":
            metadata=metadata.split(",")
            withLinks=""
            for m in metadata:
                withLinks+="<a href='https://www.ebi.ac.uk/QuickGO/term/"+m+"' target='_blank'>"+m+"</a>, "
            len_metadata=len(metadata)
            if len_metadata > 1: textString="Shared locations: "+withLinks[:-2]
            else: textString="Shared location: "+withLinks[:-2]
        elif type=="DOM":
            textString=metadata
        elif type=="PIN":
            metadata=metadata.replace("pubmed:","").split(",")
            withLinks=""
            for m in metadata:
                withLinks+="<a href='https://pubmed.ncbi.nlm.nih.gov/"+m+"' target='_blank'>"+m+"</a>, "
            len_metadata = len(metadata)
            if len_metadata > 1: textString=str(len_metadata)+" publications with pmid: "+withLinks[:-2]
            else: textString=str(len_metadata)+" publication with pmid: "+withLinks[:-2]
        elif type=="GIN":
            exp,score=metadata.split("|")
            textString="Correlation of "+score+" "+exp
        elif type=="GRG":
            exp,score=metadata.split("|")
            withLinks=""
            for m in exp.split(","):
                withLinks+="<a href='https://www.encodeproject.org/files/"+m+"' target='_blank'>"+m+"</a>, "
            textString="Score of "+score+" from "+withLinks[:-2]
        elif type=="MIR":
            metadata = metadata.split(",")
            len_metadata=len(metadata)
            max_num = 100
            if len_metadata>max_num: 
                metadata = ", ".join(metadata[0:int(np.min([len_metadata,max_num/2]))])
                textString="Co-regulated by: "+metadata+", and "+str(int(len_metadata-max_num/2))+" more"
            else:
                metadata=", ".join(metadata)
                textString="Co-regulated by: "+metadata
        elif type=="PEX":
            meta,score=metadata.split("|")
            if meta!="":
                meta=meta.split(",")
                metadata=", ".join(meta)
                len_meta=len(meta)
                if len_meta > 1: textString="Jaccard index of: "+str(score)+" with common terms: "+metadata
                else: textString="Jaccard index of: "+str(score)+" with common term: "+metadata
            else:
                textString="Correlation of: "+str(score)
        elif type=="TFB":
            metadata = metadata.split(",")
            len_metadata=len(metadata)
            max_num = 100
            if len_metadata>max_num: 
                metadata = ", ".join(metadata[0:int(np.min([len_metadata,max_num/2]))])
                textString="Co-regulated by: "+metadata+", and "+str(int(len_metadata-max_num/2))+" more"
            else:
                metadata=", ".join(metadata)
                textString="Co-regulated by: "+metadata

    else:
        textString="No data available"

    return textString

def getDetailsForEvidence(evidenceType,evidenceDetails,species,query_species,goldStandardInDb,orthologMap,orthologTransfer,host_pathogen_flag=False):
    # Getting detailed information for each evidence type
    allEvidence=[]

    for sp in orthologMap: # Iterating over all species having orthologs
        
        if evidenceType in orthologTransfer or sp==query_species: # Only doing transfer for evidence that should be transferred or for the species itself
            # Getting evidenceID to get info for
            evidenceInDb=list(Evidence.objects.filter(Q(type=evidenceType) & Q(species__tax_id=sp)).values_list("id", flat=True))
            orthoA=orthologMap[sp]['A']
            orthoB=orthologMap[sp]['B']
            for evInDb in evidenceInDb: # Iterating over evidences (since there are multiple per type for e.g. mex)
                for a in orthoA: # Iterating over ortholog ids in case there are one-to-many mappings
                    for b in orthoB:
                        evidenceModel = apps.get_model(app_label='data', model_name=evidenceType) # Getting the name of the table to query             
                        if int(a)>int(b): # Checking in which order to query protein ids
                            allEvidence.extend(evidenceModel.objects.filter(Q(proteinA_id=a) & Q(proteinB_id=b) & Q(evidence_id=evInDb)))
                        else:
                            allEvidence.extend(evidenceModel.objects.filter(Q(proteinA_id=b) & Q(proteinB_id=a) & Q(evidence_id=evInDb)))

    # Iterating over all evidence collected
    for ev in allEvidence:
        # Collecting the LLR function
        if host_pathogen_flag:
            ev_species = Evidence.objects.get(Q(type=evidenceType) & Q(species__tax_id=species))
            pLLR = PolynomialLLR.objects.get(Q(goldStandard=goldStandardInDb) & Q(evidence=ev_species))
        else:
            pLLR = PolynomialLLR.objects.get(Q(goldStandard=goldStandardInDb) & Q(evidence=ev.evidence))
        coefficients = pLLR.function
        if 'Nan' not in coefficients:
            spName=ev.evidence.species.species_name
            min_score, max_score = pLLR.scoreRange.split('|')
            coefficients = [float(i) for i in coefficients.split(',')]
            metadata=''
            if ev.metadata!=None:
                metadata=ev.metadata

            if evidenceType=='MEX':
                ### Convert to raw values and set cutoff for MEX
                scoreRange = ev.evidence.scoreRange.split('|')
                if ev.evidence.scoreRange!=None:
                    original_mex_min,original_mex_max = float(scoreRange[0]),float(scoreRange[1])
                    originalScore = np.round(ev.score * (original_mex_max - original_mex_min) + original_mex_min,3)
                    if originalScore>0.5:
                        ## fake a vector with values between 0.5 and original_mex_max
                        mex_range = np.array([0.5,originalScore,original_mex_max])
                        scaler = MinMaxScaler(feature_range=(0, 1))
                        mex_range = np.round(scaler.fit_transform(mex_range.reshape(-1, 1)), decimals=4)
                        llrValue = np.round(np.polyval(coefficients, np.clip(mex_range[1][0], float(min_score), float(max_score))), 3)
                    else: 
                        # Skip this evidence, it was not used
                        continue

                    # For MEX, the sign need to be put back to the original score as well
                    sign=""
                    if ev.sign!=None:
                        sign=ev.sign
                    originalScore=sign+str(originalScore)
                    metadata=metadata+"|"+str(originalScore)
            else:
                if evidenceType=="GRG" or evidenceType=="GIN" or evidenceType=="PEX":
                    # For some evidence types, scores are scaled and normalized, this needs to be reverted to get a proper score 
                    llrValue = np.round(np.polyval(coefficients, np.clip(ev.score, float(min_score), float(max_score))), 3)
                    if ev.evidence.scoreRange!=None:
                        original_min,original_max=ev.evidence.scoreRange.split("|")
                        originalScore=np.round(ev.score * (float(original_max) - float(original_min)) + float(original_min),3)
                        metadata=metadata+"|"+str(originalScore)
                    else:
                        metadata=metadata+"|"+str(ev.score)
                else:
                    if host_pathogen_flag:
                        # Min-Max scaling to the min,max of polyLLR in host species
                        originalScore = np.array([0,ev.score,1])
                        scaler = MinMaxScaler(feature_range=(float(min_score), float(max_score)))
                        minmax_score = scaler.fit_transform(originalScore.reshape(-1, 1))
                        llrValue = np.round(np.polyval(coefficients, np.clip(minmax_score[1][0], float(min_score), float(max_score))), 3)
                    else:
                        # Getting the LLR (from score, according to LLR function for this evidence)
                        llrValue = np.round(np.polyval(coefficients, np.clip(ev.score, float(min_score), float(max_score))), 3)

            ## TODO --> During network training, in getNetwork.orthologEvidence function, 
            ## we solve many-to-one relationships by averaging the raw score values. 
            ## This is not implemented in the website (we list all values), meaning there 
            ## can be multiple raw score for the same link from the same data source

            # Collecting info on raw score, LLR and metadata, to put int he table for each evidence in support of this link
            evidenceObject={'llr':llrValue,'llrColor':getLLRColor(llrValue),'type':evidenceType,'info':getMetadataText(metadata, evidenceType),'species':spName[0]+(spName.split()[1][:2]).upper()}
            evidenceDetails.append(evidenceObject)

    return evidenceDetails
