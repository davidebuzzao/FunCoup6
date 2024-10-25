import numpy as np
import matplotlib as mpl

def getLLRColor(llr):
    bigLLR=5.0
    smallLLR=-1.5
    if llr >bigLLR:
        llr=bigLLR
    elif llr < smallLLR:
        llr=smallLLR
    
    if llr>0:
        val=(bigLLR - llr) / bigLLR
        c1='#45b584' # Put max saturated color for high LLRs here
        c2='#ffffff' # ranging to white 
        c1=np.array(mpl.colors.to_rgb(c1))
        c2=np.array(mpl.colors.to_rgb(c2))
        return mpl.colors.to_hex((1-val)*c1 + val*c2)
    elif llr < 0:
        val=(smallLLR - llr) / smallLLR
        c1='#f05d5e' # Put max saturated color for low LLRs here
        c2='#ffffff' # ranging to white   
        c1=np.array(mpl.colors.to_rgb(c1))
        c2=np.array(mpl.colors.to_rgb(c2))
        return mpl.colors.to_hex((1-val)*c1 + val*c2)
    elif llr==0:
        return '#ffffff'
    
def getNodeColors(): # List of colors is used to color nodes for different species for comparative interactomics
    colors = ["#4584b5", "#4ca64e", "#6f4483", "#ffb22c", "#f7909f", "#fe672a", "#83ba69", "#8b8784","#000000","#e2647f", "#d1698c", "#da332a", "#ff7a34", "#4ca64e", "#2a897e", "#62b0dc", "#32436f", "#5f5395", "#5a3340", "#9e6c34", "#4b7f52", "7293a0", "#e8555f"]
    return colors 

def getQueryBorder(): # Set color for query node border here
    return "#393939"