import pandas as pd
from Bio import SeqIO
import numpy as np
import math
def main():
    config = "/root/comp561/1a_output.txt"
    myfa= list(SeqIO.parse("/root/comp561/Vibrio_vulnificus.ASM74310v1.dna.toplevel.fa", "fasta"))
    print (type(myfa[0].seq))
#     print ("read inputs done!")
#     run(myfa,config)
    
if __name__ == "__main__":
    main()
    
def viterbi(s, config):
# config = "/root/comp561/1a_output.txt"
    wq = float("-inf")
    stateSet = ["inter","start","middle","stop"]
    i = [1,0,0,0,0,0]
    transm = pd.DataFrame(columns = stateSet, index = stateSet)
    #the rows of transm should be the same as state order in column; column->row means transmission from 
    #this column state to the row state
    emis_in = []
    emis_start = []
    emis_m = []
    emis_t = []
    ein = pd.DataFrame()
    es = pd.DataFrame()
    et = pd.DataFrame()
    em = pd.DataFrame()
    #rows in order as specified in 1a_output.txt
    with open(config, "r") as c:
        tmp = []
        tmp2 = []
        lines=c.readlines()
        startindex = 0
        stopindex = 0
        midindex = 0
        for i in range(len(lines)):
            if "start" in lines[i]:
                startindex = i
            elif "stop" in lines[i]:
                stopindex = i
            elif "middle" in lines[i]:
                midindex = i
        for i in range (11,15,1):
            tmp = lines[i].split(' ')
            #this is just for level 1
            if float (tmp[1].strip('\n'))!=0:
                emis_in.append(float (tmp[1].strip('\n')))
            else: 
                emis_in.append(0.000001)
    #             emis_in.append(float (tmp[1].strip('\n')))    
        ein = pd.DataFrame(emis_in, index = ["A","T","G","C"], columns = ["intergenic emission"])
        for i in range (startindex+1,stopindex,1):
            tmp = lines[i].split(' ')
            tmp2.append(tmp[0].strip(':'))
            emis_start.append(float (tmp[1].strip('\n')))
        es = pd.DataFrame(emis_start, index = [tmp2], columns = ["start emission"])
        tmp2 = []
        for i in range (stopindex+1,midindex,1):
            tmp = lines[i].split(' ')
            tmp2.append(tmp[0].strip(':'))
            emis_t.append(float (tmp[1].strip('\n')))
        et = pd.DataFrame(emis_t, index = [tmp2], columns = ["stop emission"])
    #         print (ein)
    #         print (es)
        tmp2 = []
        for i in range (midindex+1,len(lines),1):
            tmp = lines[i].split(' ')
            tmp2.append(tmp[0].strip(':'))
            emis_m.append(float (tmp[1].strip('\n')))
        em = pd.DataFrame(emis_m,index = [tmp2],columns = ["middle emission"])
    #         print (em)
    #         print (et)
        tmp = lines[3].split(' ')
        avg_inter = (float (tmp[5].strip('\n')))
        tmp = lines[8].split(' ')
        avg_gene = (float(tmp[5].strip('\n')))
    # transm = pd.DataFrame(columns = ["inter","start","middle","stop"])
    transm["inter"] = [((avg_inter-1)/avg_inter), (1/avg_inter), 0,0]
    transm["start"] = [0.0,0.0,1.0,0.0]
    transm["middle"] = [0.0,0.0,((avg_gene-1)/avg_gene),(1/avg_gene)]
    transm["stop"] = [1.0,0.0,0.0,0.0]   
    # print (transm)
    # print (es)
    # print (em)
    # print (et)
        #construct the filling dataframe and ptr

#     print ("transition and emission matrices filled")
    #     grand_obj = ""
    #     for i in range (0,len(myfa),1):
    #         grand_obj = grand_obj+myfa[i].seq
    #     s = list(grand_obj)
    #     s = myfa[0].seq
    seqLen= len(s)
    print (seqLen)


    # s = list(myfa[0].seq[0:3000])
    prb= np.full((4,seqLen),float("-inf"))
    ptr = np.full((4,seqLen),-9)
    prb[0,0] = (ein.loc[s[0]])
    #for ptr, from state inter =0, start = 1, mid = 2, stop = 3
    ptr[0,0] = 0
    prb[:,0]
#     print ("big matrices initiated")
    for i in range(1,seqLen,1):
        if (i%50000==0):
            print ("im at column "+str(i))
        #for this is inter
        tmplst = []
        tmplst.append (prb[0,i-1]+math.log(transm.loc["inter","inter"]))#from inter to inter
        tmplst.append (prb[3,i-3]+math.log(transm.loc["inter","stop"]))#from stop to inter
        prb[0,i] = math.log(ein.loc[s[i]].values[0]) + max(tmplst)
        #store the pointer
        if (tmplst[0]>tmplst[1]) and ptr[0,i-1]==0:
            ptr[0,i]=0
        if (tmplst[0]<tmplst[1]) and (i>=3):
            if(ptr[3,i-3]==2)and (max(prb[:,i-3])==prb[3,i-3]):
                #if 3 nucs ago, stop is from middle state: stop state appeared 3 nucs ago, so this is inter
                ptr[0,i]=3

        if (i<=seqLen-3):
            triplets = s[i]+s[i+1]+s[i+2]

            #for this is stop
            #only coming from middle and when it is the first of any of our stop codon
            if (triplets in et.index)and (prb[2,i-3]!=wq):
                prb[3,i] = math.log(et.loc[triplets].values[0])+prb[2,i-3]+math.log(transm.loc["stop","middle"])
                if (ptr[2,i-2]!=1)and(prb[2,i-1]!=1):
                    ptr[3,i] = 2
                
            #start
        #     if current and next 2 combines as start codon, it is either inter->start or mid->mid
        #     to fill prb[1,i], i need inter->start
            elif (triplets in es.index)and (prb[0,i-1]!=wq):
                prb[1,i] = math.log(es.loc[triplets].values[0])+prb[0,i-1]+math.log(transm.loc["start","inter"])
                if max(prb[:,i-1])==prb[0,i-1]:#)and (ptr[2,i-3]!=2):#if previous is not middle
                    ptr[1,i] = 0

            #to fill prb[2,i], i need mid ->mid
            if (triplets in em.index) and (prb[1,i-3]!=wq or prb[2,i-3]!=wq):
                tmplst = []
                #from start to mid
                tmplst.append(prb[1,i-3]+math.log(transm.loc["middle","start"]))
                #from mid to mid
                tmplst.append(prb[2,i-3]+math.log(transm.loc["middle","middle"]))
                prb[2,i] = math.log(em.loc[triplets].values[0])+max(tmplst)
                if (tmplst[0]>tmplst[1])and(ptr[1,i-3]==0):
                    ptr[2,i] = 1
                elif(tmplst[0]<tmplst[1])and(ptr[2,i-3]==2 or ptr[2,i-3]==1):
                    ptr[2,i] = 2
            
#     print ("matrices filled")
    return prb,ptr

def traceback(s,prb,ptr):
    the_max_prob = float('-inf')
    the_index = 0
    the_stops = []
    the_starts = []

    for i in range(4):
    #i wanna first find the state of the last nucleotide
    #compare and find the max prob
        if prb[i,len(s)-1]>the_max_prob:
            the_max_prob = ptr[i,(len(s)-1)]
            the_index = i
#     print ("now i know the state stored in the_index")
    #for loop: note you need to refer different row for different states
    x = len(s)-1
    while x>=0:
        if the_index ==0: #max(prb[:,x])
            #if my current nuc is inter, then i look at [3,x-3] to see if it's 2
            if (ptr[3,x-3]==2) and (ptr[0,x]==3):#update with previous max state
    #     print(the_max_state)
                the_stops.append(x)
                the_index =3#update to the cur max state
                x-=3#you skip 2
            else:
                x-=1
                
        elif the_index ==2:
            # middle: coming from start
            if ptr[1,x-3]==0 and (ptr[2,x-3]==-9):
                the_starts.append(x-2)
                the_index =1
                x=x-3
            #else if current comes from from middle
#             elif (ptr[2,x-3]==2 or ptr[2,x-3]==1):
#                 x=x-3
            else:
                x-=3
                
        elif the_index ==1:
            if (ptr[2,x-3]==-9):
                the_index = 0 #update the state to be inter
                x-=1
            else:
                the_index = 2
                x-=3
        
        else:
            the_index =2
            x=x-3
               
                
    the_starts.reverse()
    the_stops.reverse()
    return the_starts, the_stops



def write_to_output(ind,start,stop):
    for i in range (0,len(start),1):
        with open ("/root/comp561/1c_output.gff3","a+") as f:
            newline = ind +"\t"+"ena"+"\t" +"CDS"+str(start[i]) + "\t"+str(stop[i]) + "\t"+"."+"\t"+"0"+"\t"+"."+"\t"+"\n"
            f.writelines(newline)

def run(m,c):
    for i in range (0,len(m),1):
        prb,ptr = viterbi(m[i].seq,c)
        print ("viterbi for seq " + str(i) +"done")
        start,stop = traceback(list(m[i].seq),prb,ptr)
        print ("traceback done")
        #     print (len(start))
        #     print (len(stop))
        write_to_output(m[i].id,start, stop)
        #     print (start)
        #     print (stop)
        #     print (m[i].seq[(start[9]):(stop[9])])
        #     print (m[i].seq[2240:2256])
        print (m[i].id)
        print ("write finished")
