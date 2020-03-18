#!/usr/bin/env python
# coding: utf-8

# In[2]:


#blast
#with probabilistic model
# all the import goes here
# from Bio import SeqIO
from timeit import default_timer as timer
from random import random
import pandas as pd
import numpy as np
import math
import six
import string


# In[3]:


#function splitting my query
def splitq (query,wordlength):
    wordlst = []
    if wordlength<=0:
        print ("received negative number for wordlength, returning whole query")
        print (query)
        return list(query)
    elif wordlength > len(query):
        print ("invalid word length against query: returning whole query")
        print (query)
        return list(query)
    else:
        for i in range(0,len(query)-wordlength+1,1):
            if wordlst ==[]:
                wordlst = [str(query[i:i+wordlength])]
            else:
                wordlst.append(str(query[i:i+wordlength]))
#         print (wordlst)
        return wordlst


# In[4]:


#aligning the words on sequences
def alignment(query,wordlength,prob_genome,lstA,lstT,lstC,lstG):
    gl = np.size(prob_genome,1)
#   receives the word list as well as probabilistic genome as input
    wordlst = splitq(query,wordlength)
    l = len(query)
#     chosenList = list_choice(wordlst[0][0],lstA,lstT,lstC,lstG)
    firstList = [x for x in range (0,gl)]
    score,pos = locate(wordlst[0],prob_genome,firstList)
    print ("score is" +str(score))
    print ("pos is" +str(pos))
    p = []
    s = []
    thr = (max(score)-min(score))*3/10
    threshold = -(math.log(0.9)*(wordlength-1))+1
    if len(score)>1:
        for x in score:
#             if x<=(min(score)+thr):
            if x<=threshold:
                if p ==[]:
                    p = [pos[score.index(x)]]
                    s = [x]
                else:
                    if x in s:
                        continue
                    else:
                        p.append(pos[score.index(x)])
                        s.append(x)
    
        for i in range(1,len(wordlst),1):
            chosenList = list_choice(wordlst[i][0],lstA,lstT,lstC,lstG)
            scr2,pos2 = locate(wordlst[i],prob_genome,chosenList)
            for j in range(0,len(pos2),1):
                a=pos2[j]
                if i==1:
                    if (a-i) in pos:
                        if scr2[j]<threshold:
                            threshold = scr2[j]
                        if a-i in p:
                            continue
                        else:
                            p.append(a-i)
                            s.append(score[p.index(a-i)]-scr2[j])
                    elif scr2[j]<=threshold:
                        if a-i in p:
                            continue
                        else:
                            p.append(a-i)
                            s.append(scr2[j])
                    else:
                        continue
                    if len(s)>1:
                        for x in s:
                            if x>=0.15:
                                p.remove(p[s.index(x)])
                                s.remove(x)
                else:
                    if (a-i) in p:
                        if scr2[j]<threshold:
                            threshold = scr2[j]
                        s[p.index(a-i)] = s[p.index(a-i)]-scr2[j]
                    elif scr2[j]<=threshold:
                        if a-i in p:
                            continue
                        else:
                            p.append(a-i)
                            s.append(scr2[j])
                    else:
                        continue
    #                 s,p = remove(s,p)
                if len(s)>1:
                    for x in s:
                        if x>=0.15:
                            p.remove(p[s.index(x)])
                            s.remove(x)
    else: 
        s = score
        p = pos
    if p ==[] or s ==[] or len(p)!=len(s):
        print ("i did not find a potential alignment")
        return -1
    else:
        print len(s)
        print "s is "+str(s)
        print "p is "+str(p)
        if len(p)==1:
            startpoint = p[0]
#             return startpoint
            return s,p
        else:
            if len(p)<=10:
                return s,p
            else:
#                 r = len(s)-1
#                 s,p = quicksort(s,p,0,r)
                s_filtered = s[:10]
                p_filtered = p[:10]
                return s_filtered,p_filtered

# left with an updated list
# if this list length==1, then return start point, end point, and the alignment (query, sequence) 
# else: calculate score under each pos and find the minimum


#   trivial function for finding right place
def locate (word,prob_genome,aList):
    score = []
    pos = []
    gl = np.size(prob_genome,1)
    x = -(math.log(0.7)*len(word))
    for i in range(0,len(aList),1):
        c1 = cvt(word[0])
        #traverse across genome to identify suitable pos for first instance of word
        sati = -(math.log(prob_genome[c1,i]))
        for j in range (1,len(word),1):
            c = cvt(word[j])
            if i+j<gl:
                sati -= math.log(prob_genome[c,(i+j)])
            else:
                i = gl      
        if sati <=x and i!=gl:
            if score ==[]:
                score = [sati]
                pos = [i]
            else:
                score.append(sati)
                pos.append(i)
    return score,pos


# In[5]:


#a function that reads in the genome
def create_db (fasta,conf):
#     function receving file names
#     fasta = "/root/comp561/chr22.maf.ancestors.42000000.complete.boreo.fa"
#     conf = "/root/comp561/chr22.maf.ancestors.42000000.complete.boreo.conf"
    #     fasta = "/root/comp561/200trial.fa"
    #     conf = "/root/comp561/200trial.conf"
    f = open(fasta, "r")
    genome = f.readline()
    if fasta !="/root/comp561/chr22.maf.ancestors.42000000.complete.boreo.fa":
        genome = genome[:-1]
    lstA,lstT,lstC,lstG = idx_filtering(genome)
    c = open(conf,"r")
    lst = c.readline().split(' ')
    lst.pop()
    lst = map(float,lst)
    gl = len(genome)
    print gl
    print len(lst)
    # prob_genome = pd.DataFrame(columns = list(range(0,3000)), index = ['A','T','C','G']) 
    prob_genome = np.zeros((4,gl))
    print prob_genome.shape
    start = timer()
    if len(lst)==gl:
        for i in range(0,gl,1):
            if i%100000==0:
                print ("generated prob_genome column " + str(i))
            ind = genome[i]
            if genome[i]=='A':
                prob_genome[0,i] = lst[i]
                if lst[i]!=1:
                    prob_genome[1:,i] = (1-lst[i])/3
                else:
                    prob_genome [1:,i] = 0.00001
            elif genome[i]=='T':
                prob_genome[0,i] = lst[i]
                if lst[i]!=1:
                    prob_genome[1:,i] = (1-lst[i])/3
                else:
                    prob_genome [1:,i] = 0.00001
                prob_genome[0,i],prob_genome[1,i] =prob_genome[1,i], prob_genome[0,i]
            elif genome[i]=='C':
                prob_genome[0,i] = lst[i]
                if lst[i]!=1:
                    prob_genome[1:,i] = (1-lst[i])/3
                else:
                    prob_genome [1:,i] = 0.00001
                prob_genome[0,i],prob_genome[2,i] =prob_genome[2,i], prob_genome[0,i]
            else:
                prob_genome[0,i] = lst[i]
                if lst[i]!=1:
                    prob_genome[1:,i] = (1-lst[i])/3
                else:
                    prob_genome [1:,i] = 0.00001
                prob_genome[0,i],prob_genome[3,i] =prob_genome[3,i], prob_genome[0,i]
        return prob_genome,lstA,lstT,lstC,lstG
    elif len(lst)==4*gl:
        #implement dataframe
        for i in range(0,gl*4,4):
            ind = genome[i]
        #     print (ind)
            for j in range(0,4,1):
                prob_genome[ind,i+j] = lst[i+j]
        return prob_genome,lstA,lstT,lstC,lstG
    else:
        print ("something wrong about the input")
    end = timer()
    print("time: {}".format(end - start))


# In[6]:


#converter function for nucleotides to index
def cvt(char):
    if char is 'A':
        return 0
    elif char is 'T':
        return 1
    elif char is 'C':
        return 2
    elif char is 'G':
        return 3
    else:
        print ("there exist some non-nucleotidal entries in the query")
        return -1
    


# In[7]:


#a list filtering function being called in create_db
#then you pass all these outputs into alignment function
def idx_filtering(genome):
    lstA = []
    lstT = []
    lstC = []
    lstG = []
    for x in genome:
        if x is 'A':
            if lstA ==[]:
                lstA = [genome.index(x)]
            else:
                lstA.append(genome.index(x))
        elif x is 'T':
            if lstT ==[]:
                lstT = [genome.index(x)]
            else:
                lstT.append(genome.index(x))
        elif x is 'C':
            if lstC ==[]:
                lstC = [genome.index(x)]
            else:
                lstC.append(genome.index(x))
        elif x is 'G':
            if lstG ==[]:
                lstG = [genome.index(x)]
            else:
                lstG.append(genome.index(x))
        else:
            print ("genome containing none nucleotide entries")
    return lstA,lstT,lstC,lstG
#and a function to identify which list to pass
def list_choice(char,lstA,lstT,lstC,lstG):
    if char is 'A':
        return lstA
    elif char is 'T':
        return lstT
    elif char is 'C':
        return lstC
    elif char is 'G':
        return lstG
    else:
        print ("query contains non-nucleotidal entries: program ending")
        return []


# In[8]:


#function to filter list
def remove_excess (s,p,wordlength):
#     thre = 
    for i in s:
        if i >= thre:
            p.remove(p[s.index(i)])
            s.remove(i)
    return s,p


# In[9]:


def quicksort(s,p,l,r):
#Base case: No need to sort says of length <= 1
    if l>=r or r-l==1:
        return s,p
    # Choose pivot to be the last element in the subsay
    else:
        pivot = s[r]
        cnt = l
        big = r
        for i in range(l,r+1,1):
            # If an element less than or equal to the pivot is found...
            if (s[i] <= pivot):
                tmp = s[cnt]
                s[cnt] = s[i]
                s[i] = tmp
                tmp = p[cnt]
                p[cnt] = p[i]
                p[i] = tmp
                cnt+=1
            else:
                continue
        quicksort(s,p,l, cnt-2) 
        quicksort(s,p, cnt, r)  


# In[10]:


# # tester for quicksort
# def main():
#     s = [10, 7, 8, 9, 1, 5,11,6];
#     p = [7,4,5,6,1,2,8,3]
    
#     quicksort(s,p,0, (len(s)-1));
#     print s
#     print p
# if __name__ =="__main__":
#     main()


# In[75]:


# implement Needleman-Wunch# 4 inputs: genome file; score for match; score forr mismatch; gap penalty bdef     
def SuborMatch(ms,mms,s,y):
    charIndex = cvt(y)
#     print charIndex
    if charIndex == np.where(max(s)):       
        return ms #match score  
    
    else:
        for x in range(0,4,1):
            if x==charIndex:
#                 print x
#                 print s[x]
                return math.log(s[x])
            
def mgf(genome,p,query):
    ms = 1
    b = -1
    lenToCompare = int(len(query)*1.5)
    align_lst = []
    num_mut = []
    T = query
    for i in p:
        mut = 0
        S = genome[:,i:i+lenToCompare]
        m = len(S)
        n = len(T)
#         X= np.full((m+1,n+1),float("-inf"))
#         I = np.full((m+1,n+1),float("-inf"))
#         D = np.full((m+1,n+1),float("-inf"))
        X= np.full((m+1,n+1),0.0)
        I = np.full((m+1,n+1),0.0)
        D = np.full((m+1,n+1),0.0)
        align = np.full((m+1,n+1),'n')
        X[0,0]=0
        X[0,1]=b
        X[1,0]=b
        I[0,1]=b
        I[1,0]=b
        D[0,1]=b
        D[1,0]=b
        I[0,0]=0
        D[0,0]=0
        for i in range(1,m+1,1):
            for j in range(1,n+1,1):
                I[i,j] = max(X[i-1,j]+b,D[i-1,j]+b)#b
                D[i,j] = max(X[i,j-1]+b,I[i,j-1]+b)#b
                s=SuborMatch(1,-1,S[:,i-1],T[j-1])
#                 print (s)
#                 print ("for X: " +str(X[i-1,j-1]+s))
#                 print ("for I: "+str(I[i-1,j-1]+s))
#                 print ("for D: " +str(D[i-1,j-1]+s))
                X[i,j] = max(X[i-1,j-1]+s,I[i-1,j-1]+s, D[i-1,j-1]+s)
                if (X[i,j]==X[i-1,j-1]+s):
                    align[i,j] = 'x'
                elif (X[i,j]==I[i-1,j-1]+s):
                    align[i,j] = 'd'
                else:
                    align[i,j] = 'i'
        s_inAlign = []
        t_inAlign = []
        i =m
        j =n
        print X
        print I
        print D
        print ("the final X: " +str(X[i,j]))
        print ("the final I: "+str(I[i,j]))
        print ("the final D: " +str(D[i,j]))
        score = max(X[i,j],I[i,j],D[i,j])
        while (i >=0 and j>=0):
            if align[i,j]=='x':
                continue
#         #         print ("yes")
#                 s_inAlign.append(S[i-1])
#                 t_inAlign.append(T[j-1])
#                 i=i-1
#                 j=j-1
#             elif align[i,j]=='i':
            else:
                mut
#                 t_inAlign.append(T[j-1])
#                 s_inAlign.append('-')
#                 j=j-1
#             elif align[i,j]=='d':
#                 s_inAlign.append(S[i-1])
#                 t_inAlign.append('-')
#                 i=i-1
#             else:
#                 break
#         s_inAlign.reverse()
#         t_inAlign.reverse()
#         align_s = string.join(s_inAlign)
#         align_t = string.join(t_inAlign)
#     with open ('/root/comp401/hw1_feedback_q3.txt','a') as f:
#         f.write("score for short run is"+str(score)+'\n')
#         f.write(align_s+'\n')
#         f.write(align_t+'\n')
#         print (align_s)
#         print (align_t)
        if align_lst ==[]:
            align_lst = [score]
            num_mut = [mut]
        else:
            align_lst.append(score)
            num_mut.append(mut)
    print "aligning scores are" +str(align_lst)
    print "equivalent number of mutations: " +str(num_mut)
    print "their starting pos are "+str(p)


# In[76]:


def main():
    fasta = "/root/comp561/100000trial.fa"
    conf = "/root/comp561/100000trial.conf"
    pg,lstA,lstT,lstC,lstG = create_db(fasta,conf)
    s,p = alignment("CAACTAAC",3,pg[:,:100],lstA[:100],lstT[:100],lstC[:100],lstG[:100])
    mgf(pg[:,:100],p,"CAACTAAC")
main()


# In[74]:


# main function
def main():
#     fake_prob_g = fakegenome()
#     fasta = "/root/comp561/chr22.maf.ancestors.42000000.complete.boreo.fa"
#     conf = "/root/comp561/chr22.maf.ancestors.42000000.complete.boreo.conf"
    fasta = "/root/comp561/100000trial.fa"
    conf = "/root/comp561/100000trial.conf"
    pg,lstA,lstT,lstC,lstG = create_db(fasta,conf)
#     print pg[:200,:200]
    start = timer()
#     s,p = alignment("CAACTAACCACCACCCCTGTCTCCACTCACCGGAACAGAGACTCCCCCAG",7,pg,lstA,lstT,lstC,lstG)#orignal
#     w,s = alignment("CACCTAACCACCTCCCCTGTCTGCACTCACCGGAACAGAGACTCCACCAG",7,pg,lstA,lstT,lstC,lstG)#with 4 mutations
#     w,s = alignment("AAACTAAACACCACTCCTGTCTGCACTCACCGGAAAAGAGACTCCCCCAG",7,pg,lstA,lstT,lstC,lstG)#with mutation at the first spot
    s,p = alignment("CACATAACCACCTCCCCTGTCTGCACTCACCGGAACAGAGACTCCACCAG",7,pg,lstA,lstT,lstC,lstG)# with 5 mutations, 1 length>1
#     s,p = alignment("CACATAACCACCTCCCCTGTCTGCACTCACCGGAACCTCTGAGCCACCAG",7,pg,lstA,lstT,lstC,lstG)#5 mutations, 2 length >1, 1 word off
#     s,p=alignment("CCCTACAACCACCACAAATGTCTCCACTCACCGGAACAGAGACTTTTCCAG",7,pg,lstA,lstT,lstC,lstG)#3 mutations lengths>2
#     w,s=alignment("CCCCTAACCACCACAAATGTCTCCACTCACCGGAACAGAGACTTTTCCAG",7,pg,lstA,lstT,lstC,lstG)#first mutation length 4->2, rest same as prev
#     s,p=alignment("CCCCTAACCACCACAAAACATTCCACTCACCGGAACAGAGACTTTTCCAG",7,pg,lstA,lstT,lstC,lstG)#turned second mutation into a whole word
    mgf(pg,p,"CACATAACCACCTCCCCTGTCTGCACTCACCGGAACAGAGACTCCACCAG")
    end = timer()
    print("time: {}".format(end - start))

if __name__ =="__main__":
    main()


# In[ ]:





# In[ ]:




