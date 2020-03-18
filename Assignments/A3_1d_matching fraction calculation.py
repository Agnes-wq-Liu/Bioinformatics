#1d: give fraction perfect match of both ends
# match start not end
#match end not start
#neither
#for annotated fraction
# import biopython
import pandas as pd
annot =pd.read_csv('/root/comp561/Vibrio_vulnificus.ASM74310v1.37.gff3',skiprows = 394,header = None,sep = '\t')
# annot = pd.read_csv('/root/comp561/Vibrio_vulnificus.ASM74310v1.37.gff3',header = None, sep = '\t')
mygf = pd.read_csv('/root/comp561/1c_output.gff3',header = None,sep = '\t')

annot = annot[annot.iloc[:,2] == "CDS"]
annot = annot[annot.iloc[:,6]=="+"]
annot_s = list(annot.iloc[:,3].values)
annot_t = list(annot.iloc[:,4].values)
mygf_s = list(mygf.iloc[:,3].values)
mygf_t = list(mygf.iloc[:,4].values)
ab =0
ast = 0
ae = 0
an = 0
mb = 0
mst = 0
me = 0
mn=0
for i in range (0,len(annot_s),1):
    if (int(annot_s[i]) in mygf_s) and (int(annot_t[i]) in mygf_t):
        ab+=1
    elif int(annot_s[i]) in mygf_s:
        ast+=1
    elif int(annot_t[i]) in mygf_t:
        ae+=1
    else:
        an+=1
ab = float(ab)/float((len(annot_s)))
ast = float(ast)/float((len(annot_s)))
ae = float(ae)/float((len(annot_s)))
an = float(an)/float((len(annot_s)))
with open ("/root/comp561/1d_output.txt","w") as d:
    d.writelines("Fractions based on annotation file:\n")
with open ("/root/comp561/1d_output.txt","a+") as d:
    d.writelines("\nmatch on both coordinates is: " +str(ab)+"\n")
    d.writelines("match on start only is: " +str(ast)+"\n")
    d.writelines("match on end only is: " +str(ae)+"\n")
    d.writelines("match on neither coordinates is: " +str(an)+"\n")
    
for i in range (0,len(mygf_s),1):
    if (float(mygf_s[i]) in annot_s) and (float(mygf_t[i]) in annot_t) :
        mb+=1
    elif float(mygf_s[i]) in annot_s:
        mst+=1
    elif float(mygf_t[i]) in annot_t:
        me+=1
    else:
        mn+=1
mb = float(mb)/float((len(mygf_s)))
mst = float(mst)/float((len(mygf_s)))
me = float(me)/float((len(mygf_s)))
mn = float(mn)/float((len(mygf_s)))
with open ("/root/comp561/1d_output.txt","a+") as d:
    d.writelines("\n\nthe fractions based on my prediction file:\n")
with open ("/root/comp561/1d_output.txt","a+") as d:
    d.writelines("\nmatch on both coordinates is: " +str(mb)+"\n")
    d.writelines("match on start only is: " +str(mst)+"\n")
    d.writelines("match on end only is: " +str(me)+"\n")
    d.writelines("match on neither coordinates is: " +str(mn)+"\n")
