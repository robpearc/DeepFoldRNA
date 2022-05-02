import numpy as np 
import sys

inputfile=sys.argv[1]
def np2txt(npy,txtfile):
    L,_,D=npy.shape
    lines=[]
    for i in range(L):
        for j in range(L):
            #start from the first line
            aline=str(i+1)+' '+str(j+1)+' '
            values=[str(avalue) for avalue in npy[i,j]]
            aline+=' '.join(values)+'\n'
            lines.append(aline)
    wfile=open(txtfile,'w')
    wfile.writelines(lines)
    wfile.close()      

def np2txt1d(npy,txtfile):
    L,D=npy.shape
    lines=[]
    for i in range(L):
        #start from the first line
        aline=str(i+1)+' '
        values=[str(avalue) for avalue in npy[i]]
        aline+=' '.join(values)+'\n'
        lines.append(aline)
    wfile=open(txtfile,'w')
    wfile.writelines(lines)
    wfile.close()      

a=np.load(inputfile)
all_kys=list(a.keys())
print(all_kys)
for aky in all_kys:
    print(aky)
    savefile=inputfile+'_'+aky+'.txt'
    if aky=="eta_bb" or aky=="theta_bb":
        np2txt1d(a[aky],savefile)
    else:
        np2txt(a[aky],savefile)

