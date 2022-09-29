#-----
# forked from Lin's github
#
import sys
import time
import multiprocessing as mp
import numpy as np

import fmin_rgRPA_1p_1salt as fs
import thermal_rgRPA_1p_1salt as tt
import seq_list as sl


# small ion valencies
zs, zc = 1, 1

usesql = False
sqlname = 'uuuuuu'
sql_salt_const_name = 'salt_const_paras'


# fgRPA = True if calculating fg-RPA
fgRPA = False

if sys.argv[2] == "cri_calc":
   cri_calc_end = 1
else:
   cri_calc_end = 0
   umax = float(sys.argv[2])

duD = int(2)#2
du = 10**(-duD)

# Select a sequence in seq_list.py
seqname = sys.argv[1]


# reduced ealt concentration
#if len(sys.argv) > 3:
try:
    phis = float(sys.argv[3])
#else:
except:
    phis = 0

# eh, es: now it is a must before setup pH and w2r
try: 
    tt.eh = float(sys.argv[4])
    tt.es = float(sys.argv[5])    
except:
    tt.eh = 0
    tt.es = 0

# pH value
pHchange = False
pHvalue = 0
try:
    pHvalue = float(sys.argv[6])
    sig, N, the_seq = sl.get_the_charge(seqname, pH=pHvalue )
    pHchange = True
except:
    sig, N, the_seq = sl.get_the_charge(seqname)   

# short-range repulsion
try:
    wtwo_ratio = float(sys.argv[7])
    HP = tt.Heteropolymer(sig, zc=zc,zs=zs, w2=4*np.pi/3*wtwo_ratio)
except:
    wtwo_ratio = 1
    HP = tt.Heteropolymer(sig, zc=zc,zs=zs)


# doing fg-RPA
if fgRPA:
    tt.useleff=False
    HP = tt.Heteropolymer(sig,w2=0)


#----------------------- Calculate critical point -----------------------
print('Seq:' , seqname, '=' , the_seq )

print('w2=', HP['w2'])
print('phis=', phis)
print('eh=' + str(tt.eh) + ' , es=' + str(tt.es) )

t0 = time.time()
phi_cri, u_cri = fs.cri_calc( HP, phis )

print('Critical point found in', time.time() - t0 , 's')

print( 'u_cri =', '{:.8e}'.format(u_cri) , \
       ', phi_cri =','{:.8e}'.format(phi_cri) )


cri_file = '../Critical_temp.txt'
with open(cri_file, 'a')as f: 
	print("%s %.3f %.4f %.4f"  % (seqname,1.0/u_cri, u_cri,phi_cri), file = f)
calc_info = '../_phis' + str(phis) + \
                '_w2r' + str(wtwo_ratio) + \
                '_eh' + str(tt.eh) + '_es' + str(tt.es) + \
                '.txt'
with open(calc_info, 'a')as f: 
	print("   " , file = f)
if(cri_calc_end):
    sys.exit();

#---------------------------- Set up u range ----------------------------
ddu = du/10
umin = (np.floor(u_cri/ddu)+1)*ddu
uclose = (np.floor(u_cri/du)+2)*du
if umax < u_cri:
    umax = np.floor(u_cri*1.5) 
if uclose > umax:
    uclose = umax
print(umax)
uall = np.append( np.arange(umin, uclose, ddu ), \
                  np.arange(uclose, umax+du, du ) )
total_index=len(uall)

#-------------------- calculate multiple u's -------------------

def bisp(u):
    sp1, sp2 = fs.ps_sp_solve(HP, phis, u, phi_cri)
    print( u, sp1, sp2, 'sp done!', flush=True)
    bi1, bi2 = fs.ps_bi_solve(HP, phis, u, (sp1, sp2), phi_cri )
    print( u, bi1, bi2, 'bi done!', flush=True)
    return sp1, sp2, bi1, bi2
the_index=0
sp1=np.repeat(0.0, total_index)
sp2=np.repeat(0.0, total_index)
bi1=np.repeat(0.0, total_index)
bi2=np.repeat(0.0, total_index)
for i in range(0, total_index):
    sp1[i], sp2[i], bi1[i], bi2[i] =bisp(uall[i])
     

#---------------------- Prepare for output ----------------------

ind = np.where(np.array(sp1) > fs.phi_min_calc)[0]
nout = ind.shape[0]
sp1t = np.array([ sp1[i] for i in ind])
sp2t = np.array([ sp2[i] for i in ind])
bi1t = np.array([ bi1[i] for i in ind])
bi2t = np.array([ bi2[i] for i in ind])
ut = np.array([ round(uu, duD+1) for uu in uall[ind]])
new_umax = np.max(ut)

sp_out = np.zeros((2*nout+1,2))
bi_out = np.zeros((2*nout+1,2))

sp_out[:,1] = np.concatenate((ut[::-1], [u_cri], ut))
sp_out[:,0] = np.concatenate((sp1t[::-1], [phi_cri], sp2t))
bi_out[:,1] = sp_out[:,1]
bi_out[:,0] = np.concatenate((bi1t[::-1], [phi_cri], bi2t))

calc_info = seqname + '.txt'

sp_file = './results/sp_' + calc_info
bi_file = './results/bi_' + calc_info


cri_info = "u_cri= " + str(u_cri) + " , phi_cri= " + str(phi_cri) + \
           "\n------------------------------------------------------------"           

np.savetxt(sp_file, sp_out, fmt = '%.10e', header= cri_info )
np.savetxt(bi_file, bi_out, fmt = '%.10e', header= cri_info )

    



