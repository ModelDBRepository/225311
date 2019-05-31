# calcsumcurr_manyareagsynmediumtau_parts_fixeddt.py
# Python script for running a simulation of L5PC under random synaptic inputs. The script saves transmembrane currents that can be later used for calculating
# LFP and diffusion in the extracellular domain.
# Tuomo Maki-Marttunen, 2014-2016
#
# Usage:
# python calcsumcurr_manyareagsynmediumtau_parts_fixeddt.py 20 0.025 0.000042 10000 10000 2 1 200 #Run a single pyramidal cell simulation with 20 segments/compartment,
#                                                                                                 #0.025ms time step, 0.042 nS synaptic conductance, 10000 ms biological
#                                                                                                 #time, with 10000 synapses of type 2 (scattered on both apical and
#                                                                                                 #basal dendrites), random number seed 1, and performing the simulation
#                                                                                                 #in units of 200ms
#
# Output:
#   The script saves the amounts of each current type exiting the cell into the 13 extracellular compartments. The bordering extracellular compartments do not contain any
#   membrane segments, and thus 13 compartments are enough. The results are saved into MATLAB files
#
#   currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-5_tstop10000_nseg20_dt0.025_seed1_sim0x200.mat
#   currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-5_tstop10000_nseg20_dt0.025_seed1_sim1x200.mat
#   ...
#   currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-5_tstop10000_nseg20_dt0.025_seed1_sim49x200.mat



from neuron import h
import numpy
import scipy.io
from pylab import *
import time
import sys
from os.path import exists

morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"
v0 = -80
ca0 = 0.0001

synlambda = 5.0 # frequency of synaptic inputs (Hz)
syntau = 2.0    # decay time (ms)

proximalpoint = 400
distalpoint = 620
#distalpoint = 960
BACdt = 5.0
part_threshys = [100*x+68 for x in range(-1,11)]+[inf]
if len(sys.argv) > 1:
  nsegs = int(float(sys.argv[1]))
else:
  nsegs = 20
if len(sys.argv) > 2:
  dt = float(sys.argv[2])
else:
  dt = 0.025
if len(sys.argv) > 3:
  syngmax = float(sys.argv[3])
else:
  syngmax = 0.000042
if len(sys.argv) > 4:
  tstop = float(sys.argv[4])
else:
  tstop = 10000
if len(sys.argv) > 5:
  Nsynlocs = int(sys.argv[5])
else:
  Nsynlocs = 10000
if len(sys.argv) > 6:
  synloctype = int(sys.argv[6])
else:
  synloctype = 2 # 1: apical only, 2: apical & basal, 3: basal only
if len(sys.argv) > 7:
  myseed = int(sys.argv[7])
else:
  myseed = 1
if len(sys.argv) > 8:
  singleSimT = float(sys.argv[8])
else:
  singleSimT = 200

seed(myseed)
Nsims = int(1.0*tstop/singleSimT+0.9999)
dtsave = 0.1

# If the simulation is already done, exit
if exists('currsums_parts_'+str(Nsynlocs)+'areagsynsmediumtau_fixeddt_type'+str(synloctype)+'_amp'+str(syngmax)+'_tstop'+str(tstop)+'_nseg'+str(nsegs)+'_dt'+str(dt)+'_seed'+str(myseed)+'_sim'+str(Nsims-1)+'x'+str(singleSimT)+'.mat'):
  sys.exit()

# Initialize the model
h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
dtsave = """+str(dtsave)+"""
objref st1, synlist
st1 = new IClamp(0.5)
st1.del = 700
L5PC.soma st1
synlist = new List()
objref isyn,tvec,sl
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
access L5PC.apic[siteVec[0]]
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
{vsoma = new Vector()}
{casoma = new Vector()}
{vdend = new Vector()}
{cadend = new Vector()}
{vdend2 = new Vector()}
{cadend2 = new Vector()}
access L5PC.soma
{vsoma.record(&v(0.5),dtsave)}
{casoma.record(&cai(0.5),dtsave)}
access L5PC.apic[siteVec[0]]
{vdend.record(&v(siteVec[1]),dtsave)}
{cadend.record(&cai(siteVec[1]),dtsave)}
access L5PC.soma
{isoma = new Vector()}
{isoma.record(&st1.i,dtsave)}
forsec L5PC.all {
  nseg = """+str(nsegs)+"""
}
dt = """+str(dt)+"""
""")

# Initialize the variables where the transmembrane currents will be saved
h("""
objref apicalina["""+str(109*nsegs)+"""], apicalik["""+str(109*nsegs)+"""], apicalica["""+str(109*nsegs)+"""], apicalih["""+str(109*nsegs)+"""], apicalil["""+str(109*nsegs)+"""], apicalv["""+str(109*nsegs)+"""], apicalicap["""+str(109*nsegs)+"""], apicalimemb["""+str(109*nsegs)+"""]
objref basalih["""+str(84*nsegs)+"""], basalil["""+str(84*nsegs)+"""], basalv["""+str(84*nsegs)+"""], basalicap["""+str(84*nsegs)+"""], basalimemb["""+str(84*nsegs)+"""]
objref somaticina["""+str(nsegs)+"""], somaticik["""+str(nsegs)+"""], somaticica["""+str(nsegs)+"""], somaticih["""+str(nsegs)+"""], somaticil["""+str(nsegs)+"""], somaticv["""+str(nsegs)+"""], somaticicap["""+str(nsegs)+"""], somaticimemb["""+str(nsegs)+"""]
objref axonalil["""+str(2*nsegs)+"""], axonalv["""+str(2*nsegs)+"""], axonalicap["""+str(2*nsegs)+"""], axonalimemb["""+str(2*nsegs)+"""]
""")

# Initialize the segment areas
print "objref complete"
A_apical = [0]*109*nsegs
A_basal = [0]*84*nsegs
A_somatic = [0]*nsegs
A_axonal = [0]*2*nsegs
part_apical = [0]*109*nsegs
part_basal = [0]*84*nsegs
part_somatic = [0]*nsegs
part_axonal = [0]*2*nsegs
len_apical = [0]*109*nsegs
len_basal = [0]*84*nsegs
len_somatic = [0]*nsegs
len_axonal = [0]*2*nsegs
maxx = -inf
minx = inf
maxy = -inf
miny = inf
maxz = -inf
minz = inf

# Set the recordings for the compartments in the apical dendrite. Calculate also the membrane areas of each segment
for i in range(0,109):
  h("access L5PC.apic["+str(i)+"]")
  h("tmpvarx = x3d(0)")
  h("tmpvary = y3d(0)")
  h("tmpvarz = z3d(0)")
  h("tmpvarx2 = x3d(n3d()-1)")
  h("tmpvary2 = y3d(n3d()-1)")
  h("tmpvarz2 = z3d(n3d()-1)")
  coord1 = [h.tmpvarx,h.tmpvary,h.tmpvarz]
  coord2 = [h.tmpvarx2,h.tmpvary2,h.tmpvarz2]
  for j in range(0,nsegs):
    thispos = 0.5/nsegs+1.0/nsegs*j
    if nsegs == 1:
      thispos3d = [(x + y)/2 for x,y in zip(coord1,coord2)]
    else:
      thispos3d = [x + j*(y-x)/(nsegs-1) for x,y in zip(coord1,coord2)]
    part_apical[i*nsegs+j] = next(i for i,x in enumerate(part_threshys) if thispos3d[1] <= x)
    h("{apicalina["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalik["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalica["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalih["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalil["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalv["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalicap["+str(i*nsegs+j)+"] = new Vector()}")
    h("{apicalimemb["+str(i*nsegs+j)+"] = new Vector()}")
    h("L5PC.apic["+str(i)+"] insert extracellular")
    h("access L5PC.apic["+str(i)+"]")
    h("{apicalina["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].ina("+str(thispos)+    "),dtsave)}")
    h("{apicalik["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].ik("+str(thispos)+     "),dtsave)}")
    h("{apicalica["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].ica("+str(thispos)+    "),dtsave)}")
    h("{apicalih["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].ihcn_Ih("+str(thispos)+"),dtsave)}")
    h("{apicalil["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].i_pas("+str(thispos)+  "),dtsave)}")
    h("{apicalv["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].v("+str(thispos)+      "),dtsave)}")
    h("{apicalicap["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].i_cap("+str(thispos)+  "),dtsave)}")
    h("{apicalimemb["+str(i*nsegs+j)+"].record(&L5PC.apic["+str(i)+"].i_membrane("+str(thispos)+  "),dtsave)}")
    h("L5PC.apic["+str(i)+"].nseg = " + str(nsegs))
    h("tmpvar = area("+str(thispos)+")")
    A_apical[i*nsegs+j] = h.tmpvar
  h("tmpvar = L")
  len_apical[i] = h.tmpvar
print "apical complete"

# Set the recordings for the compartments in the basal dendrite. Calculate also the membrane areas of each segment
for i in range(0,84):
  h("access L5PC.dend["+str(i)+"]")
  h("tmpvarx = x3d(0)")
  h("tmpvary = y3d(0)")
  h("tmpvarz = z3d(0)")
  h("tmpvarx2 = x3d(n3d()-1)")
  h("tmpvary2 = y3d(n3d()-1)")
  h("tmpvarz2 = z3d(n3d()-1)")
  coord1 = [h.tmpvarx,h.tmpvary,h.tmpvarz]
  coord2 = [h.tmpvarx2,h.tmpvary2,h.tmpvarz2]
  for j in range(0,nsegs):
    thispos = 0.5/nsegs+1.0/nsegs*j
    if nsegs == 1:
      thispos3d = [(x + y)/2 for x,y in zip(coord1,coord2)]
    else:
      thispos3d = [x + j*(y-x)/(nsegs-1) for x,y in zip(coord1,coord2)]
    part_basal[i*nsegs+j] = next(i for i,x in enumerate(part_threshys) if thispos3d[1] <= x)
    h("{basalih["+str(i*nsegs+j)+"] = new Vector()}")
    h("{basalil["+str(i*nsegs+j)+"] = new Vector()}")
    h("{basalv["+str(i*nsegs+j)+"] = new Vector()}")
    h("{basalicap["+str(i*nsegs+j)+"] = new Vector()}")
    h("{basalimemb["+str(i*nsegs+j)+"] = new Vector()}")
    h("L5PC.dend["+str(i)+"] insert extracellular")
    h("access L5PC.dend["+str(i)+"]")
    h("{basalih["+str(i*nsegs+j)+"].record(&L5PC.dend["+str(i)+"].ihcn_Ih("+str(thispos)+"),dtsave)}")
    h("{basalil["+str(i*nsegs+j)+"].record(&L5PC.dend["+str(i)+"].i_pas("+str(thispos)+  "),dtsave)}")
    h("{basalv["+str(i*nsegs+j)+"].record(&L5PC.dend["+str(i)+"].v("+str(thispos)+      "),dtsave)}")
    h("{basalicap["+str(i*nsegs+j)+"].record(&L5PC.dend["+str(i)+"].i_cap("+str(thispos)+  "),dtsave)}")
    h("{basalimemb["+str(i*nsegs+j)+"].record(&L5PC.dend["+str(i)+"].i_membrane("+str(thispos)+  "),dtsave)}")
    h("L5PC.dend["+str(i)+"].nseg = " + str(nsegs))
    h("tmpvar = area("+str(thispos)+")")
    A_basal[i*nsegs+j] = h.tmpvar
  h("tmpvar = L")
  len_basal[i] = h.tmpvar
print "basal complete"

# Set the recordings for the compartments in the soma. Calculate also the membrane area
for i in range(0,1):
  h("access L5PC.soma["+str(i)+"]")
  h("tmpvarx = x3d(0)")
  h("tmpvary = y3d(0)")
  h("tmpvarz = z3d(0)")
  h("tmpvarx2 = x3d(n3d()-1)")
  h("tmpvary2 = y3d(n3d()-1)")
  h("tmpvarz2 = z3d(n3d()-1)")
  coord1 = [h.tmpvarx,h.tmpvary,h.tmpvarz]
  coord2 = [h.tmpvarx2,h.tmpvary2,h.tmpvarz2]
  for j in range(0,nsegs):
    thispos = 0.5/nsegs+1.0/nsegs*j
    if nsegs == 1:
      thispos3d = [(x + y)/2 for x,y in zip(coord1,coord2)]
    else:
      thispos3d = [x + j*(y-x)/(nsegs-1) for x,y in zip(coord1,coord2)]
    part_somatic[i*nsegs+j] = next(i for i,x in enumerate(part_threshys) if thispos3d[1] <= x)
    h("{somaticina["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticik["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticica["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticih["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticil["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticv["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticicap["+str(i*nsegs+j)+"] = new Vector()}")
    h("{somaticimemb["+str(i*nsegs+j)+"] = new Vector()}")
    h("L5PC.soma["+str(i)+"] insert extracellular")
    h("access L5PC.soma["+str(i)+"]")
    h("{somaticina["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].ina("+str(thispos)+    "),dtsave)}")
    h("{somaticik["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].ik("+str(thispos)+     "),dtsave)}")
    h("{somaticica["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].ica("+str(thispos)+    "),dtsave)}")
    h("{somaticih["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].ihcn_Ih("+str(thispos)+"),dtsave)}")
    h("{somaticil["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].i_pas("+str(thispos)+  "),dtsave)}")
    h("{somaticv["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].v("+str(thispos)+      "),dtsave)}")
    h("{somaticicap["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].i_cap("+str(thispos)+  "),dtsave)}")
    h("{somaticimemb["+str(i*nsegs+j)+"].record(&L5PC.soma["+str(i)+"].i_membrane("+str(thispos)+  "),dtsave)}")
    h("L5PC.soma["+str(i)+"].nseg = " + str(nsegs))
    h("tmpvar = area("+str(thispos)+")")
    A_somatic[i*nsegs+j] = h.tmpvar
    if j == nsegs/2:
      somapos3d = thispos3d[:]
  h("tmpvar = L")
  len_somatic[i] = h.tmpvar
print "somatic complete"

# Set the recordings for the compartments in the axon initial segment. Calculate also the membrane areas of each segment
for i in range(0,2):
  for j in range(0,nsegs):
    thispos = 0.5/nsegs+1.0/nsegs*j
    part_axonal[i*nsegs+j] = next(i for i,x in enumerate(part_threshys) if somapos3d[1] <= x)
    h("{axonalil["+str(i*nsegs+j)+"] = new Vector()}")
    h("{axonalv["+str(i*nsegs+j)+"] = new Vector()}")
    h("{axonalicap["+str(i*nsegs+j)+"] = new Vector()}")
    h("{axonalimemb["+str(i*nsegs+j)+"] = new Vector()}")
    h("L5PC.axon["+str(i)+"] insert extracellular")
    h("access L5PC.axon["+str(i)+"]")
    h("{axonalil["+str(i*nsegs+j)+"].record(&L5PC.axon["+str(i)+"].i_pas("+str(thispos)+"),dtsave)}")
    h("{axonalv["+str(i*nsegs+j)+"].record(&L5PC.axon["+str(i)+"].v("+str(thispos)+    "),dtsave)}")
    h("{axonalicap["+str(i*nsegs+j)+"].record(&L5PC.axon["+str(i)+"].i_cap("+str(thispos)+    "),dtsave)}")
    h("{axonalimemb["+str(i*nsegs+j)+"].record(&L5PC.axon["+str(i)+"].i_membrane("+str(thispos)+    "),dtsave)}")
    h("L5PC.axon["+str(i)+"].nseg = " + str(nsegs))
    h("tmpvar = area("+str(thispos)+")")
    A_axonal[i*nsegs+j] = h.tmpvar
  h("tmpvar = L")
  len_axonal[i] = h.tmpvar
print "axonal complete"

synbranch = [0]*Nsynlocs
syniseg = [0]*Nsynlocs
synx = [0.0]*Nsynlocs

# Calculate the probability of synapse being found in the basal dendrite.
if synloctype == 1:
  basalprob = 0.0
if synloctype == 2:
  basalprob = sum(A_basal)/(sum(A_basal) + sum(A_apical))
if synloctype == 3:
  basalprob = 1.0
print "Basal area: "+str(sum(A_basal))
print "Apical area: "+str(sum(A_apical))
print "basalprob = "+str(basalprob)

# Calculate the probabilities for the synapse being in each segment
ps_basal = [1.0*x/sum(A_basal) for x in A_basal]
cumps_basal = cumsum(ps_basal)
ps_apical = [1.0*x/sum(A_apical) for x in A_apical]
cumps_apical = cumsum(ps_basal)

# Draw the random numbers, one to decide which branch, one to decide which segment, and one to determine the distance x from 0-end
rs_branch = rand(Nsynlocs)
rs_seg = rand(Nsynlocs)
rs_x = rand(Nsynlocs)

ts_syn = []
seg_syn = [-1]*Nsynlocs
seg_syn_accurate = [-1]*Nsynlocs
part_syn = [-1]*Nsynlocs

# For each synapse, determine to which extracellular compartment it outputs the currents, and randomize the synapse activation times
# To do: Might become faster if the set of AlphaSynapses at a single synaptic location is replaced by a single point process that
# is activated at the time instants drawn here.
for isyn in range(0,Nsynlocs):
  if rs_branch[isyn] <= basalprob:
    synbranch[isyn] = 1
  if synbranch[isyn] == 0:
    seg_syn_accurate[isyn] = next((i for i,x in enumerate(cumps_apical) if x > rs_seg[isyn]))
    seg_syn[isyn] = seg_syn_accurate[isyn]/nsegs
    segnum = seg_syn_accurate[isyn] % nsegs
    h("access L5PC.apic["+str(seg_syn[isyn])+"]")
    mystr = "L5PC.apic["+str(seg_syn[isyn])+"]"
    part_syn[isyn] = part_apical[seg_syn[isyn]*nsegs+segnum]
  else:
    seg_syn_accurate[isyn] = next((i for i,x in enumerate(cumps_basal) if x > rs_seg[isyn]))
    seg_syn[isyn] = seg_syn_accurate[isyn]/nsegs
    segnum = seg_syn_accurate[isyn] % nsegs
    h("access L5PC.dend["+str(seg_syn[isyn])+"]")
    mystr = "L5PC.dend["+str(seg_syn[isyn])+"]"
    part_syn[isyn] = part_basal[seg_syn[isyn]*nsegs+segnum]
  ts = []
  t = 0
  secx = 1.0*segnum/nsegs+1.0/nsegs*rs_x[isyn]
  while t < tstop:
    t = t - 1000.0/synlambda*log(1-rand())
    ts.append(t)
    h("{synlist.append(new AlphaSynapse("+str(secx)+"))}")
    h("syni = synlist.count()-1")
    h("synlist.o[syni].tau = "+str(syntau))
    h("synlist.o[syni].gmax = "+str(syngmax))
    h("synlist.o[syni].e = 0")
    h("synlist.o[syni].onset = "+str(t))
  ts_syn.append(ts[:])


Nsyns = h.syni+1

h("""
tstop = """+str(tstop)+"""
v_init = """+str(v0)+"""
cai0_ca_ion = """+str(ca0)+"""
st1.amp = 0
st1.dur = 0
""")

print "Initializing..."
h("stdinit()")
print "Init complete"
tfin = 0


# Repeat the simulation of singleSimT milliseconds Nsims times. Save each simulation to a MATLAB file.
for isim in range(0,Nsims):
  h("{tvec.resize(0)}")
  h("{vsoma.resize(0)}")
  h("{vdend.resize(0)}")
  h("{casoma.resize(0)}")
  h("{cadend.resize(0)}")
  for i in range(0,109*nsegs):
    h("{apicalina["+str(i)+"].resize(0)}")
    h("{apicalik["+str(i)+"].resize(0)}")
    h("{apicalica["+str(i)+"].resize(0)}")
    h("{apicalih["+str(i)+"].resize(0)}")
    h("{apicalil["+str(i)+"].resize(0)}")
    h("{apicalv["+str(i)+"].resize(0)}")
    h("{apicalicap["+str(i)+"].resize(0)}")
    h("{apicalimemb["+str(i)+"].resize(0)}")
  for i in range(0,84*nsegs):
    h("{basalih["+str(i)+"].resize(0)}")
    h("{basalil["+str(i)+"].resize(0)}")
    h("{basalv["+str(i)+"].resize(0)}")
    h("{basalicap["+str(i)+"].resize(0)}")
    h("{basalimemb["+str(i)+"].resize(0)}")
  for i in range(0,nsegs):
    h("{somaticina["+str(i)+"].resize(0)}")
    h("{somaticik["+str(i)+"].resize(0)}")
    h("{somaticica["+str(i)+"].resize(0)}")
    h("{somaticih["+str(i)+"].resize(0)}")
    h("{somaticil["+str(i)+"].resize(0)}")
    h("{somaticv["+str(i)+"].resize(0)}")
    h("{somaticicap["+str(i)+"].resize(0)}")
    h("{somaticimemb["+str(i)+"].resize(0)}")
  for i in range(0,2*nsegs):
    h("{axonalil["+str(i)+"].resize(0)}")
    h("{axonalv["+str(i)+"].resize(0)}")
    h("{axonalicap["+str(i)+"].resize(0)}")
    h("{axonalimemb["+str(i)+"].resize(0)}")
  tfin = tfin+singleSimT
  print "Starting run "+str(isim)+" until "+str(tfin)+" ms"
  h("continuerun("+str(tfin)+")")
  print "Run "+str(isim)+" complete, tfin = "+str(tfin)

  Vsoma=np.array(h.vsoma)
  Vdend=np.array(h.vdend)
  Casoma=np.array(h.casoma)
  Cadend=np.array(h.cadend)
  times=np.array([tfin-singleSimT+dtsave*i for i in range(0,len(Vsoma))])

  # Here, we calculate the extracellular compartment -wise currents of each species. We will need the previously saved data on the areas of each dendritic compartment and the 
  # extracellular compartment to which each dendritic compartment belongs to (in addition to the transmembrane currents saved during the simulation) in order to do this.
  ina_all = [[]]*len(part_threshys)
  ik_all = [[]]*len(part_threshys)
  ica_all = [[]]*len(part_threshys)
  ih_all = [[]]*len(part_threshys)
  il_all = [[]]*len(part_threshys)
  v_all = [[]]*len(part_threshys)
  icap_all = [[]]*len(part_threshys)
  imemb_all = [[]]*len(part_threshys)
  A_tot_all = [[]]*len(part_threshys)
  for ipart in range(0,len(part_threshys)): # Loop through the extracellular compartments
    ina = [0.0]*len(times)
    ik = [0.0]*len(times)
    ica = [0.0]*len(times)
    ih = [0.0]*len(times)
    il = [0.0]*len(times)
    v = [0.0]*len(times)
    icap = [0.0]*len(times)
    imemb = [0.0]*len(times)
    A_tot = sum([x for x,y in zip(A_apical,part_apical) if y==ipart]) + sum([x for x,y in zip(A_basal,part_basal) if y==ipart]) + sum([x for x,y in zip(A_somatic,part_somatic) if y==ipart]) + sum([x for x,y in zip(A_axonal,part_axonal) if y==ipart])
    for i in range(0,109*nsegs): # Loop through apical dendritic compartments
      if part_apical[i]==ipart:
        ina = [x+y for x,y in zip(ina,[A_apical[i]*x for x in np.array(h.apicalina[i])])]
        ik = [x+y for x,y in zip(ik,[A_apical[i]*x for x in np.array(h.apicalik[i])])]
        ica = [x+y for x,y in zip(ica,[A_apical[i]*x for x in np.array(h.apicalica[i])])]
        ih = [x+y for x,y in zip(ih,[A_apical[i]*x for x in np.array(h.apicalih[i])])]
        il = [x+y for x,y in zip(il,[A_apical[i]*x for x in np.array(h.apicalil[i])])]
        v = [x+y for x,y in zip(v,[A_apical[i]*x for x in np.array(h.apicalv[i])])]
        icap = [x+y for x,y in zip(icap,[A_apical[i]*x for x in np.array(h.apicalicap[i])])]
        imemb = [x+y for x,y in zip(imemb,[A_apical[i]*x for x in np.array(h.apicalimemb[i])])]
    for i in range(0,84*nsegs): # Loop through basal dendritic compartments
      if part_basal[i]==ipart:
        ih = [x+y for x,y in zip(ih,[A_basal[i]*x for x in np.array(h.basalih[i])])]
        il = [x+y for x,y in zip(il,[A_basal[i]*x for x in np.array(h.basalil[i])])]
        v = [x+y for x,y in zip(v,[A_basal[i]*x for x in np.array(h.basalv[i])])]
        icap = [x+y for x,y in zip(icap,[A_basal[i]*x for x in np.array(h.basalicap[i])])]
        imemb = [x+y for x,y in zip(imemb,[A_basal[i]*x for x in np.array(h.basalimemb[i])])]
    for i in range(0,nsegs): # Loop through somatic compartment
      if part_somatic[i]==ipart:
        ina = [x+y for x,y in zip(ina,[A_somatic[i]*x for x in np.array(h.somaticina[i])])]
        ik = [x+y for x,y in zip(ik,[A_somatic[i]*x for x in np.array(h.somaticik[i])])]
        ica = [x+y for x,y in zip(ica,[A_somatic[i]*x for x in np.array(h.somaticica[i])])]
        ih = [x+y for x,y in zip(ih,[A_somatic[i]*x for x in np.array(h.somaticih[i])])]
        il = [x+y for x,y in zip(il,[A_somatic[i]*x for x in np.array(h.somaticil[i])])]
        v = [x+y for x,y in zip(v,[A_somatic[i]*x for x in np.array(h.somaticv[i])])]
        icap = [x+y for x,y in zip(icap,[A_somatic[i]*x for x in np.array(h.somaticicap[i])])]
        imemb = [x+y for x,y in zip(imemb,[A_somatic[i]*x for x in np.array(h.somaticimemb[i])])]
    for i in range(0,2*nsegs): # Loop through axonal compartments
      if part_axonal[i]==ipart:
        il = [x+y for x,y in zip(il,[A_axonal[i]*x for x in np.array(h.axonalil[i])])]
        v = [x+y for x,y in zip(v,[A_axonal[i]*x for x in np.array(h.axonalv[i])])]
        icap = [x+y for x,y in zip(icap,[A_axonal[i]*x for x in np.array(h.axonalicap[i])])]
        imemb = [x+y for x,y in zip(imemb,[A_axonal[i]*x for x in np.array(h.axonalimemb[i])])]
    ina_all[ipart] = [x/100.0 for x in ina[:]]
    ik_all[ipart] = [x/100.0 for x in ik[:]]
    ica_all[ipart] = [x/100.0 for x in ica[:]]
    ih_all[ipart] = [x/100.0 for x in ih[:]]
    il_all[ipart] = [x/100.0 for x in il[:]]
    v_all[ipart] = v[:]
    icap_all[ipart] = [x/100.0 for x in icap[:]]
    imemb_all[ipart] = [x/100.0 for x in imemb[:]]
    A_tot_all[ipart] = A_tot

  dict = {'ina': numpy.array(ina_all), 'ik': numpy.array(ik_all),
          'ica': numpy.array(ica_all), 'ih': numpy.array(ih_all),
          'il': numpy.array(il_all), 'VtimesA': numpy.array(v_all), 'imemb': numpy.array(imemb_all), 'Vsoma': numpy.array(Vsoma),
          'icap': numpy.array(icap_all), 'times': numpy.array(times), 'A': numpy.array(A_tot_all),
          'ts_syn': numpy.array(ts_syn), 'part_syn': numpy.array([1+x for x in part_syn])}
  scipy.io.savemat('currsums_parts_'+str(Nsynlocs)+'areagsynsmediumtau_fixeddt_type'+str(synloctype)+'_amp'+str(syngmax)+'_tstop'+str(tstop)+'_nseg'+str(nsegs)+'_dt'+str(dt)+'_seed'+str(myseed)+'_sim'+str(isim)+'x'+str(singleSimT)+'.mat', dict)
