#!/usr/bin/env desres-exec
#{
# exec desres-cleanenv -m Python/2.7.1-06A/bin -m simplejson/2.1.5-1/lib-python -- python $0 "$@"
#}

import os
import optparse
import sys
from decimal import Decimal
import simplejson
from collections import OrderedDict
pythonpath=os.path.join(os.path.dirname(sys.argv[0]), '..', 'python')
sys.path.append(pythonpath)
import DesFF

convaname={
  "DISU" : {"1CB":"CB","1SG":"SG","1HG1":"HG1","2CB":"BAD1","2SG":"BAD2","2HG1":"BAD3" },
  "TP1"  : {"1HH":"HH"},
  "TP2"  : {"1HH":"HH"},
  "SP1"  : {"1HG1":"HG1"},
  "SP2"  : {"1HG1":"HG1"},
  "THP1" : {"1HG1":"HG1"},
  "THP2" : {"1HG1":"HG1"},
  "Zn2+" : {"ZN" : "Zn2+"},
  "Cs+"  : {"CES" : "Cs+"},
  "Mg2+" : {"MG"  : "Mg2+"},
  "Cl-"  : {"CLA" : "Cl-"},
  "Ca2+" : {"CAL" : "Ca2+"},
  "Cd2+" : {"CD" : "Cd2+"},
  "Ba2+" : {"BAR" : "Ba2+"},
  "Li+"  : {"LIT" : "Li+"},
  "Na+"  : {"SOD" : "Na+"},
  "K+"   : {"POT" : "K+"},
  "Rb+"  : {"RUB" : "Rb+"}
  }

# Better residue names
convres={
  "HSE" : "HIE",
  "HSP" : "HIP",
  "HSD" : "HIS",
  "ZN2" : "Zn2+",
  "MG"  : "Mg2+",
  "CD2" : "Cd2+",
  "CLA" : "Cl-",
  "CAL" : "Ca2+",
  "BAR" : "Ba2+",
  "LIT" : "Li+",
  "SOD" : "Na+",
  "POT" : "K+",
  "RUB" : "Rb+",
  "CES" : "Cs+"
  }

# Better atom type names
convatype={
  "ZN"  : "Zn2+",
  "MG"  : "Mg2+",
  "CLA" : "Cl-",
  "CAL" : "Ca2+",
  "BAR" : "Ba2+",
  "CAD" : "Cd2+",
  "LIT" : "Li+",
  "SOD" : "Na+",
  "POT" : "K+",
  "RUB" : "Rb+",
  "CES" : "Cs+",
  "CT1x" : "CT1X",
  "CT2x" : "CT2X",
  "CT3x" : "CT3X",
  "ON2b" : "ON2B"
}

vfiles={}

FFconv     = ()
skipres    = ["TIP3","TP3M", "DUM"]

cmapnread  =0
cmapntotal =0
cmapdict   =dict()
_has_cmap  =False
atype2atomicnumber=dict()

#resname    =''
#residue    =''
#resq       =''
#residue_is_patch = 0
presdict   =dict()

# Container for unsatisfied valences.
# this is now populated from charmm DECL statements
bsplit={
  #  '-C' :'$1', # N termini amino acid incomplete valence
  #  '+N' :'$2', # C termini amino acid incomplete valence
  #  '+SG':'$3', # disulfide bond (CYS) incomplete valence
  #  '-O3':'$1', # 5' termini nucleic acid incomplete valence
  #  '+P' :'$2', # 3' termini nucleic acid incomplete valence
  #  '-C2':'$3'  # PEG monomer incomplete valence
  }

# conversion factors for charmm ff units into viparr ff units

conv_lj       = ( Decimal(repr(2.0 / 2.0**(1./6.))), -1 ) # sig,eps
conv_nbfix_lj = ( Decimal(repr(1.0 / 2.0**(1./6.))), -1 ) # sig,eps

def clean_line(l):
  l=l.strip()
  cstart=l.find('!')
  comment=""
  if(cstart<0):
    words=l.split()
  else:
    comment=l[cstart+1:].strip()
    words=l[:cstart].split()

  # viparr1 uses python to eval the json data. 
  # Python comments in eval'd strings cause problems
  comment=comment.replace('#','_') 
  return(words,comment)

def apply_noatom(atoms):
  for i,a in enumerate(atoms):
    if('+' == a[0] or '-' == a[0]):
      atoms[i]=bsplit[a]
  return atoms
  
def fix_aname(an,res):
  global convaname

  if (isinstance(an,str)):
    names=[an,]
  else:
    names=an

  newnames={}
  if convaname.has_key(res):
    newnames=convaname[res]
  for i in xrange(len(names)):
    names[i]=names[i].upper()
    if newnames.has_key(names[i]):
      names[i]=newnames[names[i]]
  
  if (isinstance(an,str)):
    an=names[0]
  else:
    an=names
  return an


def preload_merge(viparrdir):
  import elementdata
  global atype2atomicnumber
  
  if not os.path.exists(viparrdir) or not os.path.isdir(viparrdir):
    raise UserWarning("Bad base directory for viparr merge: "+viparrdir)
  file=os.path.join(viparrdir,"mass")
  if not os.path.exists(file):
    raise UserWarning("Couldnt find mass table: "+file)    
  fh=open(file,"r")
  lines=fh.read()
  fh.close()
  data=eval(lines)
  if(len(data)==0): return
  for d in data:
    if isinstance(d,list):
      atype=d[0]
      mass=d[1]
    elif isinstance(d,dict):
      atype=d["type"][0]
      mass=d["params"]["amu"]
    else:
      raise UserWarning("Bad mass file: %s \n  contents: %s"%(file,d))

    if(mass==0.0): continue
    element=atype[0]
    if(len(atype)>1): element+=atype[1].lower()
    if(element not in elementdata.element2number or abs(elementdata.data[element]["mass"][0]-float(mass))>0.1): element=element[0]
    if(element not in elementdata.element2number or abs(elementdata.data[element]["mass"][0]-float(mass))>0.1): 
      raise "Couldnt convert atomtype/mass to element: %s"%(str(d))

    atype2atomicnumber[atype]=elementdata.element2number[element]


def parse_stretch(lines):
  for l in lines[1:]:
    (words,comment)=clean_line(l)
    if(len(words)==0): continue

    atypes=[ (convatype[a] if convatype.has_key(a) else a) for a in words[0:2]] 
    params=OrderedDict([("r0",Decimal(words[3])),("fc",Decimal(words[2]))])

    FFconv.add_parameter(vfiles["bond"],atypes,params,comment)
  return 

def parse_angle(lines):
  for l in lines[1:]:
    
    (words,comment)=clean_line(l)
    if(len(words)==0): continue
    
    atypes=[ (convatype[a] if convatype.has_key(a) else a) for a in words[0:3]]
    params = OrderedDict([("theta0",Decimal(words[4])),("fc",Decimal(words[3]))])

    FFconv.add_parameter(vfiles["angle"],atypes,params,comment)

    if(len(words) == 7 ):
      params = OrderedDict([("r0",Decimal(words[6])),("fc",Decimal(words[5]))])

      FFconv.add_parameter(vfiles["urey"],atypes,params,comment)
    
  return

# Im adding this here because its hard to get right and at least charmm and amber need it
def collect_proper_trig(lastcollect,atypes, params, comment):
  from decimal import Decimal

  order=["phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"]
 
  phase=params[0]
  fc=params[1]
  pn=params[2]
  zero=Decimal("0.0")
   
  if (lastcollect!=atypes):
    data=[zero]*8
    data[0]=zero
    data[1]=fc
    if (phase==zero):
      data[1+pn]=fc
    elif (phase==Decimal("180.0")):
      data[1+pn]=-fc
    else:
      data[0]   =phase
      data[1+pn]=fc
    data=OrderedDict([ (k,v) for k,v in zip(order,data) ])
    FFconv.add_parameter(vfiles["proper"],atypes,data,comment)
  else:
    idxlist=FFconv.find_parameters(vfiles["proper"],atypes)
    idxlist.reverse()

    found=False
    lastidx=idxlist[0]+1
    for idx in idxlist:
      if idx != lastidx-1: break
      lastidx=idx
      data=FFconv.parameters[vfiles["proper"]].parameters[idx]["params"].values()
      if( (phase == zero or phase == Decimal("180.0")) and data[0] == zero ):
        found=True
        if data[1+pn] != zero:
          #if( (phase == zero and data[1+pn]==-fc) or (phase == Decimal("180.0") and data[1+pn]==-fc)): return
          raise UserWarning("multiplicity has already been assigned a value: %s %d %f %f"%(str(atypes),pn,data[1+pn],fc))
        data[1]+=fc
        if (phase==zero):
          data[1+pn]=fc
        elif (phase==Decimal("180.0")):
          data[1+pn]=-fc
        break
      # merge for phase != {0,180}
      elif(phase == data[0]):
        found=True
        if data[1+pn] != zero:
          #if(data[1+pn]==fc): return
          raise UserWarning("multiplicity has already been assigned a value: %s %d %f %f"%(str(atypes),pn,data[1+pn],fc))
        data[1]+=fc      
        data[1+pn]=fc
        break
    if(found):
      data=OrderedDict([ (k,v if v != zero else zero) for k,v in zip(order,data) ])
      FFconv.parameters[vfiles["proper"]].parameters[idx]["params"]=data
    else:
      data=[zero]*8
      data[0]=phase
      data[1]+=fc      
      data[1+pn]=fc
      data=OrderedDict([ (k,v) for k,v in zip(order,data) ])
      FFconv.add_parameter(vfiles["proper"],atypes,data,comment,check=False) 
      #data=[atypes,data,comment]
      #FFconv.parameters[vfiles["proper"]].append(data)
  return 

def parse_proper(lines):
  lastcollect=[]
  for l in lines[1:]:

    (words,comment)=clean_line(l)    
    if(len(words)==0): continue

    for i in xrange(4):
      words[i]=words[i].replace('X','*')

    fc     = Decimal(words[4])
    pn     = int(words[5])               
    phase  = Decimal(words[6])
    params = [phase, fc, pn]               

    atypes=[ (convatype[a] if convatype.has_key(a) else a) for a in words[0:4]]

    collect_proper_trig(lastcollect,atypes,params,comment)
    lastcollect=atypes
  return

def parse_improper(lines):
  for l in lines[1:]:
    (words,comment)=clean_line(l)    
    if(len(words)==0): continue    
    
    for i in xrange(4):
      words[i]=words[i].replace('X','*')

    fc    = Decimal(words[4])
    pn    = int(words[5])  # This is ignored
    phase = Decimal(words[6])
    params= OrderedDict([("phi0",phase),("fc", fc)])

    atypes=[ (convatype[a] if convatype.has_key(a) else a) for a in words[0:4]]

    FFconv.add_parameter(vfiles["improperH"],atypes,params,comment)
  return


def parse_cmap(lines):
  global cmapnread
  global cmapntotal
  global cmapdict

  for l in lines[1:]:
    (words,comment)=clean_line(l)    
    if(len(words)==0): continue   

    if(cmapntotal==0):

      if(len(words) !=9 ):
        raise UserWarning("Can't parse cmap term")

      npts=int(words[8])
      cmapntotal = npts*npts
      cmapnread=0
      cmapdict =dict()
      t1=words[0:4]
      t2=words[4:8]
      spacing=360.0/npts
      cmapdict["torsions"]= t1+t2
      cmapdict["comment"] = ""
      cmapdict["phi"] = -180.0
      cmapdict["psi"] = -180.0
      cmapdict["energy"]  = list()
      continue

    for t in words:
      cmapdict["energy"].append([cmapdict["phi"],cmapdict["psi"],Decimal(t)])
      cmapdict["psi"] += 15.0
      if(cmapdict["psi"] == 180.0):
        cmapdict["psi"] = -180.0
        cmapdict["phi"] += 15.0    
      cmapnread=cmapnread+1

    if(cmapnread>=cmapntotal):
      FFconv.add_databrick("cmap",cmapdict["energy"])
      cmaptableid=OrderedDict([("cmapid",len(FFconv.databrick["cmap"]))])
      FFconv.add_parameter(vfiles["torsiontorsion"],cmapdict["torsions"],cmaptableid,cmapdict["comment"])
      cmapntotal=0
  
  return


def parse_lj(lines):
  global conv_lj

  for l in lines[1:]:

    (words,comment)=clean_line(l)    
    if(len(words)==0): continue   

    atypes=[ (convatype[a] if convatype.has_key(a) else a) for a in words[0:1]]

    params=OrderedDict([("sigma", Decimal(words[3])), ("epsilon", Decimal(words[2])) ])
    for idx, pname in enumerate(["sigma","epsilon"]): 
      params[pname]*=conv_lj[idx]

    FFconv.add_parameter(vfiles["vdw1_126"],atypes,params,comment)
    
    if(len(words) == 7 ):
      params = OrderedDict([("sigma", Decimal(words[6])), ("epsilon",Decimal(words[5]))])
      for idx, pname in enumerate(["sigma","epsilon"]):
        params[pname]*=conv_lj[idx]

      FFconv.add_parameter(vfiles["vdw1_14_126"],atypes,params,comment)

  return

def parse_nop(lines):
  return

def parse_NBFIX(lines):
  for l in lines[1:]:
    if(len(l)==0): continue
    words,comment=clean_line(l)
    nwords=len(words)
    if nwords==0: continue
    if(nwords!=4):
      print "Invalid NBFIX statement (nwords=%d): %s" %(nwords,str(words))
    if(not opts.viparr3):
      print "NBFIX can only be used with viparr3 forcefields. aborting"
      sys.exit()
      
    atypes=[ (convatype[a] if convatype.has_key(a) else a) for a in words[0:2]]
      
    params=OrderedDict([("sigma", Decimal(words[3])), ("epsilon", Decimal(words[2])) ])
    for idx, pname in enumerate(["sigma","epsilon"]):
      params[pname]*=conv_nbfix_lj[idx]
      
    FFconv.add_parameter("vdw2",atypes,params,comment)
  return

def parse_mass(lines):
  import elementdata
  global atype2atomicnumber
  for l in lines:
    (words,comment)=clean_line(lines[0])
    nwords=len(words)
    if(nwords==0): continue
    if(nwords!=4 and nwords !=5):
      print "Invalid MASS statement (nwords=%d): %s"%(nwords,str(words))
      return  

    idx=int(words[1])
    atype=words[2]
    if(convatype.has_key(atype)):
      atype=convatype[atype]
    mass  = OrderedDict([("amu",Decimal(words[3]))])
    possibleElem=[]
    if(nwords==5):
      element = words[4][0]
      if(len(words[4])>1): element+=words[4][1:].lower()
      possibleElem.append(element)
    else:
      element=atype[0].upper()
      cmass=float(words[3])
      if(element in elementdata.element2number and abs(elementdata.data[element]["mass"][0]-cmass)<=0.1):
        possibleElem.append(element)
      
      for i in xrange(1,len(atype)):
        element+=atype[1].lower()
        if(element in elementdata.element2number and abs(elementdata.data[element]["mass"][0]-cmass)<=0.1):
          possibleElem.append(element)        
    
    if (len(possibleElem)!=1):
      if(atype != "DUM"):
        print "Couldnt convert MASS statement: %s"%(str(words))
      continue
    element=possibleElem[0]
    atype2atomicnumber[atype]=elementdata.element2number[element]

    comment = ""
    FFconv.add_parameter(vfiles["mass"],[atype],mass,comment)
  return

def parse_residue(lines):
  global skipres
  global presdict
  global _has_cmap
  global currenttfile
  global atype2atomicnumber

  otherres=["+","-"]
  resname=''
  resq=Decimal("0.0")
  aqsum=Decimal("0.0")
  for l in lines:

    (words,comment)=clean_line(l)    
    if(len(words)==0): continue   

    if(words[0][0:4].upper() == "RESI" or words[0][0:4].upper() == "PRES"):
      resname = words[1]
      if(resname in skipres): return
      residue = DesFF.Residue()
      if(words[0][0:4].upper() == "PRES"):
        res_is_patch = True
      else:
        res_is_patch = False

      if(convres.has_key(resname)):
        resname=convres[resname]
      residue.resname=resname
      aqsum=Decimal("0.0")
      if(len(words)>2): 
        resq=Decimal(words[2])
      else:
        print "Warning: No charge found for residue %s. assuming 0.0"%(resname)
        resq=Decimal("0.0")
    elif(words[0][0:4].upper() == "GROU"):
      continue
    elif(words[0][0:4].upper() == "ATOM"):
      aname=fix_aname(words[1],resname)
      atype=words[2].upper()
      if(convatype.has_key(atype)):
        atype=convatype[atype]
      acharge=Decimal(words[3])
      aqsum+=acharge
      residue.add_atom(aname,comment)

      if(aname != 'DUM'):
        try:
          residue.add_atom_atomicnumber(aname,atype2atomicnumber[atype])
        except KeyError:
          print "Error while generating template for residue: "+resname
          print "  Atom type '%s'->'%s' wasnt defined in a MASS statement"%(words[2],atype)
          print "  You are most likely trying to convert an incomplete forcefield."
          print "  Try using '--merge' or additional '-p' options to this conversion utility"
          print "  Currently known atomtypes and atomic numbers:"
          for k,v in atype2atomicnumber.iteritems(): print "    %s: %d"%(k,v)
          sys.exit(1)
          
      residue.add_atom_types(aname,[atype])
      residue.add_atom_charge(aname,acharge)
    elif(words[0][0:4].upper() == "BOND" or words[0][0:4].upper() == "DOUB" or words[0][0:4].upper() == "TRIP"):
      words[1:]=fix_aname(words[1:],resname)
      i=1
      while(1):
        if (i>= len(words)):
          break

        if words[i][0] not in otherres and words[i+1][0] not in otherres:
          residue.add_bondedterm("bonds",words[i:i+2],res_is_patch)
        elif(words[i][0] in otherres and words[i+1][0] in otherres):
          raise UserWarning("Cant have bond that isnt part of this residue")
        elif(words[i][0] in otherres):
          noatoms=apply_noatom(words[i:i+2])
          noatoms.reverse() # real atoms 1st is nice
          residue.add_bondedterm("bonds",noatoms,res_is_patch)
          if(not res_is_patch):
            nxt=otherres[(otherres.index(words[i][0])+1)%2]
            noatoms=apply_noatom([words[i][1:],nxt+words[i+1]])
            residue.add_bondedterm("bonds",noatoms,res_is_patch)
        elif(words[i+1][0] in otherres):
          noatoms=apply_noatom(words[i:i+2])
          residue.add_bondedterm("bonds",noatoms,res_is_patch)
          if(not res_is_patch):
            nxt=otherres[(otherres.index(words[i+1][0])+1)%2]
            noatoms=apply_noatom([nxt+words[i],words[i+1][1:]])
            noatoms.reverse() # real atoms 1st is nice 
            residue.add_bondedterm("bonds",noatoms,res_is_patch)
        i=i+2
    elif(words[0][0:4].upper() in ["ANGL",'THET']):
      #These can be automatically generated
      pass
      #words[1:]=fix_aname(words[1:],resname)
      #words[1:]=apply_noatom(words[1:])
      #i=1
      #while(1):
      #  if (i>= len(words)):
      #    break
      #  residue.add_bondedterm("angles",words[i:i+3],res_is_patch)
      #  i=i+3    
    elif(words[0][0:4].upper() == "IMPR" or words[0][0:4].upper() == "IMPH"):
      words[1:]=fix_aname(words[1:],resname)
      words[1:]=apply_noatom(words[1:])
      i=1
      while(1):
        if (i>= len(words)):
          break
        residue.add_bondedterm("impropers",words[i:i+4],res_is_patch)
        i=i+4
    elif(words[0][0:4].upper() == "CMAP"):
      words[1:]=fix_aname(words[1:],resname)
      words[1:]=apply_noatom(words[1:])
      residue.add_bondedterm("cmap",words[1:9],res_is_patch)
      _has_cmap=True
    elif(words[0][0:4].upper() in ["DONO","ACCE", "BILD","PATC", "DIHE"] ):
      continue
    elif(words[0][0:2].upper() == "IC"):      
      continue
    elif(words[0][0:4].upper() == "DELE"):
      if not res_is_patch: raise UserWarning("DELETE directive inside of non-patch")
      if (words[1][0:4].upper() == "ATOM" and len(words)==3):
        words[2:]=fix_aname(words[2:],resname)
        residue.add_bondedterm("delete",[words[2]],res_is_patch)
      elif(words[1][0:4].upper() == "IMPR" and len(words)==6):
        words[2:6]=fix_aname(words[2:6],resname)
        residue.add_bondedterm("delete_impr",words[2:6],res_is_patch)        
      elif(words[1][0:4].upper() != "ACCE"):
        s="Unknown DELETE directive "+l.strip()
        raise UserWarning(s)
    else:
      raise UserWarning("I dont know what to do with: "+str(words))

  if (resname!=''):
    if(res_is_patch):
      presdict[resname]=residue
    else:
      assert aqsum==resq
      FFconv.add_template(currenttfile,resname,residue)
      # print "%s %s" %(currenttfile,FFconv.templates_order["templates."+currenttfile])

  return

def do_patch(pres,residue,del_loc):
  # delete comes first
  if(pres.bondedterms.has_key("delete")):
    atoms = pres.bondedterms["delete"]
    del(pres.bondedterms["delete"])
    for a in atoms:
      rmlist=residue.remove_atom(a[0])
      for type in rmlist:
        for i in rmlist[type]:
          if i<del_loc[type]:
            del_loc[type] = i

  if(pres.bondedterms.has_key("delete_impr")):
    impr = pres.bondedterms["delete_impr"]
    del(pres.bondedterms["delete_impr"])
    for i in impr: residue.remove_bondedterm("impropers",i)
            
  # atom addition comes next (special case of atom replacement)
  for at in pres.atoms:
    found=False
    for i,t in enumerate(residue.atoms):
      if at[0] == t[0]:
        #rint "Rep: "+resname+str(t)+" with "+presname+str(a)
        found=True
        residue.atoms[i]=at
    if (not found):
      # print "Add: "+str(at)+" at "+str(del_loc["atoms"])+" ( "+str(len(residue.atoms))+" )"
      residue.atoms.insert(del_loc["atoms"],at)
      del_loc["atoms"]=del_loc["atoms"]+1

  for tp in pres.bondedterms:
    for i in pres.bondedterms[tp]:
      # print "Add: "+str(i)+" at "+str(del_loc[tp])+" ( "+str(len(residue.bondedterms[tp]))+" )"
      if(tp not in del_loc):
        del_loc[tp]=0
      if(tp not in residue.bondedterms):
        residue.bondedterms[tp]=list()
      residue.bondedterms[tp].insert(del_loc[tp],i)
      del_loc[tp]=del_loc[tp]+1


def apply_pres(tfname, prefix, pname, resname, presname):
  from copy import deepcopy
  global presdict

  lendel={"atoms":     1,
          "bonds":     2,
          "angles":    3,
          "propers"  : 4,
          "impropers": 4,
          "cmap":      8,
          "exclusions":2}
  
  res=FFconv.find_template(tfname,resname)
  if(res == None or not presdict.has_key(pname) or FFconv.find_template(tfname,presname) != None):
    return

  del_loc=dict()
  
  residue=deepcopy(res)
  pres   =deepcopy(presdict[pname])
  #print pres.pprint_residue()

  # clean up patch 
  if(pname=="CT1"):
    atom0=residue.atoms[residue.find_atom("CA")[0]]
    pres.atoms[pres.find_atom("CA")[0]][2]=atom0[2]+Decimal('0.10')
    pres.remove_atom("HN")
    pres.remove_atom("N")
    pres.remove_atom("HA")
  elif(pname=="NNEU"):
    atom0=residue.atoms[residue.find_atom("CA")[0]]
    pidx=pres.find_atom("CA")[0]
    pres.atoms[pidx][2]=atom0[2]+Decimal('0.12')
    if(resname=="GLY"):
      pres.atoms[pidx][3]=atom0[3]
    pres.remove_atom("HA")

  # clean up residue
  if(prefix=='N'):
    residue.remove_atom(bsplit["-C"])
    del_loc["atoms"]=1
    
    for type in residue.bondedterms:
      del_loc[type]=1
  elif(prefix=='C'):
    residue.remove_atom(bsplit["+N"])
    del_loc["atoms"]=len(residue.atoms)
    
    for type in residue.bondedterms:
      del_loc[type]=len(residue.bondedterms[type])
  elif(prefix=='5'):
    residue.remove_atom(bsplit["-O3'"])
    del_loc["atoms"]=1

    for type in residue.bondedterms:
      del_loc[type]=1
  elif(prefix=='3'):
    residue.remove_atom(bsplit["+P"])
    del_loc["atoms"]=len(residue.atoms)

    for type in residue.bondedterms:
      del_loc[type]=len(residue.bondedterms[type])
  else:
    del_loc["atoms"]=len(residue.atoms)
    for type in residue.bondedterms:
      del_loc[type]=len(residue.bondedterms[type])    

  do_patch(pres,residue,del_loc)

  # additional cleanups
  if(pname=="DISU"):
    if("+SG" not in bsplit):
      bsplit["+SG"]="$%d"%(len(bsplit)+1)
    residue.remove_atom("BAD1")
    residue.remove_atom("BAD2")
    residue.remove_atom("BAD3") 
    residue.add_bondedterm("bonds",["SG", bsplit["+SG"]])


  #add new residue
  residue.resname=presname
  FFconv.add_template(tfname,presname,residue)
  return

def apply_pres2(tfname, resname1, resname2, pname, rsuffix):
  from copy import deepcopy  
  res1=FFconv.find_template(tfname,resname1)
  res2=FFconv.find_template(tfname,resname2)  

  if(res1 == None or res2 == None or not presdict.has_key(pname)):
    return None

  residue1=deepcopy(res1)
  residue2=deepcopy(res2)  
  pres   =deepcopy(presdict[pname])

  if(len(rsuffix) != 2 ): raise UserWarning("bad rename_suffix")
  if(rsuffix != ("","")):
    if(rsuffix[0]!=""):
      for atom in residue1.atoms: residue1.rename_atom(atom[0],atom[0]+rsuffix[0])
    if(rsuffix[1]!=""):
      for atom in residue2.atoms: residue2.rename_atom(atom[0],atom[0]+rsuffix[1])
    for atom in pres.atoms+pres.bondedterms.get("delete",[]):
      if(atom[0][0]=="1"):
        pres.rename_atom(atom[0],atom[0][1:]+rsuffix[0])
      elif(atom[0][0]=="2"):
        pres.rename_atom(atom[0],atom[0][1:]+rsuffix[1])
    residue1.merge(residue2)

  del_loc=dict()
  del_loc["atoms"]=len(residue1.atoms)
  for type in residue1.bondedterms:
    del_loc[type]=len(residue1.bondedterms[type])    
  do_patch(pres,residue1,del_loc)

  return residue1

def fixup_amino_acid_patches(tfname):
  from copy import deepcopy
  global presdict
  
  pconv={"ACE":"ACE","CT3":"NMA"}

  for key in pconv:
    if(key not in presdict): continue
    pname=key
    rname=pconv[key]
    pres   =deepcopy(presdict[pname])
    
    if(pres.bondedterms.has_key("delete")):
      raise UserWarning("Cant convert a patch to a residue with atoms to delete: "+pname+" -> "+rname)

    # clear cmap (its handled by the main aa residues)
    if("cmap" in pres.bondedterms):
      del(pres.bondedterms["cmap"])
    
    if(rname=='ACE'):
      # This is handled by main aa residues
      pres.remove_bondedterm("impropers",["N", "CY", "CA", "HN"])
      pres.remove_bondedterm("bonds",["CY","N"]) 
      pres.add_bondedterm("bonds",["CY",bsplit['+N']])
      pres.rename_atom("N",bsplit['+N'])
      
    if(rname=='NMA'):
      # This is handled by main aa residues
      pres.remove_atom("C")
      pres.remove_atom("O")      
      
      pres.add_bondedterm("bonds",["NT",bsplit['-C']])
      pres.add_bondedterm("impropers",["NT",bsplit['-C'],"CAT","HNT"])

    pres.resname=rname
    del(presdict[pname])
    FFconv.add_template(tfname,rname,pres)


def patch_amino_acids():
  global currenttfile
  #This handles ACE/ACP (same patch) and CT3
  fixup_amino_acid_patches(currenttfile)
  
  apply_pres(currenttfile,"", "DISU", "CYS", "CYX" )

  # Protonation States (maestro format):
  # GLH = GLU + H  [GLUP]
  # ASH = ASP + H  [ASPP]
  # LYN = LYS - H  [LSN]
  # ARN = ARG - H  [????]
  
  apply_pres(currenttfile,"", "ASPP", "ASP", "ASH")
  apply_pres(currenttfile,"", "GLUP", "GLU", "GLH")
  apply_pres(currenttfile,"", "LSN",  "LYS", "LYN" )  

  #apply NTER and CTER pres to aa residues:
  aa=['ALA','ARG','ASN','ASP','CYS',
      'GLN','GLU','GLY','HIS','HIE',
      'HIP','ILE','LEU','LYS','MET',
      'PHE','PRO','SER','THR','TRP',
      'TYR','VAL',
      'ASH','GLH','LYN','CYX',
      'CYSF', 'CYSP', 'LYSR']

  for r in aa:
    if(r == 'GLY'):
      apply_pres(currenttfile,"N","GLYP", r,"N"+r)
    elif(r=='PRO'):
      apply_pres(currenttfile,"N","PROP", r,"N"+r) 
    else:
      apply_pres(currenttfile,"N","NTER", r,"N"+r)
    apply_pres(currenttfile,"C","CTER", "N"+r, "NC"+r)

  for tmp in ["CTER","CT1","CT2","CNEU"]:
    for r in aa:
      apply_pres(currenttfile,"C",tmp, r,r+"_"+tmp)

  for r in aa:
    if(r=="PRO"): continue
    apply_pres(currenttfile,"N","NNEU", r,"NNEU_"+r)

  # phosphorylated amino acid compounds
  apply_pres(currenttfile,"", "TP1", "TYR", "MPTYR")
  apply_pres(currenttfile,"", "TP2", "TYR", "DPTYR")
  apply_pres(currenttfile,"", "SP1", "SER", "MPSER")
  apply_pres(currenttfile,"", "SP2", "SER", "DPSER")  
  apply_pres(currenttfile,"", "THP1", "THR", "MPTHR")
  apply_pres(currenttfile,"", "THP2", "THR", "DPTHR")
  

  return


def patch_nucleic_acids():
  global currenttfile
  pyrimidines=["CYT","URA","THY"]
  purines=["ADE","GUA"]

  na=[]
  for r in pyrimidines:
    apply_pres(currenttfile,"","DEO1", r,"D"+r) # oldstyle
    apply_pres(currenttfile,"","DEOX", r,"D"+r) # newstyle
    na.extend([r,"D"+r])
    
  for r in purines:
    apply_pres(currenttfile,"","DEO2", r,"D"+r) # oldstyle
    apply_pres(currenttfile,"","DEOX", r,"D"+r) # newstyle
    na.extend([r,"D"+r])
    
  fiveprime =["5TER","5MET","5PHO","5POM","5DP"]
  threeprime=["3TER","3PO3","3PHO","3POM"]
  for r in na:
    for p5 in fiveprime:
      name=p5+"_"+r
      apply_pres(currenttfile,"5",p5, r,name)
      for p3 in threeprime:
        apply_pres(currenttfile,"3",p3, name, name+"_"+p3)
    for p3 in threeprime:
      apply_pres(currenttfile,"3",p3, r,r+"_"+p3)

def patch_lipids():
  global currenttfile

  def build_lipid(name,res1,res2,res3,patch12,patch13):
    if FFconv.find_template(currenttfile,name) is not None: return
    tmpres=apply_pres2(currenttfile, res1, res2, patch12, ("","_A"))
    if(tmpres is not None): 
      tmpres.resname=name
      FFconv.add_template(currenttfile,name,tmpres)  
      tmpres=apply_pres2(currenttfile, name, res3, patch13, ("","_B"))
      if(tmpres is not None):
        FFconv.replace_template(currenttfile,name,tmpres)
      else:
        FFconv.remove_template(currenttfile,name)

  # build DPPC from PCGL,PALM,PALM,EST1,EST2
  build_lipid("DPPC", "PCGL", "PALM","PALM","EST1","EST2")    
  # build DOPC from PCGL,OLEO,OLEO,EST1,EST2
  build_lipid("DOPC", "PCGL", "OLEO","OLEO","EST1","EST2")    
  # build DSPC from PCGL,STEA,STEA,EST1,EST2
  build_lipid("DSPC", "PCGL", "STEA","STEA","EST1","EST2")  
  # build SDPC from PCGL,STEA,DHA,EST1,EST2
  build_lipid("SDPC", "PCGL", "STEA","DHA","EST1","EST2")


  
def parse_decl(lines):
  for l in lines:
    (words,comment)=clean_line(lines[0])
    if(len(words)==0): continue
    if(len(words)!=2):
      print "Invalid DECL statement (2): %s"%(str(words))
      continue
    decl=words[1]
    if(decl in bsplit): continue
    i=len(bsplit)
    bsplit[decl]="$%d"%(i+1)

def loadparam_and_template(filename, paramonly):
  ffsections = {
    # param sections
    ("BOND", 0) : parse_stretch,
    ("ANGL", 0) : parse_angle,
    ("THET", 0) : parse_angle,
    ("DIHE", 0) : parse_proper,
    ("PHI" , 0) : parse_proper,
    ("IMPR", 0) : parse_improper,
    ("IMPH", 0) : parse_improper,
    ("CMAP", 0) : parse_cmap,
    ("NONB", 1) : parse_lj,
    ("NBON", 1) : parse_lj,
    ("NONB", 0) : parse_lj,
    ("NBON", 0) : parse_lj,
    ("HBON", 1) : parse_nop,
    ("NBFI", 0) : parse_NBFIX,
    # template sections
    ("DECL", 1) : parse_decl,
    ("MASS", 1) : parse_mass,
    ("RESI", 1) : parse_residue,
    ("PRES", 1) : parse_residue,    
    ("DEFA", 1) : parse_nop,
    ("AUTO", 1) : parse_nop,
    ("END",  0) : parse_nop
    }

  file = open(filename,'r')
  pdata=file.readlines()
  file.close()

  sections=[]
  sstop=[]
  for idx,l in enumerate(pdata):
    if len(l)==0: continue
    # strip comments and equals
    l=l.replace('=',' ').replace('!',' ! ')
    pdata[idx]=l
    words,comment=clean_line(l)
    if(len(words)==0): continue
    name=words[0][0:4].upper()
    if (len(words)>1):
      key=(name,1)
    else:
      key=(name,0)
    if ffsections.has_key(key):
      sstop.append(idx)
      while(pdata[idx].split()[-1] == '-'): idx+=1
      sections.append((key,idx))
  sstop.pop(0)
  sstop.append(len(pdata))
  
  for (sn,sb),se in zip(sections,sstop):
    if paramonly:
      if (sn[0]!="RESI" and sn[0]!="PRES"): 
        ffsections[sn](pdata[sb:se])
    else:
      if (sn[0]=="RESI" or sn[0]=="PRES"):
        ffsections[sn](pdata[sb:se])


def get_charmm_rules(viparr3):

  rules= OrderedDict([
    ("info"          , ["info1","info2"]),
    ("vdw_func"      , "LJ12_6_sig_epsilon"),  
    ("vdw_comb_rule" , "ARITHMETIC/GEOMETRIC"),
    # We need to use exclusion=4 because some of the
    # 1-4 vdw interactions use modified parameters
    ("exclusions"  , 4),  
    ("es_scale"    , [ 0.0, 0.0, 1.0 ]),
    ("lj_scale"    , [ 0.0, 0.0, 1.0 ]),    
    ("plugins"     , [
    ["bonds",      0 ],
    ["angles",     0 ],
    ["ureybradley",0 ],
    ["propers",    0 ],
    ["impropers",  0 ],
    ["vdw1",       0 ],
    ["vdw1_14",    0 ],
    ["exclusions", None ],
    ["pairs_es",   None ],
    ["cmap",       None ],
    ["torsion_torsion", 0 ],
    ["mass",       0 ],
    ["atoms",      None ]
    ])
    ])
  if(viparr3):
    rules["plugins"] = ["bonds", "angles", "ureybradley", "propers", "impropers", "vdw1", "exclusions", "cmap", "mass"]
    rules["nbfix_identifier"]="charmm"
  return rules

######################################################################
if __name__ == '__main__':
  import shutil

  usage = '''
  %prog [options] outdir

Description:
  ff_charmm_to_viparr is a force field conversion program.
  -t and -p are used to specify force field topology and parameter files
  respectivly; the order of topology and parameter files are important: 
  If there are conflicts, earlier topologies/parameters take precedence over later.

Simple example:
  %prog [-n template_group_name] -t top_all27_prot_lipid.rtf -p par_all27_prot_lipid.prm charmm27

'''

  templatefiles=OrderedDict()
  def parser_add_template_file(option, opt_str, value, parser):
    templatefiles.setdefault(parser.values.tname[-1],[]).append(value)


  opt = optparse.OptionParser(usage=usage)
  opt.add_option('-n', type='string', dest='tname', action='append',default=[''],
                 help='''Group subsequent templates within this template file.
                 Multiple template groups can be specified, each with its own -n''')
  opt.add_option('-t', type='string', action='callback', dest='tfile',
                 callback=parser_add_template_file,
                 help='''template file to parse; several can be listed,
                 each preceded by its own -t''')  
  opt.add_option('-p', action='append', dest='parameters', default=[],
                 help='''parameter file to parse; several can be listed,
                 each preceded by its own -p''')
  opt.add_option('-r', '--remove-cmap',action='store_true', dest='remove_cmap', default=False,
                 help='''force removal of cmap terms''')
  opt.add_option('-s', '--skip-res',action='append', dest='skipres', default=[],
                 help='''skip residue with this name''')
  opt.add_option('--merge',type='string', dest='merge', default='',
                 help='''base viparr directory for merging''')
  opt.add_option('--viparr3',action='store_true', dest='viparr3', default=False,
                 help='''output files in viparr3 format''')


  opts, args = opt.parse_args()
  if len(args) != 1:
    opt.error('incorrect number of arguments: %d  %s'%(len(args),str(args)))

  print (opts.tname,templatefiles,opts.parameters)
  skipres+=opts.skipres

  if(opts.viparr3):
    DesFF.newparameterstyle=True
  else:
    DesFF.inlineParameterComments=True
  
  vfiles=DesFF.set_viparr_files(opts.viparr3)
  FFconv = DesFF.FF()

  if(opts.merge != ''):
    preload_merge(opts.merge)
    DesFF.inlineTemplateComments=True


  for pload in [True, False]:
    if(pload):
      for f in opts.parameters:
        loadparam_and_template(f,pload)    

    currenttfile=''
    for k in templatefiles:
      currenttfile=k
      for t in templatefiles[k]:
        loadparam_and_template(t,pload)
      if(not pload):
        patch_amino_acids()
        patch_nucleic_acids()
        patch_lipids()

  FFconv.rules=get_charmm_rules(opts.viparr3)
  if(not _has_cmap or opts.remove_cmap):
    if(_has_cmap):
      print "User requested cmap removal"
    else: 
      print "No CMAP terms found"
    if(opts.viparr3):
      FFconv.rules["plugins"].remove("cmap")
    else:
      FFconv.rules["plugins"].remove(["cmap",       None ])
      FFconv.rules["plugins"].remove(["torsion_torsion", 0])
    if("cmap" in FFconv.databrick): del(FFconv.databrick["cmap"])
    if("torsiontorsion_cmap" in FFconv.parameters): del(FFconv.parameters["torsiontorsion_cmap"])
    for tf in FFconv.templates:
      for res in FFconv.templates[tf]: 
        if "cmap" in FFconv.templates[tf][res].bondedterms:
          del(FFconv.templates[tf][res].bondedterms["cmap"])
          FFconv.templates[tf][res].bondedtermkeys.remove("cmap")
  if(opts.viparr3 and "vdw2" in FFconv.parameters):      
    FFconv.rules["plugins"].append("vdw2")

  dirname=args[0]
  if(os.path.exists(dirname)):
    shutil.rmtree(dirname)
  os.mkdir(dirname)
  FFconv.write_forcefield(dirname)

  if(opts.merge): os.unlink(os.path.join(dirname,"rules"))
