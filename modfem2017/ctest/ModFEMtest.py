### ModFEM test runner - common for all problem modules
# Define some variables:
import sys
import subprocess
import os
import glob
import shutil
import zipfile

# function for comparing values
def TestOutputForValue(config_file, file_to_search, text_to_find):
  #wczytanie configa
  # k='conf.txt'
  name=open(config_file,'r')

  for i in range(9):
    w=name.readline()
    while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  p2=w.replace('\n','')
  p3=float(p2)
  name.close()


  #plik1='./output.txt'
  plik1 = file_to_search
  if os.path.isfile(plik1)==False :
      plik1=plik1+'_1'

  name=open(plik1,'r')
  p1=name.read()

  t1=p1.find (text_to_find)
  t2=p1.find ('\n')
  p1=p1[t1:t2]

  p1=p1.replace(' ','')
  t=p1.find ("=")
  p1=p1[t+1:]
  pn=float(p1)

  if pn <p3 :
      return ("Test passed. L2= "+p1+"<"+p2)
  else :
      return ("Fail. Too big Error L2= "+p1+">"+p2)

default_testing_config_file = 'conf.txt'
default_testing_config_file_mpi = 'conf.txt'

def GetCorrectVersionOfConfigFile() :
  if( ModFEMtestRunner.is_mpi ):
    #print("MPI")
    return default_testing_config_file_mpi
  else :
    if( ModFEMtestRunner.config_file_name != ""):
      #print("Z configiem: "+ModFEMtestRunner.config_file_name)
      return ModFEMtestRunner.config_file_name
    else :
      #print("Bez configa: "+default_testing_config_file)
      return default_testing_config_file

default_testing_output_file = 'output.txt'
default_testing_output_file_mpi = 'output.txt_1'

def GetCorrectVersionOfOutputFile() :
  if( ModFEMtestRunner.is_mpi ):
    return default_testing_output_file_mpi
  else :
    return default_testing_output_file

def TestL2(  ) :
  # testing for L2 Error
  #wczytanie configa
  k=GetCorrectVersionOfConfigFile()
  conf_name=open(k,'r')

  for i in range(8):
    w=conf_name.readline()
    while w.startswith('\n') or w.startswith('#'):
      w=conf_name.readline()
  p2=w.replace('\n','')
  conf_name.close()


  plik1=GetCorrectVersionOfOutputFile()
  if os.path.isfile(plik1)==False :
    plik1='./output.txt_1'

  name=open(plik1,'r')
  p1=name.read()
  t1=p1.find (" L2 norm of error ")
  p1=p1[t1:]
  t2=p1.find ('\n')
  p1=p1[:t2]

  p1=p1.replace(' ','')
  t=p1.find ("=")
  p1=p1[t+1:]
  pn=float(p1)

  
  try:
    p3=float(p2)
  except ValueError:
    pt=p2.split(" ")
    pmin=float(pt[0])
    pmax=float(pt[1])
    if pn <pmax and pn> pmin:
      msg = "Test passed. L2 is in interval: "+p1+" = ["+p2+"]"
    else :
      msg = "Fail. L2 isnt in interval: "+p1+" != ["+p2+"]"
  else:
    if pn <p3 :
      msg = "Test passed. L2= "+p1+"<"+p2
    else :
      msg = "Fail. Too big Error L2= "+p1+">"+p2
  return msg


def TestH1(  ) :
  # testing for H1 error
  #wczytanie configa
  k=GetCorrectVersionOfConfigFile()
  #print(ModFEMtestRunner.config_file_name)
  conf_name=open(k,'r')

  for i in range(7):
    w=conf_name.readline()
    while w.startswith('\n') or w.startswith('#'):
      w=conf_name.readline()
  p2=w.replace('\n','')
  conf_name.close()

  plik1=GetCorrectVersionOfOutputFile()
  if os.path.isfile(plik1)==False :
    plik1='./output.txt_1'

  name=open(plik1,'r')
  p1=name.read()

  t=p1.find ("H1 seminorm of error")
  p1=p1[t:]
  t2=p1.find('\n')
  p1=p1[:t2]
  p1=p1.replace(' ','')
  t=p1.find ("=")
  p1=p1[t+1:]
  #p1=p1.split('\n', 1)[0]
  pn=float(p1)
  

  try:
    p3=float(p2)
  except ValueError:
    pt=p2.split(" ")
    pmin=float(pt[0])
    pmax=float(pt[1])
    if pn <pmax and pn> pmin:
      msg = "Test passed. H1 is in interval: "+p1+" = ["+p2+"]"
    else :
      msg = "Fail. H1 isnt in interval: "+p1+" != ["+p2+"]"
  else:
    if pn <p3 :
      msg = "Test passed. H1= "+p1+"<"+p2
    else :
      msg = "Fail. Too big Error H1= "+p1+">"+p2

  return msg

def TestZZ(  ) :
  #testing for ZZ
  #wczytanie configa
  k=GetCorrectVersionOfConfigFile()
  conf_name=open(k,'r')

  for i in range(9):
    w=conf_name.readline()
    while w.startswith('\n') or w.startswith('#'):
      w=conf_name.readline()
  p2=w.replace('\n','')
  conf_name.close()


  plik1='./zz'
  name=open(plik1,'r')
  p1=name.read()

  p1=p1.replace(' ','')
  t=p1.find ("=")
  p1=p1[t+1:]
  pn=float(p1)

  for i in glob.glob(plik1):
    os.unlink (i)
  
  try:
    p3=float(p2)
  except ValueError:
    pt=p2.split(" ")
    pmin=float(pt[0])
    pmax=float(pt[1])
    if pn <pmax and pn> pmin:
      msg = "Test passed. ZZ is in interval: "+p1+" = ["+p2+"]"
    else :
      msg = "Fail. ZZ isnt in interval: "+p1+" != ["+p2+"]"
  else:
    if pn <p3 :
      msg = "Test passed. ZZ= "+p1+"<"+p2
    else :
      msg = "Fail. Too big Error ZZ= "+p1+">"+p2
  return msg


def TestDiff( ) :
  #wczytanie configa
  k=GetCorrectVersionOfConfigFile()
  name=open(k,'r')
  w=name.readline()

  while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  plik1=w.replace('\n','')
  w=name.readline()
  while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  plik2=w.replace('\n','')
  w=name.readline()
  while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  sc=w.replace('\n','')
  w=name.readline()
  while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  cyfr=str(w.replace('\n',''))
  w=name.readline()
  while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  jump=str(w.replace('\n',''))
  w=name.readline()
  while w.startswith('\n') or w.startswith('#'):
      w=name.readline()
  limit=str(w.replace('\n',''))
  name.close()

  zip_ref = zipfile.ZipFile(plik1+'.zip', 'r')
  zip_ref.extractall('.')
  zip_ref.close()
  zip_ref = zipfile.ZipFile(plik2+'.zip', 'r')
  zip_ref.extractall('./comparison')
  zip_ref.close()

  #wczytanie porownywanych plikow
  name=open(plik1,'r')
  namec=open(plik2,'r')
  p1=name.read()
  p2=namec.read()
  name.close()
  namec.close()

  sum=len(p1.split())


  ile = subprocess.Popen([sys.executable, sc+'fileparser.py', cyfr, sc, plik1, plik2, jump, limit], stdout=subprocess.PIPE).communicate()[0]
  ile=int(ile)
  ile=sum-ile

  #usuniecie plikow *.dmp
  for i in glob.glob('*.dmp'):
      os.unlink (i)
  for i in glob.glob('*.dmp.zip'):
      os.unlink (i)
  for i in glob.glob(plik2):
      os.unlink (i)

  if sum == ile :
      msg = "Test passed. "+str(ile)+" out of "+str(sum)
  else :
      msg = "Failed. "+str(ile)+" out of "+str(sum)

  return msg

# define a class
class ModFEMtestRunner:
  """A class for running tests with ModFEM problem binaries"""
  is_mpi = False
  config_file_name = ""

  def __init__(self, bin, target_dir, ppath, ld_paths, n_proc, config_name):
    self.bin = bin
    self.target_dir = target_dir
    self.ppath = ppath
    self.ld_paths = ld_paths
    self.n_proc = int(n_proc)
    ModFEMtestRunner.config_file_name = config_name
    ModFEMtestRunner.is_mpi = (self.n_proc > 0)
    self.SetEnvVars(self.target_dir,self.ld_paths)
    self.CleanFiles()
    

    #print("\n\n!!! "+config_name+" "+n_proc)
    
  def SetEnvVars(self,target_dir,lib_paths):
    #ustawienie zmiennych srodowiskowych
    en1=''
    en2=''
    if os.environ.has_key('LIBRARY_PATH'):
      en1=os.environ["LIBRARY_PATH"]
    if os.environ.has_key('LD_LIBRARY_PATH'):
      en2=os.environ["LD_LIBRARY_PATH"]
    os.environ["LIBRARY_PATH"]=en1+":"+en2+":"+lib_paths
    os.environ["LD_LIBRARY_PATH"]=en1+":"+en2+":"+lib_paths
    #przejscie do katalogu ze skryptem
    os.chdir(target_dir)

  def CleanFiles(self):
    #usuwanie plikow
    for i in glob.glob('*.dmp'):
      os.unlink (i)
    for i in glob.glob('output.txt'):
      os.unlink (i)

  def Run(self):
    if os.path.isfile(self.ppath+"/"+self.bin)==False :
      self.bin=self.bin+"_d"
      if os.path.isfile(self.ppath+"/"+self.bin)==False :
        print('Fail. Binary file doesnt exist!')
        exit()
    if self.n_proc > 0 :
      #aktualizacja plikow
      shutil.copy("input_interactive.txt_mpi", "input_interactive.txt")
      self.output_file = 'output.txt_1'
      
      #aktualizacja plikow
      if( ModFEMtestRunner.config_file_name != "") :
        k=GetCorrectVersionOfConfigFile()
        conf_name=open(k,'r')

        #wczytanie nazwy pliku problemu
        for i in range(10):
          w=conf_name.readline()
          while w.startswith('\n') or w.startswith('#'):
            w=conf_name.readline()
        problemdat=w.replace('\n','')
        problemdatnew=problemdat[0:problemdat.rfind(".dat")+4]
        
        #wczytanie nazwy pliku slvera
        w=conf_name.readline()
        while w.startswith('\n') or w.startswith('#'):
          w=conf_name.readline()
        solverdat=w.replace('\n','')
        solverdatnew=solverdat[0:solverdat.rfind(".dat")+4]
        
        #wczytanie nazwy pliku input_interactive
        w=conf_name.readline()
        while w.startswith('\n') or w.startswith('#'):
          w=conf_name.readline()
        inputdat=w.replace('\n','')
        
        #print('Fail. Binary file doesnt exist!'+problemdat+" "+solverdat+" "+inputdat)
        #print(solverdatnew+"ZZ"+solverdat+"ZZ"+ModFEMtestRunner.config_file_name)
        conf_name.close()
        
        #shutil.copy(inputdat, "input_interactive.txt")
        
        ############## tymczasowe  ####################
        #if(k != "conf.txt" ) :
        
        shutil.copy(ModFEMtestRunner.config_file_name, "conf.txt")
        #shutil.copy(inputdat, "input_interactive.txt")
        ###############################################
        
        if(problemdat != problemdatnew ) :
          shutil.copy(problemdat, problemdatnew)
        if(solverdat != solverdatnew ) :
          shutil.copy(solverdat, solverdatnew)
        #shutil.copy(inputdat, "input_interactive.txt")
      
      self.retcode = subprocess.Popen(["mpirun","-np",str(self.n_proc),self.ppath+"/"+self.bin,'.'],stdout=subprocess.PIPE).communicate()[0]
      self.output = open(self.output_file,'r').read()
    else :
      #aktualizacja plikow
      if( ModFEMtestRunner.config_file_name != "") :
        k=GetCorrectVersionOfConfigFile()
        conf_name=open(k,'r')

        #wczytanie nazwy pliku problemu
        for i in range(10):
          w=conf_name.readline()
          while w.startswith('\n') or w.startswith('#'):
            w=conf_name.readline()
        problemdat=w.replace('\n','')
        problemdatnew=problemdat[0:problemdat.rfind(".dat")+4]
        
        #wczytanie nazwy pliku slvera
        w=conf_name.readline()
        while w.startswith('\n') or w.startswith('#'):
          w=conf_name.readline()
        solverdat=w.replace('\n','')
        solverdatnew=solverdat[0:solverdat.rfind(".dat")+4]
        
        #wczytanie nazwy pliku input_interactive
        w=conf_name.readline()
        while w.startswith('\n') or w.startswith('#'):
          w=conf_name.readline()
        inputdat=w.replace('\n','')
        
        #print('Fail. Binary file doesnt exist!'+problemdat+" "+solverdat+" "+inputdat)
        #print(solverdatnew+"ZZ"+solverdat+"ZZ"+ModFEMtestRunner.config_file_name)
        conf_name.close()
        
        #shutil.copy(inputdat, "input_interactive.txt")
        
        ############## tymczasowe  ####################
        #if(k != "conf.txt" ) :
        
        shutil.copy(ModFEMtestRunner.config_file_name, "conf.txt")
        #shutil.copy(inputdat, "input_interactive.txt")
        ###############################################
        
        if(problemdat != problemdatnew ) :
          shutil.copy(problemdat, problemdatnew)
        if(solverdat != solverdatnew ) :
          shutil.copy(solverdat, solverdatnew)
        shutil.copy(inputdat, "input_interactive.txt")
      else :
        shutil.copy("input_interactive.txt_run", "input_interactive.txt")
      self.output_file = 'output.txt'
      self.retcode = subprocess.Popen([self.ppath+"/"+self.bin,'.'],stdout=subprocess.PIPE).communicate()[0]
      # zz working only in not distributed env.
      #t=self.retcode.find ("Zienkiewicz-Zhu error estimator =")
      self.output = open(self.output_file,'r').read()
      t=self.output.find ("Zienkiewicz-Zhu error estimator =")
      p1=self.output[t:]
      t=p1.find ("\n")
      p1=p1[:t]
      #print(p1)
      plik = open('zz', 'w')
      plik.write(p1)
      plik.close()

    for i in glob.glob('input_interactive.txt'):
      os.unlink (i)

    zm=int(self.output.count("Total solver times") )

    #print('>>>>>>>>>>>>>>> self.output='+self.output)
    #print('>>>>>>>>>>>>>>> self.retcode='+self.retcode)

    if zm < 1 :
      zm=int(self.output.count("Time for solving the problem") )

    if zm < 1 :
      zm=int(self.retcode.count("Total solver times") )

    if zm < 1 :
      zm=int(self.retcode.count("Time for solving the problem") )
      
    if zm > 0 :
      self.run_msg = 'Test passed.' +self.retcode
    else :
      self.run_msg = 'Fail. ' +self.retcode

    return self.run_msg





