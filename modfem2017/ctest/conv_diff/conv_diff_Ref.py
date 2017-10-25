#!/usr/bin/env python


import sys
import subprocess
import os
import glob
import shutil


arg1=sys.argv[1]
ppath=sys.argv[2]
arg3=sys.argv[3]


#ustawienie zmiennych srodowiskowych
en1=''
en2=''
if os.environ.has_key('LIBRARY_PATH'):
    en1=os.environ["LIBRARY_PATH"]
if os.environ.has_key('LD_LIBRARY_PATH'):
    en2=os.environ["LD_LIBRARY_PATH"]
os.environ["LIBRARY_PATH"]=en1+":"+en2+":"+arg3
os.environ["LD_LIBRARY_PATH"]=en1+":"+en2+":"+arg3


#przejscie do katalogu ze skryptem
os.chdir(arg1)

#print(os.listdir("."))

#usuwanie plikow
for i in glob.glob('error'):
    os.unlink (i)
for i in glob.glob('output.txt'):
    os.unlink (i)

#aktualizacja plikow
shutil.copy("input_interactive.txt_ref", "input_interactive.txt")

k='conf.txt'
name=open(k,'r')

for i in range(7):
    w=name.readline()
    while w.startswith('\n') or w.startswith('#'):
        w=name.readline()
bin=w.replace('\n','')
name.close()
#print(os.getcwd())

retcode = subprocess.Popen([ppath+"/"+bin,'.'], stdout=subprocess.PIPE).communicate()[0]

for i in glob.glob('input_interactive.txt'):
    os.unlink (i)
    
zm=retcode[1:26]

plik = open('output.txt')
try:
	tekst = plik.read()
finally:
	plik.close()


t=tekst.find ("After refinement test.")
p1=tekst[t:]
t=p1.find ("Elements")
p1=p1[t:]

plik = open('comparison/ref.txt')
try:
	tekst = plik.read()
finally:
	plik.close()


if tekst == p1 :
    print('Test passed.'+p1)
else :
    print('Fail. '+p1)
