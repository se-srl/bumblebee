from os.path import basename
import sys, getopt
import csv
import re
from collections import defaultdict
from difflib import SequenceMatcher
from Levenshtein import *
from warnings import warn
import difflib
import time
import gc
import os
import collections

stringVal=""
from tasks import searchFunc

threadFound=False
speciesName=""

def usage(arg) :
  
  sys.exit(2)

def matchFunc(matchVal, myukdict):
  distHM=defaultdict(list)

  for name,species in myukdict.items():
    # print name
    for spec in species:
      #print specCount
      
      matchV=isMatch(matchVal, spec)
      if matchV[0]>-1:
        if not name in distHM[matchV[0]]: distHM[matchV[0]].append(name)
      
      #if matchVal in spec:
        #print name
        # endTime= time.time()
        # speciesName=name
  return distHM



def reverseComplement(val):
  val=val[::-1]
  val=val.replace('A','B')
  val=val.replace('T','A')
  val=val.replace('B','T')
  val=val.replace('C','B')
  val=val.replace('G','C')
  val=val.replace('B','G')
  return val



def isMatch1(a,b):
  dist=distance(str(a), str(b))
  distVal= float("{0:.2f}".format(dist))    
  matchDist=100-((len(a)+distVal-len(b))/len(a)*100)
  if matchDist>50:
    return matchDist
  return -1



def isMatch(a,b,bases):
  i=0
  j=0
  position=0
  count=0
  maxVal=0
  if len(b)>len(a):
    while i<=(len(b)-len(a)):
      substr=b[i:i+len(a)]
      t=0
      j=0
      count=0
      while j<len(a):
        if a[j]==substr[t] or substr[t]=='N' or a[j]=='N':
          
          count+=1
        j+=1
        t+=1
      if len(a) - count <= bases:
        return True
      i+=1
  return False

def preProcess(forwardPrimer, reversePrimer, val):
  start=val.find(forwardPrimer)
  if start >=0:
    stringStart=start+len(forwardPrimer)
    stringEnd=stringStart+len(forwardPrimer)
    
    val=val[stringStart:]


   
   
  else:
    

    start=val.find(reversePrimer)

    if start>=0:
      stringStart=start+len(reversePrimer)
      
      val=val[stringStart:]
      
      val=reverseComplement(val)
    else:
      print("Fuck")
  return val

def findLocation(startList,endList,val):
  listVal=[val]*len(startList)
  res=searchFunc.chunks(zip(startList, endList,listVal),10)()

  searchVals=[]
  for vals in res:
    if vals != "NOT":
      searchVals.append(val)
  return searchVals

def resolveGenus(speciesList):
  genusList=[]
  for species in speciesList:
    genusList.append(species.split(" ")[0])
  return  (len(set(genusList)) == 1)

def findLocation(val, motifdict, keys,blastdb):
  potentialSpeciesList=[]
  for key in motifdict.keys():
    start=key.split(',')[0]
    # if start=="TATACTCCTGAATACGAAACCAAAGATACTGATATCTTGGCAGCA" or start=="TATACTCCTGAATACGAAACCAAAGATACTGATATCTTGGCAGCA" or start =="CAAGTATGGTCGTCCCCTGTTGGGATGTACTATTAAACCTAAATT" or start=="TATACTCCTGAATATGAAACCAAGGATACTGATATCTTGGCAGCA" or start=="CAAGTATGGTCGTCCTCTGTTGGGATGTACTATTAAACCTAAATT": 
    #   continue
    if val.find(start)>=0 :
      potentialSpeciesList.extend(motifdict[key])
    else:
      if val.find(reverseComplement(start)) >= 0:
        potentialSpeciesList.extend(motifdict[key])
      
  #if it doesn't resolve to genus
  if os.path.exists("OutFile"):
    os.remove("OutFile")

  fm="1"
  resolved = resolveGenus(potentialSpeciesList)
  dictMatch=defaultdict(list)
  if not resolved:
    for species in potentialSpeciesList:
      bc_strings = []
      bc_names = []
      for string in keys[species]:
        if len(string)>0:
          bc_strings = []
          bc_names = []
          speciesL=[]
          speciesN=[]
          string_names=""
          string_number=""
          if os.path.exists("query"):
            os.remove("query")
          with open('query', 'w') as the_file:
            print(string)
            the_file.write(string[45:-45])
            print("Querying")
            print(string[45:-45])
          if os.path.exists("OutFile"):
            os.remove("OutFile")
          import subprocess
          bashCommand1 = "blastn -query query -db "+ blastdb +" -out OutFile"
          subprocess.run(bashCommand1,shell=True)
          bc = " cat OutFile  | grep \">\" | awk NR==1 | awk '{print $1,$4}'  | cut -d \">\" -f 2 | tr '_' ' '"
          bc_strings = " sed -ne '/Value/,$ p' OutFile |  sed -n '/>/q;p' |sed '2d' | awk '{print $1}'"
          bc_numbers = " sed -ne '/Value/,$ p' OutFile |  sed -n '/>/q;p' |sed '2d' | awk '{print $(NF-1)}'"
          bc_accuracy = "cat OutFile | grep Identities | awk '{print $4}'|grep -o '[0-9]\+'"

          output=subprocess.run(bc_strings, shell=True,stdout=subprocess.PIPE)
          string_names= output.stdout.decode('utf-8')
          output=subprocess.run(bc_numbers, shell=True,stdout=subprocess.PIPE)

          string_number = output.stdout.decode('utf-8')

          output = subprocess.run(bc_accuracy, shell=True,stdout=subprocess.PIPE)
          string_accuracy = output.stdout.decode('utf-8')
          speciesL = string_names.split('\n')[1:]
          speciesN = [int(float((item))) for item in string_number.split('\n')[1:-2]]
          accuracyN = [int(item) for item in string_accuracy.split('\n')[:-1]]
          if len(speciesN) > 0:
            maxVal= max(speciesN)
            
            countVals = speciesN.count(maxVal)
            if (countVals == 1 and accuracyN[0]>=97):
              return speciesL[0].replace('_',' ')+"!1"
            else:
              if (accuracyN[0]<97):
                break
              else:
                print("NOT")
                 
                retval=(';'.join(speciesL[0:countVals])).replace('_',' ')+"!1"
                return retval
          else:
            break

          

       
        
    print("using second motif")
    if os.path.exists("OutFile"):
      os.remove("OutFile")
    potentialSpeciesList=[]
    for key in motifdict.keys():
      end=key.split(',')[1]
      if val.find(end)>=0 :
        potentialSpeciesList.extend(motifdict[key])
      else:
        if val.find(reverseComplement(end)) >= 0:
          potentialSpeciesList.extend(motifdict[key])
        

    
    resolved = resolveGenus(potentialSpeciesList)
    dictMatch=defaultdict(list)
    if not resolved:
      for species in potentialSpeciesList:
        bc_strings = []
        bc_names = []
        for string in keys[species]:
          bc_strings = []
          bc_names = []
          string_names=""
          string_number=""
          speciesL=[]
          speciesN=[]
          if os.path.exists("query"):
            os.remove("query")
          if os.path.exists("OutFile"):
            os.remove("OutFile")
          with open('query', 'w') as the_file:
            the_file.write(string[45:-45])
          print("running blast second")
          print(string[45:-45])
          import subprocess

          bashCommand1 = "blastn -query query -db "+ blastdb +" -out OutFile"
          subprocess.run(bashCommand1,shell=True)
          bc = " cat OutFile "
          output=subprocess.run(bc, shell=True,stdout=subprocess.PIPE)
          bc_strings = "cat OutFile | sed -ne '/Value/,$ p' OutFile |  sed -n '/>/q;p' |sed '2d' | awk '{print $1}'"
          bc_numbers = "cat OutFile | sed -ne '/Value/,$ p' OutFile |  sed -n '/>/q;p' |sed '2d' | awk '{print $(NF-1)}'"
          bc_accuracy = "cat OutFile | grep Identities | awk '{print $4}'|grep -o '[0-9]\+'"

          output=subprocess.run(bc_strings, shell=True,stdout=subprocess.PIPE)
          string_names= output.stdout.decode('utf-8')
          output=subprocess.run(bc_numbers, shell=True,stdout=subprocess.PIPE)
          output = subprocess.run(bc_accuracy, shell=True,stdout=subprocess.PIPE)
          string_accuracy = output.stdout.decode('utf-8')
          string_number = output.stdout.decode('utf-8')
          speciesL = string_names.split('\n')[1:]
          speciesN = [int(item) for item in string_number.split('\n')[1:-2]]
          accuracyN = [int(item) for item in string_accuracy.split('\n')[:-1]]
          maxVal= max(speciesN)
          countVals = speciesN.count(maxVal)
          
          if (countVals == 1 and accuracyN[0]>=97):
            return speciesL[0].replace('_',' ')+"!2"
          else:
            if (accuracyN[0]<97):
              return " !2"
            else:
              print("NOT")
              return (';'.join(speciesL[0:countVals])).replace('_',' ')+"!2"
    else:
      return ';'.join(set(potentialSpeciesList))+"!2"      
  
  print("resolved bitch")
  return ';'.join(set(potentialSpeciesList))+"!1" 

def findLocationLong(startList,endList,val):
  searchVals=[]
  maxVal=0
  positionLeft=0
  positionRight=0
  found=False
  for i in range(0,len(startList)):
    foundLeft=isMatch(startList[i],val)
    foundRight=isMatch(endList[i],val)
    if foundLeft[1]>0 and foundRight[1]>0 and foundRight[1]>foundLeft[1]:
      if min(foundLeft[0],foundRight[0]) > maxVal:
        maxVal=min(foundLeft[0],foundRight[0])
        positionLeft=foundLeft[1]
        positionRight=foundRight[1]
        found=True
  if found:
    print(maxVal)
    print(positionLeft)
    print(positionRight)
    return val[positionLeft:positionRight]   
  return None

def fullMatchSpecies(val, speciesDict):

  for species,vals in speciesDict.items():
    for v in vals:
      dist=distance(v, val)
      distVal= float("{0:.2f}".format(dist))    
      matchDist=100-((len(v)+distVal-len(val))/len(v)*100)
      matchDist2=100-((len(val)+distVal-len(v))/len(val)*100)
      if matchDist>80 or matchDist2>80:
        print(val)
        print(species)
        return species

  return ''

def main(argv) :  

  fastafile='input/RhosReference-fasta.csv'
  ukfastafile='input/ukspeciesrbcl.csv'
  pollenfile='input/out.csv'
  forwardPrimer='ATGTCACCACAAACAGAGACTAAAGC'
  reversePrimer='AGTCCACCGCGTAGACATTCAT'
  nucleotideFile='input/nucleotide.txt'
  matchList=[]
  
  ukstartList=[]
  ukendList=[]
  pairList=[]
  startList=[]
  endList=[]

  motifdict=defaultdict(list)
  startmotifdict=defaultdict(list)
  endmotifdict=defaultdict(list)
  speciesDict = defaultdict(list)

  with open(fastafile) as csvfile:
    csv_reader= csv.DictReader(csvfile)
    for row in csv_reader:
      species=row['Species']
      rbcl=row['RBCL'].replace('-','')
      start=rbcl.strip()[0:45]
      end=rbcl.strip()[-45:]
      startList.append(start)
      endList.append(end)
      if species not in speciesDict:
        speciesDict[species].append(rbcl)
      else:
        if rbcl not in speciesDict[species]:
          speciesDict[species].append(rbcl)
      if species not in motifdict[start+","+end]:
        motifdict[start+","+end].append(species)
      if species not in startmotifdict[start]:
        startmotifdict[start].append(species)
      if species not in endmotifdict[end]:
        endmotifdict[end].append(species)
  
  for key in motifdict:
    with open('motif-all.csv', 'a') as the_file:
      the_file.write(key+"," +  ','.join(motifdict[key])+"\n")
  #reversing the motif dict -- finding out which species is exactly identified by which motif
  reversemotifDict = defaultdict(list)
  for key,value in motifdict.items():
    for species in value:
      if key not in reversemotifDict[species]:
        reversemotifDict[species].append(key)

  #if there is a species for which there is a motif that uniquely identifies it, remove that species from the other motifs
  for species,motifs in reversemotifDict.items():
    if len(motifs) == 1:
      #look for the species name in all the motifs and remove it
      for key, value in motifdict.items():
        if key not in motifs:
          if species in value:
            motifdict[key].remove(species)


  for key in motifdict:
    with open('motif.csv', 'a') as the_file:
      the_file.write(key+"," +  ','.join(motifdict[key])+"\n")
   
  


 
  with open(pollenfile) as csvfile:
    csv_reader= csv.DictReader(csvfile)
    mypollen = defaultdict(list)
    next(csv_reader)
    for row in csv_reader:
      mypollen[row['ID']].append(row['Value'])


  for key in mypollen.keys():
    total_length=len(mypollen[key])
    for val in mypollen[key]:
      if len(val) < 400:
        mypollen[key].remove(val)
    total_removed_perc = len(mypollen[key])*1.0/total_length
   
 
  print(mypollen.keys())
  
  t=0
  for key in list(mypollen.keys())[t:]:
    print(key)
    count=0
    length=0
    for val in mypollen[key]:
      if len(val) >=400:
        length+=1
        print("searching Rhos")
        # remove forward/reverse primer and do reverse complement
        val=preProcess(forwardPrimer,reversePrimer,val)
        #are these exact matches?
        #do an exact match of the species to see if we can find it
        searchVal = findLocation(val, motifdict,speciesDict,"nucl5")
        print(searchVal)
        firstM=searchVal.split("!")[1]
        searchVal = searchVal.split("!")[0].strip()
        if len(searchVal)>0:  
          print("found in Rhos! "+ searchVal)
          count+=1
          with open('bee-results-alluk.csv', 'a') as the_file:
            print(key+"," + val + ", 1 , 0 ," + firstM+','+searchVal+"\n")
            the_file.write(key+"," + val + ", 1 , 0 ," + firstM+','+searchVal+"\n")
            if "Cardamine" in searchVal:
              return
        else:
          print("searching in uk!")
          #search for it in uk
          #create the speciesDict 
          motifdict=defaultdict(list)
          startmotifdict=defaultdict(list)
          endmotifdict=defaultdict(list)
          speciesDict = defaultdict(list)

          with open(ukfastafile) as csvfile:
            csv_reader= csv.DictReader(csvfile)
            for row in csv_reader:
              species=row['Species']
              rbcl=row['RBCL'].replace('-','')
              start=rbcl.strip()[0:45]
              end=rbcl.strip()[-45:]
              startList.append(start)
              endList.append(end)
              if species not in speciesDict:
                speciesDict[species].append(rbcl)
              else:
                if rbcl not in speciesDict[species]:
                  speciesDict[species].append(rbcl)
              if species not in motifdict[start+","+end]:
                motifdict[start+","+end].append(species)
              if species not in startmotifdict[start]:
                startmotifdict[start].append(species)
              if species not in endmotifdict[end]:
                endmotifdict[end].append(species)
          
          for key1 in motifdict:
            with open('motif-all-uk.csv', 'a') as the_file:
              the_file.write(key1+"," +  ','.join(motifdict[key1])+"\n")
          #reversing the motif dict -- finding out which species is exactly identified by which motif
          reversemotifDict = defaultdict(list)
          for key1,value in motifdict.items():
            for species in value:
              if key1 not in reversemotifDict[species]:
                reversemotifDict[species].append(key1)

          #if there is a species for which there is a motif that uniquely identifies it, remove that species from the other motifs
          for species,motifs in reversemotifDict.items():
            if len(motifs) == 1:
              #look for the species name in all the motifs and remove it
              for key3, value in motifdict.items():
                if key3 not in motifs:
                  if species in value:
                    motifdict[key3].remove(species)


          for key2 in motifdict:
            with open('motif-uk.csv', 'a') as the_file:
              the_file.write(key2+"," +  ','.join(motifdict[key2])+"\n")

          searchVal = findLocation(val, motifdict,speciesDict,"uk")
          firstM=searchVal.split("!")[1]
          searchVal = searchVal.split("!")[0].strip()
          #create motifDict
          if len(searchVal)>0:  
            print("found in uk!"+searchVal)
            count+=1
            with open('bee-results-alluk.csv', 'a') as the_file:
              the_file.write(key+"," + val + ", 0 , 1 ," + firstM+','+searchVal+"\n")
          else:
            with open('bee-results-alluk.csv', 'a') as the_file:
              the_file.write(key+"," + val + ", 0 , 0 , +"+firstM+",-\n")

    if length>0:
      perc=count/length*100
    else:
      perc=-1
    with open('perc.csv','a') as the_file:
      the_file.write(key+","+str(perc)+"\n")
    print(str(perc))
  return
  
  
  

if __name__ == "__main__" : 
  main(sys.argv[1:])


