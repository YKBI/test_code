# -*- coding: utf-8 -*-
"""
Created on Thu May  7 11:16:19 2020

@author: syntekabio
"""
import os
#lists
rec_n = []
lig_n = []
ac = []
nac = []
sk = []
hh = []
rmsd = []
rfeats = []
rfeats1 = []
#dictionaries
rec_dic = {}
lig_dic = {}
cc_dic = {}
class docka:
    def native(self,x,y): #generate native file, rosetta input file 'native'
        x = self.x
        y = self.y
        z = self.z
        recl = []
        pepl = []
        with open(x,'r') as R:
            for line in R.readlines():
                if line.startswith('ATOM'):
                    recl.append(line.strip())
        with open(y,'r') as P:
            for line in P.readlines():
                if line.startswith('ATOM' or 'HETATM'):
                    pepl.append(line.strip())
        with open(z,'w') as N:
            for line in recl:
                N.write(line + '\n')
            N.write('TER\n')
            for line in pepl:
                N.write(line + '\n')
            N.write('END\n')
            

    def mutate(self,x,y): #peptide mutation part, mmtsb
        x = self.x
        y = self.y
        mtd = []
        os.environ['MMTSBDIR'] = '/lwork01/mmtsb/'
        mmtsb = os.environ['MMTSBDIR']
        os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + mmtsb + '/perl:' + mmtsb + '/bin'
        
        os.system('mutate.pl -seq 1:' + y + ' ' + x + ' > ' + y + '_mtd.pdb')
    #enva working functions
    '''    
    def aopt: #new accessibility option
        
    def bopt: #buffer option
        
    def mopt: #new accessibility option
        
    def eopt: #original accessibility option
    
    #text writing after enva working
    def acc:
        
    def new_acc:
        
    def rmsd:
        
    def skew:
        
    def Buff:
    '''  
    
    
class atoma:
    def reT(self,x): #reduce trimming part
        Olines = []
        recl = []
        pepl = []
        with open(x,'r') as F:
            for line in F.readlines():
                fline = line.strip()
                if fline.startswith('ATOM' or 'HETATM'):
                    if fline[-1] != 'H' :
                        Olines.append(fline)
        with open(x.split('.pdb')[0] + '_red.pdb','w') as rd:
            for line in Olines:
                if line.startswith('ATOM' or 'HETATM'):
                    if line[21] == 'A':
                        rd.write(line + '\n')
                    elif line[21] != 'A' :
                        rd.write(line[:21] + 'B ' + line[23:])
                elif line.strip() == '':
                    rd.write('TER\n')
        
    def orev(self,x,y): #reading original file for output rearrange
        x = self.x
        y = self.y
        rec_ser = []
        lig_ser = []
        rnn = 1
        lnn = 1
        with open(x,'r') as R:
            for line in R.readlines():
                rec_n.append(line[7:13].strip())
                if line[13:17].strip() == 'OXT' or line[13:17].strip() == '':
                    pass
                else:
                    if rnn == int(line[22:26].strip()):
                        rec_ser.append(line[13:17].strip())
                        rec_dic[rnn] = rec_ser
                    else:
                        rec_dic[rnn] = rec_ser
                        rnn = int(line[22:26].strip())
                        rec_ser = []
                        rec_ser.append(line[13:17].strip())
        with open(y,'r') as P:
            for line in P.readlines():
                lig_n.append(line[7:13].strip())
                if line[13:17].strip() == 'OXT' or line[13:17].strip() == '':
                    pass
                else:
                    if lnn == int(line[22:26].strip()):
                        lig_ser.append(line[13:17].strip())
                        lig_dic[lnn] = lig_ser
                    else:
                        lig_dic[lnn] = lig_ser
                        lnn = int(line[22:26].strip())
                        lig_ser = []
                        lig_ser.append(line[13:17].strip())
                        
        return rec_dic,lig_dic,rec_ser,lig_ser

    #def atrev: #output rearrange part