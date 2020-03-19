import pandas as pd
import os,sys,glob,psutil,shutil
import multiprocessing
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
def reading_DB(x):
    hladb = pd.read_csv(x,sep='\t')
    indb = hladb[hladb['PDB'] == pdbid.upper()].dropna()
    print indb
    o_seq = indb['Peptide Sequence'].values[0]
    #hlaclass = indb['HLA CLASS'].values[0]
    atype = indb['Alpha TYPE'].values[0]
    btype = indb['Beta TYPE'].values[0]
    aclass = indb['Alpha CLASS'].values[0]
    bclass = indb['Beta CLASS'].values[0]
    ffreq = indb['Freq'].values[0]
    seqlen = len(o_seq)
    hla = aclass + ''.join(atype.split(':')) + '_' + bclass + ''.join(btype.split(':'))
    return ffreq,hla,seqlen,o_seq
def reT(x):
    nam = x.split('.')
    num = nam[3]
    nam1 = '.'.join([nam[0], nam[3], nam[2]]))
    Olines = []
    recl = []
    pepl = []
    with open(x,'r') as F:
        for line in F.readlines():
            tline = line.strip()
            if tline.startswith('ATOM') :
                if tline[-1] != 'H':
                    Olines.append(tline)
    for line in Olines:
        if line[21] == 'A':
            recl.append(line)
        elif line[21] != 'A':
            pepl.append(line[:21] + 'B ' + line[23:])
    with open(nam1 + '.red.pdb','w') as W:
        for line in recl:
            W.write(line + '\n')
        W.write('TER\n')
        for line in pepl:
            W.write(line + '\n')

def revise_origin(x):
    rec_ser = []
    lig_ser = []
    rnn = 1
    nnn = 1
    with open(RECLIB + x + '_rec.pdb','r') as R: #receptor reading
        for line in R.readlines():
            rec_n.append(line[7:13].strip())
            if line[13:17].strip() == 'OXT' or line[13:17].strip() == '' :
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
    with open('../' + x,'r') as P: #peptide reading
        for line in P.readlines():
            lig_n.append(line[7:13].strip())
            if line[13:17].strip() == 'OXT' or line[13:17].strip() == '':
                pass
            else:
                if nnn == int(line[22:26].strip()):
                    lig_ser.append(line[13:17].strip())
                    lig_dic[nnn] = lig_ser
                else:
                    lig_dic[nnn] = lig_ser
                    nnn = int(line[22:26].strip())
                    lig_ser = []
                    lig_ser.append(line[13:17].strip())

def atom_revise(x):
    nam = x.split('.')
    revnam = '.'.join([nam[0],nam[1],'rev',nam[3]]
    rnn = 1
    nnn = 1
    olig_ser = {}
    orec_ser = {}
    lig_dic2 = {}
    rec_dic2 = {}
    lig_lines = []
    out_lines = []

    with open(x,'r') as F:
        for line in F.readlines():
            if line[21:23].strip() == 'A':
                if rnn == int(line[23:26].strip()):
                    orec_ser[line[13:17].strip()] = line
                    rec_dic2[rnn] = orec_ser
                else:
                    rec_dic2[rnn] = orec_ser
                    rnn = int(line[23:26].strip())
                    orec_ser = {}
                    orec_ser[line[13:17].strip()] = line
            elif line[21:23].strip() == 'B':
                if nnn == int(line[23:26].strip()):
                    olig_ser[line[13:17].strip()] = line
                    lig_dic2[nnn] = olig_ser
                else:
                    lig_dic2[nnn] = olig_ser
                    nnn = int(line[23:26].strip())
                    olig_ser = {}
                    olig_ser[line[13:17].strip()] = line
    with open('feat.list','r') as feat:
        for i in feat.readlines():
            for j in lig_dic[i]:
                lig_lines.append(lig_dic2[i][j].strip())
        for i in rec_dic:
            for j in rec_dic[i]:
                out_lines.append(rec_dic[i][j].strip())
    with open(revnam,'w') as W:
        for i,j in zip(out_lines,rec_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:21] + 'A ' + i[23:] + '\n')
        W.write('TER\n')
        for i,j in zip(lig_lines,lig_n):
            W.write(i[:7] + j.rjust(4) + '  ' + i[13:21] + 'B ' + i[23:] + '\n')
        W.write('END\n')

def txt_writing(x):
    ff_list = []
    hh_list = []
    tt_acc = []
    tt_phi = []
    tt_psi = []
    mat_lst = []
    ave_acc = []
    tot_lig_acc = 0
    ave_lig_acc = 0
    tot_lig_num = 0
    pdbl = []
    tttt = []
    newl = []

    for i in sorted(glob.glob('*rev.pdb')):
        os.system('csplit -f \'' + i + '_\' ' + i + ' \'/TER/\' --quiet')
        os.system('rmsd_total ' + PEPLIB + '/' + pdbid + '_pep.pdb ' + i + '_01 bb >> ' + x + '_rmsd.txt')
    for f, g in zip(sorted(glob.glob('*a.out')), sorted(glob.glob('*m.out'))):
        sidx = []
        pdb = '_'.join([pdbid, str(x),f.split('.')[0].split('_')[6]])
        pdbl.append(pdb)
        adf = pd.read_csv(f, sep='\s+', skiprows=[0, 0], header=None)
        mdf = pd.read_csv(g, sep='\s+', header=None)
        mdf_filter = mdf[mdf[17] == 1].iloc[:, 9:16]
        for i in mdf_filter[9]:
            tot_lig_acc = tot_lig_acc + float(i)
            tot_lig_num += 1
            if tot_lig_num > 0:
                ave_lig_acc = tot_lig_acc / float(tot_lig_num)
            else:
                ave_lig_acc = 0
        ave_acc.append(ave_lig_acc)

        new_acc = {}
        num = 0
        for i in rfeats:
            new_acc['AA_' + cc_dic[str(i)]] = sum([adf[adf[1] == int(i.split('_')[0])][adf[2] == i.split('_')[2]][adf[3] == i.split('_')[1]][11].values,1])
        for i in mdf_filter.T:
            tidx = str(mdf_filter.T[i][13]) + '_' + str(mdf_filter.T[i][14]) + '_' + str(mdf_filter.T[i][15])
            if tidx in rfeats:
                tttt.append(tidx)
                zdx = rfeats.index(tidx) + 1
                sidx.append(zdx)
                num += 1
        newl.append(pd.DataFrame(new_acc))
        if not os.path.exists(pdb.split('.')[0] + '.ser'):
            with open(pdb.split('.')[0] + '.ser', 'w') as W:
                for sid in list(set(sidx)):
                    W.write(str(sid) + '\n')
        iskew_cmd = 'iskew ' + pdb.split('.')[0] + '.ser >> total_' + x + '_sk.txt'
        subprocess.call(iskew_cmd, shell=True)
        ratio = float(num) / float(len(rfeats))
        mat_lst.append(ratio)
    new_ac = pd.concat(newl)
    new_ac['total_rec_acc'] = new_ac.sum(axis=1)
    new_ac['ave_lig_acc'] = ave_acc
    new_ac['%Match'] = mat_lst
    new_ac['PDB'] = pdbl
    new_ac.set_index('PDB').reset_index()

    new_ac.to_csv(x + '_nac.txt', sep='\t', index=False)
    nn_list = []
    for f in sorted(glob.glob('*b.out')):
        with open(f, 'r') as F:
            n = 0
            for line in F.readlines():
                tt = line[7:13].strip() + '_' + line[17:21].strip() + '_' + line[13:17].strip()
                if tt in rfeats:
                    n += 1
            nn_list.append(n)
        with open(f, 'r') as F:
            hh_list.append(len(F.readlines()))
        ff_list.append(pdbid + '_' + str(x) + '_' + f.split('.')[0].split('_')[6])
    with open(x + '_hh.txt', 'w') as W:
        W.write('PDB\tN.of.BB_full\tN.of.BB_feat\n')
        for i, j, k in zip(ff_list, hh_list, nn_list):
            W.write(str(i) + '\t' + str(j) + '\t' + str(k) + '\n')
    acc_columns = ['P%d' % i for i in range(1, seqlen + 1)]
    phi_columns = ['PHI%d' % i for i in range(1, seqlen + 1)]
    psi_columns = ['PSI%d' % i for i in range(1, seqlen + 1)]
    pdbl = []
    for f in sorted(glob.glob('*.env')):
        pdb = '_'.join([pdbid, str(x), f.split('.')[0].split('_')[6]])
        pdbl.append(pdb)
        with open(f, 'r') as F:
            acc = []
            psi = []
            phi = []
            for line in F.readlines():
                if line.find('chain') and line.startswith('ATOM'):
                    envs = line.split()[9:]
                    acc.append(envs[2])
                    phi.append(envs[4])
                    psi.append(envs[5])
            tt_acc.append(pd.DataFrame(acc, index=acc_columns).T)
            tt_phi.append(pd.DataFrame(phi, index=phi_columns).T)
            tt_psi.append(pd.DataFrame(psi, index=psi_columns).T)
    pdbdf = pd.DataFrame(pdbl, columns=['PDB'])
    acc_df = pd.concat(tt_acc, ignore_index=True).reindex()
    phi_df = pd.concat(tt_phi, ignore_index=True).reindex()
    psi_df = pd.concat(tt_psi, ignore_index=True).reindex()
    act_df = pd.concat([pdbdf, acc_df, phi_df, psi_df], axis=1)
    act_df.to_csv(x + '_ac_ct.txt', sep='\t', index=False)

def sheba_enva(x):
    nam = x.split('.')
    num = nam[1]
    nana = pdbid + '_' + num
    namtrf = '.'.join([nam[0],nam[1],nam[2],'trf'])
    outform1 = pdbid + '_' + num + '.pdb'
    outform2 = pdbid + '_' + num + '_new.pdb'
    outform3 = pdbid + '_' + num + '_het.pdb'
    os.system('sheba_01 -x ' + refer + ' ' + nam)
    os.system('sheba_01 -t ' + namtrf + ' ' + nam)
    shutil.move(nam + '.pdb',outform1)
    os.system('python clean_pdb.py ' + outform1 + ' B')
    os.system('csplit -f \'%s_\' %s.pdb \'/TER/\' --quiet')
    os.system('cat %s_rec.pdb %s_B.pdb > %s_new.pdb')
    os.system('sed \'s/ATOM  /HETATM/\' %s_01 > %s_02')
    os.system('cat %s_rec.pdb %s_02 > %s_het.pdb')

    if not os.path.exists(nana + '_a.out'):
        os.system('enva.v2 -a ' + outform2 + '> ' + nana + '_a.out')
    if not os.path.exists(nana + '_b.out'):
        os.system('enva.v2 -b ' + outform3 + '> ' + nana + '_b.out')
    if not os.path.exists(nana + '_m.out'):
        os.system('enva.v2 -m ' + outform2 + ' ' + pdbid +  'B > ' + nana + '_m.out' )
    os.system('enva.v2 -e ' + outform2 + ' B')

os.environ['neogear'] = '/awork06-1/neoscan_gear/'
gear = os.environ['neogear']
os.environ['PATH'] += ':' + os.environ['PATH'] + ':' + gear
PDBLIB = '/awork06-1/YKLee/hla_class2/' + ffreq + '/' +
RECLIB = '/awork06-1/YKLee/hla_class2/' + ffreq + '/' +
PEPLIB = '/awork06-1/YKLee/hla_class2/' + ffreq + '/' +
pdbid = sys.argv[1]
infile = sys.argv[2]
(ffreq,hla,seqlen,o_seq) = reading_DB(infile)
refer = PDBLIB + '/' + inpdb + '_native.pdb'
ncpu = str(psutil.cpu_count() - 1)
feats = pd.read_csv(PDBLIB + '/' + hla + '_native/' + pdbid + '.out',sep='\t',header=None)
feats_filt = feats[feats[17] == 1].iloc[:,13:16]
for i in feats_filt.values.tolist():
    rfeats.append(str(i[0]) + '_' + '_'.join(i[1:]))
feats1_filt = feats[feats[17] == 1].iloc[:,12:16]
for i in feats1_filt.values.tolist():
    rfeats1.append(str(i[0]) + '_' + '_'.join(i[2:]))
for i,j in zip(rfeats,rfeats1):
    cc_dic[i] = j
