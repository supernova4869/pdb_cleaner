import re

import pandas as pd
import numpy as np
import datetime


class PDB:
    def __init__(self, pdb_file=''):
        if pdb_file == '':
            self.title = []
            self.primary = []
            self.heterogen = []
            self.secondary = []
            self.connectivity_annotation = []
            self.miscellaneous = []
            self.crystallographic = []
            self.coordinate = []
            self.models = []
            self.atoms = []
            self.ters = []
            self.endmdls = []
            return
        with open(pdb_file) as f:
            lines = f.readlines()
            self.title = [l for l in lines if l.startswith('HEADER')
                          or l.startswith('OBSLTE')
                          or l.startswith('TITLE')
                          or l.startswith('CAVEAT')
                          or l.startswith('COMPND')
                          or l.startswith('SOURCE')
                          or l.startswith('KEYWDS')
                          or l.startswith('EXPDTA')
                          or l.startswith('AUTHOR')
                          or l.startswith('REVDAT')
                          or l.startswith('SPRSDE')
                          or l.startswith('JRNL')
                          or l.startswith('REMARK')]
            self.primary = [l for l in lines if l.startswith('DBREF')
                            or l.startswith('SEQADV')
                            or l.startswith('SEQRES')
                            or l.startswith('MODRES')]
            self.heterogen = [l for l in lines if l.startswith('HET')
                              or l.startswith('HETNAM')
                              or l.startswith('HETSYN')
                              or l.startswith('FORMUL')]
            self.secondary = [l for l in lines if l.startswith('HELIX')
                              or l.startswith('HELIX')]
            self.connectivity_annotation = [l for l in lines if l.startswith('SSBOND')
                                 or l.startswith('LINK')
                                 or l.startswith('CISPEP')]
            self.miscellaneous = [l for l in lines if l.startswith('SITE')]
            self.crystallographic = [l for l in lines if l.startswith('CRYST1')
                                     or l.startswith('MTRIXn')
                                     or l.startswith('ORIGXn')
                                     or l.startswith('SCALEn')]
            # 坐标部分
            self.coordinate = [l for l in lines if l.startswith('MODEL')
                               or l.startswith('ATOM')
                               or l.startswith('ANISOU')
                               or l.startswith('TER')
                               or l.startswith('HETATM')
                               or l.startswith('ENDMDL')]
            self.models = []
            self.atoms = []
            self.ters = []
            self.endmdls = []
            self.__parse_coordinates__()
            self.connectivity = [l for l in lines if l.startswith('CONECT')]
            self.bookkeeping = [l for l in lines if l.startswith('MASTER')
                                or l.startswith('END')]

    def __parse_coordinates__(self):
        for c in self.coordinate:
            if c.startswith('MODEL'):
                self.models.append(MODEL(c))
            elif c.startswith('ATOM') or c.startswith('HETATM'):
                self.atoms.append(ATOM(c))
            elif c.startswith('TER'):
                self.ters.append(TER(c))
            elif c.startswith('ENDMDL'):
                self.endmdls.append(ENDMDL(c))
            else:
                print('混进了一个奇怪的东西')

    def to_pdb(self, pdb_file):
        pdb_file = open(pdb_file, 'w')
        pdb_file.write(
            'REMARK   Created by supernova/PaperTools (https://gitee.com/supernova_bingbing/paper-tools)\n'
            'REMARK   Created: ' + str(datetime.datetime.now()) + '\n')
        for line in self.coordinate:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                pdb_file.write(str(ATOM(line)))
            elif line.startswith('TER'):
                pdb_file.write(str(TER(line)))
        pdb_file.write('END\n')
        pdb_file.close()


class ATOM:
    def __init__(self, line):
        self.typ = line[:6].strip()
        self.atnum = int(line[6:11].strip())
        self.atid = self.__atid_for_human__(line[12:16].strip())
        self.altloc = line[16].strip() or ' '
        self.resname = line[17:20].strip()
        self.chainid = line[21].strip() or 'X'
        self.resid = int(line[22:26].strip())
        self.icode = line[26].strip() or ' '
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.occupy = float(line[54:60].strip())
        self.bf = float(line[60:66].strip())
        self.atname = line[76:78].strip()
        self.charge = line[78:80].strip()

    def __atid_for_pdb__(self, atid):
        # 从方便理解的atid到pdb使用的atid
        if len(atid) == 4:
            if re.compile(r'H[ABGDEZH][1-9][1-9]').match(atid):
                return atid[-1] + atid[:-1]
            else:
                print('Warning: atom name seems incorrect:', atid)
        else:
            return ' ' + atid

    def __atid_for_human__(self, atid):
        # 从pdb使用的atid到方便理解的atid
        if len(atid) == 4:
            # 理论上应该是HG11, HG12这样
            if re.compile(r'[1-9]H[ABGDEZH][1-9]').match(atid):
                return atid[1:] + atid[0]
            else:
                print('Warning: atom name seems incorrect:', atid)
        else:
            return atid

    def __repr__(self):
        return '{0:6s}{1:5d} {2:4s}{3:1s}{4:3s} {5:1s}{6:4d}{7:1s}   {8:8.3f}{9:8.3f}{10:8.3f}{11:6.2f}{12:6.2f}' \
               '          {13:>2s}{14:>2s}\n' \
            .format(self.typ, self.atnum, self.__atid_for_pdb__(self.atid), self.altloc, self.resname, self.chainid,
                    self.resid, self.icode, self.x, self.y, self.z, self.occupy, self.bf, self.atname, self.charge)


class TER:
    def __init__(self, line):
        self.typ = line[:6].strip()
        self.atnum = int(line[6:11].strip())
        self.resname = line[17:20].strip()
        self.chainid = line[21].strip()
        self.resid = int(line[22:26].strip())
        self.icode = line[26].strip()
    
    def __repr__(self):
        return '{0:6s}{1:5d}      {2:3s} {3:1s}{4:4d}{5:1s}\n'.\
            format(self.typ, self.atnum, self.resname, self.chainid, self.resid, self.icode)


class MODEL:
    def __init__(self, line):
        self.typ = line[:6].strip()
        self.serial = int(line[10:14].strip())

    def __repr__(self):
        return '{0:6s}    {1:4d}\n'.\
            format(self.typ, self.serial)


class ENDMDL:
    def __init__(self, line):
        self.typ = line[:6].strip()

    def __repr__(self):
        return '{0:6s}\n'.\
            format(self.typ)


def csv2pdb(csv_file, pdb_file):
    pdb = PDB()
    a = np.array(pd.read_csv(csv_file, header=None))
    for line in a:
        if line[0].startswith('HETATM') or line[0].startswith('ATOM'):
            atomline = '{0:6s}{1:5d} {2:4s} {3:3s} {4:1s}{5:4d}    {6:8.3f}{7:8.3f}{8:8.3f}' \
                       '{9:6.2f}{10:6.2f}\n' \
            .format(line[0], line[1], line[2], line[3], line[4],
                    line[5], line[6], line[7], line[8],
                    float(line[9]), float(line[10]))
            pdb.atoms.append(ATOM(atomline))
            pdb.coordinate.append(atomline)
        else:
            terline = '{0:6s}{1:5d}      {2:3s} {3:1s}{4:4d}{5:1s}\n' \
            .format(line[0], int(line[1]), line[2], line[3], int(line[4]), ' ')
            ter = TER(terline)
            pdb.ters.append(TER(terline))
            pdb.coordinate.append(terline)
    pdb.__parse_coordinates__()
    pdb.to_pdb(pdb_file)


if __name__ == '__main__':
    pdb_zen = PDB('test/pdb/MOF.pdb')
    pdb_zen.to_pdb('test/pdb/MOF_short.pdb')
