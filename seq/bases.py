import os
import sys
import multiprocessing as mp
import pandas as pd
import re

from loguru import logger
from adducts import Adducts

element = {'H':1.007825035,
           'D':2.014101779,
           '[2H]':2.014101779,
           'B':11.0093055,
           'C':12.0000000,
           '[13C]':13.00335483,
           'N':14.003074,
           '[15N]':15.00010897,
           'O':15.99491463,
           '[17O]':16.999131,
           '[18O]':17.9991603,
           'F':18.99840322,
           'P':30.973762,
           'S':31.9720707,
           'Cl':34.96885272,
           'As':74.9215942,
           'Br':78.9183361,
           'Se':79.9165196,
           'I':126.904473}

class Base:

    def __init__(self, row):
        self._name = str(row[0]).strip()
        self._longname = str(row[1]).strip()

        self._start = bool(int(row[2]))
        self._internal = bool(int(row[3]))
        self._end = bool(int(row[4]))
        self._formula = self.str2formula(str(row[5]))
        self._formulastr = str(row[5])
        self._mass = self.formula2mass(self._formula)
        #logger.info("mass {}", self._mass)

    def __str__(self):
        msg = "mass {} start {} internal {} end {}".format(self._mass, self._start, self._internal, self._end)
        return msg

    @property
    def mass(self):
        return self._mass

    @property
    def name(self):
        return self._name

    @property
    def start(self):
        return self._start

    @property
    def internal(self):
        return self._internal

    @property
    def end(self):
        return self._end

    @property
    def formulastr(self):
        return self._formulastr

    @property
    def longname(self):
        return self._longname

    def str2formula(self, s):
        formuladict = {p[0]:p[1] for p in re.findall("([\+\-]\[\w+\]|[\+\-]\w|\[\w+\]|\w)([0-9]*)",s)}
        for k in formuladict:
            if len(formuladict[k]) < 1:
                formuladict[k] = 1
            else:
                formuladict[k] = int(formuladict[k])

        keys = [k for k in formuladict.keys() if "+" in k or "-" in k]
        for key in keys:
            atom = re.findall("(\w)",key)[0]
            if "+" in key:
                if atom in formuladict:
                    formuladict[atom] += formuladict[key]
                else:
                    formuladict[atom] = formuladict[key]
            if "-" in key:
                if atom in formuladict:
                    formuladict[atom] -= formuladict[key]
                else:
                    formuladict[atom] = -1*formuladict[key]
            formuladict.pop(key, None)

        return formuladict

    def formula2mass(self, formuladict):
        return sum(element[e]*formuladict[e] for e in formuladict)

class DirectBase(Base):

    def __init__(self, row):
        self._name = str(row[0]).strip()
        self._mass = row[1]
        if row.shape[0] == 3:
            self._blast_name = str(row[2]).strip()
        else:
            self._blast_name = ''

    def __str__(self):
        msg = "name {} mass {}".format(self._name, self._mass)
        return msg

    @property
    def blast_name(self):
        return self._blast_name

class Bases:

    def __init__(self, csv_file=None, start=None):
        #self.bases_list = self.load_csv(csv_file)
        self._start = start
        self._max = None
        self._min = None
        self._adds = Adducts()
        self.bases_list = self.load_csv_directly(csv_file)
        self.terminal_list = self.load_terminal(csv_file)

    def load_csv_directly(self, csv_file):
        if not csv_file:
            cur_dir = os.path.dirname(os.path.abspath(__file__))
            csv_file = os.path.join(cur_dir, "statics/bases_modif.csv")
        logger.info("type {} name {}", type(csv_file), csv_file)
        #df = pd.read_csv(csv_file, nrows=6)
        df = pd.read_csv(csv_file)
        df = df[['Name', 'Exact Mass', 'BlastName']].dropna()
        df = df.sort_values(by='Exact Mass')
        logger.info("dataframe shape {}", df.shape)
        logger.info("dataframe {}", df)
        bases_list = list()
        for idx, row in df.iterrows():
            #logger.info("idx {} row {}", idx, row)
            base = DirectBase(row)
            #print(base)
            bases_list.append(base)

        return bases_list

    def load_terminal(self, csv_file):
        if not csv_file:
            cur_dir = os.path.dirname(os.path.abspath(__file__))
            csv_file = os.path.join(cur_dir, "statics/bases_terminal.csv")
        logger.info("type {} name {}", type(csv_file), csv_file)
        df = pd.read_csv(csv_file, nrows=4)
        df = df[['Name', 'Exact Mass']].dropna()
        df = df.sort_values(by='Exact Mass')
        logger.info("dataframe shape {}", df.shape)
        logger.info("dataframe {}", df)
        bases_list = list()
        for idx, row in df.iterrows():
            #logger.info("idx {} row {}", idx, row)
            base = DirectBase(row)
            #print(base)
            bases_list.append(base)

        return bases_list

    def load_csv(self, csv_file=None):
        if not csv_file:
            cur_dir = os.path.dirname(os.path.abspath(__file__))
            csv_file = os.path.join(cur_dir, "statics/bases.csv")

        df = pd.read_csv(csv_file, skiprows=6)
        logger.info("dataframe shape {}", df.shape)

        bases_list = list()
        for idx, row in df.iterrows():
            #logger.info("idx {} row {}", idx, row)
            #if idx > 5:
                #break
            if row[0] not in ['A', 'C', 'G', 'U', 'M']:
                continue
            base = Base(row)
            bases_list.append(base)
            #print("{},{},{},{},{},{},{}".format(base.name, base.longname, base.start, base.internal, base.end, base.formulastr, base.mass))
            #print(base)

        return bases_list

    @property
    def adds(self):
        return self._adds

    @property
    def max(self):
        if not self._max:
            m = 0
            for base in self.bases_list:
                if base.mass > m:
                    m = base.mass
            self._max = m
            logger.info("base MAX {}", self._max)

        return (self._max + self.adds.max) * 1.1

    @property
    def min(self):
        if not self._min:
            m = sys.maxsize
            for base in self.bases_list:
                if base.mass < m:
                    m = base.mass
            self._min = m

        return (self._min - self.adds.max) * 0.9

    @property
    def start(self):
        if self._start:
            return self._start
        #return 694.2397
        return 0

    def ppm2dm(self, m, ppm=10.0):
        return (m * ppm) / 1.0E6

    def dm2ppm(self, m, stringency):
        return abs((stringency * 1.0E6) / m)

    def match_start(self, df, start_mass):
        base_mass = start_mass
        stringency = self.ppm2dm(base_mass)
        df = df.iloc[(df.Mass-base_mass).abs().argsort()[:1]]
        df_match = df[abs(df.Mass - base_mass) <= stringency]
        if df_match.empty:
            logger.info("Created a start point.")
            df.loc[:, 'Mass'] = start_mass

        return df.iloc[0,:]

    def match_start_edge(self, diff, mass):
        #logger.info("stringency {} diff {}", stringency, diff)
        stringency = self.ppm2dm(mass)
        adds_combs = self.adds.combs()
        base_mass = self.start
        #for adduct in adds_combs:

        adduct = 0
        comb_mass = base_mass + adduct
        if (diff >= (comb_mass - stringency)) and (diff <= (comb_mass + stringency)):
            return True
        comb_mass = base_mass - adduct
        if (diff >= (comb_mass - stringency)) and (diff <= (comb_mass + stringency)):
            return True

        return False

    def match_base(self, diff, mass):
        #logger.info("stringency diff {} mass {}", diff, mass)
        adds_combs = self.adds.combs()

        for base in self.bases_list:
            #for adduct in adds_combs:
            for adduct in [0.0]:
                comb_mass = base.mass + adduct
                mass_theoretical = mass + comb_mass
                stringency = self.ppm2dm(mass_theoretical)
                if abs(diff - comb_mass) <= stringency:
                #if (diff >= (comb_mass - stringency)) and (diff <= (comb_mass + stringency)):
                    ppm = self.dm2ppm(mass_theoretical, diff-comb_mass)
                    #logger.info("base {} ppm {} adduct {} add_name {}", base, ppm, adduct, self.adds.adduct_name(adduct))
                    return base, ppm, self.adds.adduct_name(adduct)
                #comb_mass = base.mass - adduct
                #if (diff >= (comb_mass - stringency)) and (diff <= (comb_mass + stringency)):
                    #ppm = self.dm2ppm(mass, diff-comb_mass)
                    #return base, ppm
        return None, None, None

    def match_terminal_np(self, mass, df):
        idxs = list()
        base_names = dict()
        blast_names = dict()
        base_ppms = dict()
        for base in self.terminal_list:
            mass_theoretical = mass + base.mass
            stringency = self.ppm2dm(mass_theoretical)
            df_base = df[(abs(df.Mass - mass - base.mass) <= stringency) & (df.Vol > 6E6)].copy()
            #df_base = df[df.Vol > 1e7]
            func = self.dm2ppm
            #df_base['PPM'] = func(mass_theoretical, df_base['Mass'] - mass_theoretical)
            df_base.loc[:,'PPM'] = func(mass_theoretical, df_base['Mass'] - mass_theoretical)
            idxs.extend(list(df_base.index))
            for item in df_base.index:
                base_names[item] = base.name
                blast_names[item] = base.blast_name
                base_ppms[item] = round(df_base.loc[item, 'PPM'], 2)

        return list(set(idxs)), base_names, base_ppms, blast_names
    
    def match_bases_np(self, mass, df):
        adds_combs = self.adds.combs()
        idxs = list()
        base_names = dict()
        blast_names = dict()
        base_ppms = dict()
        for base in self.bases_list:
            for adduct in [0.0]:
            #for adduct in adds_combs:
                comb_mass = base.mass + adduct
                mass_theoretical = mass + comb_mass
                stringency = self.ppm2dm(mass_theoretical)
                df_base = df[abs(df.Mass - mass - comb_mass) <= stringency].copy()
                func = self.dm2ppm
                df_base.loc[:,'PPM'] = func(mass_theoretical, df_base['Mass'] - mass_theoretical)
                #idx = -1
                #if not df_base.empty:
                    #idx = df_base['Quality Score'].idxmax()
                    #idxs.append(idx)
                idxs.extend(list(df_base.index))
                for item in df_base.index:
                    base_names[item] = base.name
                    blast_names[item] = base.blast_name
                    base_ppms[item] = round(df_base.loc[item, 'PPM'], 2)

        return list(set(idxs)), base_names, base_ppms, blast_names

    def match_conjunct(self, t, h):
        conjs = self.bases_list[-2:]
        for conj in conjs:
            diff = t - (h + conj.mass)
            if diff <= 0:
                continue

            ret = self.match_base(diff, h + conj.mass)
            if ret:
                return True

        return False

if __name__ == '__main__':
    bs = Bases()
    logger.info("bases max {} min {}", bs.max, bs.min)
