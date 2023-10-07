import os
import numpy as np
import pandas as pd
from itertools import starmap, combinations

from loguru import logger

class Adducts:

    def __init__(self, csv_file=None):
        self._mass = np.array([], dtype=float)
        self._name = np.array([])
        self._action = np.array([], dtype=int)
        self.adducts_list = self.load_csv(csv_file)

    def load_csv(self, csv_file=None):
        if not csv_file:
            cur_dir = os.path.dirname(os.path.abspath(__file__))
            csv_file = os.path.join(cur_dir, "statics/adducts.csv")

        df = pd.read_csv(csv_file, skiprows=3, nrows=3)
        logger.info("dataframe shape {}", df.shape)

        adducts_list = list()
        for idx, row in df.iterrows():
            self._name = np.append(self._name, str(row[0]).strip())
            self._mass = np.append(self._mass, float(row[1]))
            self._action = np.append(self._action, int(row[2]))

        return adducts_list

    @property
    def mass(self):
        return self._mass

    @property
    def max(self):
        return max(self.combs())

    @property
    def min(self):
        return np.min(self._mass)
    
    @property
    def name(self):
        return self._name

    def find(self, m, stringency):
        return (self.name[np.nonzero((self.mass >= m - stringency) & (self.mass <= m + stringency))],
                self.mass[np.nonzero((self.mass >= m - stringency) & (self.mass <= m + stringency))])

    def combs(self):
        cbs = list()
        cbs.append(0.0)
        masses = list(self.mass[:2])
        cbs.extend(masses)

        _l2 = combinations(masses, 2)
        l2 = starmap(lambda x,y:x+y, _l2)
        #cbs.extend(l2)

        _l3 = combinations(masses, 3)
        l3 = starmap(lambda x,y,z:x+y+z, _l3)
        #cbs.extend(l3)

        cbs = list(set(cbs))
        return cbs

    def adduct_name(self, mass):
        if mass == 0.0:
            return None
        for idx, m in enumerate(self.mass):
            if m == mass:
                return self.name[idx]

    def adduct_value(self, name):
        if not name:
            return 0
        for idx, n in enumerate(self.name):
            if n == name:
                return self.mass[idx]

if __name__ == '__main__':
    bs = Adducts()
    cb = bs.combs()
    logger.info("len {} {}", len(cb), cb)
    logger.info("max {} min {}", bs.max, bs.min)
