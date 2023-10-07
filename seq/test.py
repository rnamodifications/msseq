import os
import sys
import itertools
import tempfile
import pandas as pd
from process import action
import multiprocessing as mp

DATASETS = [
        ('19.csv', 'AAAACCAGUCAGUCUACGC'),
        ('20.csv', 'CGCAUCUGACUGACCCGAUA'),
        ('21.csv', 'AGGGUUGACUCGAUUUAGGCG'),
        ('201.csv', 'GAGUCAUUACCAUUGCCAAA'),
        ('202.csv', 'CUCUCACAUCCUACAAAUGU'),
        ('203.csv', 'AGAACUCACAUUGAACUUAU'),
        ('204.csv', 'UAUUUCCCCUUCUACAUGCG'),
        ('206.csv', 'AUAAGUCAGGUACAGUCACA')#,
        #('205.csv', 'GCCAGGCCCUAGUGUACCGC'),
        #('30.csv', 'CAUAGCAUGCACCCCUAGGAGAUCUCAGCU'),
        #('207.csv', 'GGGUUGACUCGAUUUAGGCG')
        ]
COMBINED_SAMPLES_NUM = 2

class SeqDataset:
    def __init__(self, ds):
        self.filename = ds[0]
        self.seq = ds[1]
        self._df = self.get_dataframe(self.filename)

    @property
    def df(self):
        return self._df

    def get_dataframe(self, name):
        fpath = os.path.join('data/samples/', self.filename)
        read_data = pd.read_excel
        if name.split('.')[-1] == 'csv':
            read_data = pd.read_csv
        df = read_data(fpath)
        #df = pd.read_csv(fpath)
        df = df[['Mass', 'RT', 'Vol', 'Width', 'Quality Score']]
        return df

class MergedDataset:
    def __init__(self, df, seqs):
        self._df = df
        self._seqs = seqs
        self._fpath = None

    @property
    def df(self):
        return self._df

    @property
    def seqs(self):
        return self._seqs

    @property
    def fpath(self):
        if not self._fpath:
            self._fpath = self.write_tmp(self._df)
        return self._fpath

    def write_tmp(self, df):
        _, tmp_path = tempfile.mkstemp(suffix='.csv')
        df.to_csv(tmp_path)
        return tmp_path

def make_combination(N):
    if N < 1 or N > len(DATASETS):
        print("cannot handle {} length combinations".format(N))
        return

    combs = list(itertools.combinations(DATASETS, N))
    print(combs)
    return combs

def test_combined_datasets(combs):
    if not combs:
        print("cannot handle None combinations")
        return
    merged_datasets = [combine_datasets(comb) for comb in combs]
    return merged_datasets

def evaluate_merged_dataset(merged_dataset, a=True):
    def common(a, b):
        a_set = set(a)
        b_set = set(b)
        return a_set & b_set

    seqs = merged_dataset.seqs
    fpath = merged_dataset.fpath
    data = action(fpath, 3)
    exp_seqs = data['data'].get('sequences')
    fw_crosstalk = data['data']['fw_crosstalk']
    bw_crosstalk = data['data']['bw_crosstalk']
    fw_modification = data['data']['fw_modification']
    bw_modification = data['data']['bw_modification']
    if exp_seqs and seqs:
        print("seqs {} experimental_seqs {}".format(seqs, exp_seqs))
        exp_seqs.sort()
        seqs.sort()
        if exp_seqs == seqs:
            return (len(seqs), len(exp_seqs), fw_crosstalk, bw_crosstalk, fw_modification, bw_modification)

    print("Failed seqs {} experimental_seqs {}".format(seqs, exp_seqs))
    common_set = common(seqs, exp_seqs)
    return (len(seqs), len(common_set), fw_crosstalk, bw_crosstalk, fw_modification, bw_modification)

def process_merged_datasets(merged_datasets):
    with mp.Pool(mp.cpu_count()) as p:
            evals = p.starmap(evaluate_merged_dataset, [(ds, True) for ds in merged_datasets])
    #evals = [evaluate_merged_dataset(ds) for ds in merged_datasets]
    print(evals)
    theory = sum([item[0] for item in evals])
    practice = sum([item[1] for item in evals])
    result = practice * 1.0 / theory
    print("evaluating result {:.2f}%".format(result * 100))
    fw_crosstalks = [item[2] for item in evals]
    bw_crosstalks = [item[3] for item in evals]
    fw_ct_mean = sum(fw_crosstalks) / len(fw_crosstalks)
    bw_ct_mean = sum(bw_crosstalks) / len(bw_crosstalks)
    print("forward crosstalks mean {} {}\nbackward crosstalks mean {} {}".format(fw_ct_mean, fw_crosstalks, bw_ct_mean, bw_crosstalks))
    fw_modification = [item[4] for item in evals]
    bw_modification = [item[5] for item in evals]
    fw_md_mean = sum(fw_modification) / len(fw_modification)
    bw_md_mean = sum(bw_modification) / len(bw_modification)
    print("forward modification mean {} {}\nbackward modification mean {} {}".format(fw_md_mean, fw_modification, bw_md_mean, bw_modification))
def combine_datasets(comb):
    ds_num = len(comb)
    ds = list()
    seq_datasets = [SeqDataset(comb[i]) for i in range(ds_num)]
    dfs = [seq_ds.df for seq_ds in seq_datasets]
    seqs = [seq_ds.seq for seq_ds in seq_datasets]

    merged_dfs = pd.concat(dfs)
    merged_dataset = MergedDataset(merged_dfs, seqs)
    return merged_dataset

def main(num=0):
    if num == 0:
        num = COMBINED_SAMPLES_NUM
    combs = make_combination(num)
    merged_datasets = test_combined_datasets(combs)
    process_merged_datasets(merged_datasets)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        num = int(sys.argv[1])
    else:
        num = 0
    main(num)
