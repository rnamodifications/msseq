import os
import sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from statsmodels.nonparametric.smoothers_lowess import lowess
import click

from loguru import logger
from bases import Bases

class Compound:

    def __init__(self, row):
        self.special = False
        if 'Special' in row:
            self.special = row['Special']
        self.rt = row['RT']
        self.mass = row['Mass']
        self.vol = row['Vol']
        if 'Width' in row:
            self.width = row['Width']
        else:
            self.width = 0.8
        if 'Quality Score' in row:
            self.score = row['Quality Score']
        else:
            self.score = 100

        self.hit = False
        self.terminal_hit = False
        self.refine_hit = False

    def __str__(self):
        msg = "RT {} Mass {} Vol {} Width {} Score".format(self.rt, self.mass, self.vol, self.width, self.score)
        return msg

class Compounds:

    def __init__(self, csv_file=None, label=None):
        self._cpd = self.load_csv(csv_file)
        self._bases = Bases(start=label)

    @property
    def cpd(self):
        return self._cpd

    @property
    def bases(self):
        return self._bases

    def load_csv(self, csv_file=None):
        if not csv_file:
            cur_dir = os.path.dirname(os.path.abspath(__file__))
            csv_file = os.path.join(cur_dir, "statics/compounds.csv")

        self.df = None
        func = pd.read_excel
        if csv_file[-3:] == 'csv':
            func = pd.read_csv

        try:
            cols = ['Mass', 'RT', 'Vol', 'Width', 'Quality Score']
            self.df = func(csv_file, usecols=cols)
        except ValueError as ve:
            logger.error("load csv error {}", ve)
            cols = ['Mass', 'RT', 'Vol']
            self.df = func(csv_file, usecols=cols)

        #self.df = self.df[self.df.Mass < 8000]
        self.df = self.df.sort_values(by='Mass')
        logger.info("dataframe shape {}", self.df.shape)

        compounds_list = list()
        for idx, row in self.df.iterrows():
            compound = Compound(row)
            compounds_list.append(compound)

        return compounds_list

    def _cpds_nearby(self, idx):
        cpd = self.cpd[idx]
        cpds = list()
        j = idx - 1
        while j >= 0:
            diff = cpd.mass - self.cpd[j].mass 
            if diff > self.bases.max:
                break

            '''
            if diff < self.bases.min:
                j -= 1
                continue
            '''

            diff_rt = abs(cpd.rt - self.cpd[j].rt)
            if diff_rt > 0 and diff_rt < cpd.width * 5.0:
                cpds.append(self.cpd[j])
            j -= 1

        return cpds

    def _find_cpd_match(self, cpd, nearby_cpds):
        cpds = list()
        max_vol = 0
        _cpd = None
        chose_ppm = 0
        base = None
        if not nearby_cpds:
            return _cpd, base

        for nearby_cpd in nearby_cpds:
            if not cpd or not nearby_cpd:
                continue

            diff = cpd.mass - nearby_cpd.mass
            base_match, ppm = self.bases.match_base(diff, cpd.mass)
            if base_match:
                if nearby_cpd.vol > max_vol:
                    max_vol = nearby_cpd.vol
                    base = base_match
                    _cpd = nearby_cpd
                    chose_ppm = ppm
                #cpds.append(nearby_cpd)

        #cpd = max(cpds, key=lambda item:item.vol)
        return _cpd, base

    def find_nearby_cpd(self, cpd, idx):
        nearby_cpds = self._cpds_nearby(idx)
        if not nearby_cpds:
            return None, None
        #cpd = self.cpds[idx]
        _cpd, _base = self._find_cpd_match(cpd, nearby_cpds)
        return _cpd, _base

    def generate_edges(self):
        edges = list()
        i = len(self.cpd) - 1
        while i > 0:
            cpd, base = self.find_nearby_cpd(self.cpd[i], i)
            if cpd:
                edge = (self.cpd[i], cpd, base)
                edges.append(edge)
            i -= 1

        return edges

    def made_edges(self, G, edges):
        logger.info("adding edges to Graph...")
        for edge in edges:
            G.add_edge(edge[0], edge[1], base=edge[2])

        return G
    
    def check_graph_edges(self, G):
        for edge in G.edges():
            s = edge[0]
            e = edge[1]
            if s.mass < e.mass:
                logger.warn("--->WRONG---{} {}", edge[0], edge[1])
                return False

        return True

    def clean_invalid_edges(self, G):
        logger.info("removing invalid edges...")
        invalid_edges = list()
        for level1node in G.nodes():
            for level2node in G.successors(level1node):
                edge = G[level1node][level2node]
                base = edge['base']
                if base.end and (not base.start) and (not base.internal):
                    invalid_edges.append((level1node, level2node))

                for level3node in G.successors(level2node):
                #if G.successors(level2node):
                    edge = G[level2node][level3node]
                    base = edge['base']
                    if (not base.internal) and (not base.end):
                        invalid_edges.append((level2node, level3node))
        logger.info("{} edges, {} invalid edges to remove", len(list(G.edges())), len(invalid_edges))
        invalid_edges = list(set(invalid_edges))
        for edge in invalid_edges:
            G.remove_edge(edge[0], edge[1])
        return G

    def lead_tail_nodes(self, G):
        non_leading_nodes = list()
        tail_nodes = list()
        leading_nodes = list()
        for node in G.nodes():
            pre_nodes = list(G.predecessors(node))
            sub_nodes = list(G.successors(node))
            non_leading_nodes.extend(sub_nodes)
            if (not pre_nodes) and sub_nodes:
                leading_nodes.append(node)
            #if (not pre_nodes) and (not sub_nodes):
                #non_leading_nodes.append(node)
            elif pre_nodes and (not sub_nodes):
                tail_nodes.append(node)
        '''
        non_leading_nodes = list(set(non_leading_nodes))
        leading_nodes_filter = filter(lambda node: node not in non_leading_nodes, G.nodes())
        leading_nodes = list(leading_nodes_filter)
        '''
        logger.info("{} leading nodes {} tail nodes", len(leading_nodes), len(tail_nodes))
        return leading_nodes, tail_nodes

    def longest_paths(self, G, leading_nodes, tail_nodes):
        logger.info("generating LONGEST edges...")
        full_paths = list()
        for lead in leading_nodes:
            for tail in tail_nodes:
                paths = nx.all_simple_paths(G, lead, tail)
                paths = list(paths)
                if paths:
                    path_size_maps = map(lambda x: len(x), paths)
                    max_size = max(list(path_size_maps))
                    max_paths = [path for path in paths if len(path) == max_size]
                    full_paths.extend(max_paths)

        logger.info("total paths {}", len(full_paths))
        paths_len = [len(path) for path in full_paths]
        max_size = max(paths_len)
        full_paths = [path for path in full_paths if len(path) > max_size - 1]
        logger.info("paths max_len {} full_paths {}", max_size, len(full_paths))
        return full_paths
                
    def longest_path(self, G, suppress=3):
        path = nx.dag_longest_path(G)
        path_len = nx.dag_longest_path_length(G)
        logger.info("currently longest path {}", path_len)
        if path_len <= 3:
            return None
        return path

    def conjunct_paths(self, paths):
        heads = [path[0] for path in paths]
        tails = [path[-1] for path in paths]
        matches = list()
        for idx_t, t in enumerate(tails):
            for idx_h, h in enumerate(heads):
                if idx_t == idx_h or t.mass <= h.mass:
                    continue
                if self.bases.match_conjunct(t.mass, h.mass):
                    matches.append((idx_t, idx_h))

        logger.info("matches {}", matches)
        new_paths = list()
        for m in matches:
            path = list()
            path.extend(paths[m[0]])
            path.extend(paths[m[1]])
            new_paths.append(path)

        return new_paths

    def top_n_longest_path(self, G, N):
        idx = 0
        paths = list()
        seqs = list()
        while idx < N:
            idx += 1
            if G.__len__() == 0:
                break

            logger.info("G size to be handle {}", G.number_of_nodes())
            path = self.longest_path(G)
            if not path:
                break

            paths.append(path)

            seq = self.sequences_of_path(G, path)
            if seq:
                seqs.append(seq)

            for node in path:
                if G.has_node(node):
                    G.remove_node(node)

        logger.info("paths size {}", len(paths))
        return paths, seqs

    def dfs_edges(self, G):
        logger.info("generating DFS edges...")
        edges = list()

        # calculate max length
        max_size = 0
        for node in G.nodes():
            dfs_edge = nx.dfs_edges(G, node)
            dfs_edge_list = list(dfs_edge)
            length = len(dfs_edge_list)
            if length > max_size:
                max_size = length

        # collect internal nodes
        non_leading_nodes = list()
        for node in G.nodes():
            pre_nodes = list(G.predecessors(node))
            sub_nodes = list(G.successors(node))
            non_leading_nodes.extend(sub_nodes)
            if (not pre_nodes) and (not sub_nodes):
                non_leading_nodes.append(node)

        non_leading_nodes = list(set(non_leading_nodes))
        logger.info("collected {} internal nodes from {} total nodes", len(non_leading_nodes), G.size())
        leading_nodes_filter = filter(lambda node: node not in non_leading_nodes, G.nodes())
        leading_nodes = list(leading_nodes_filter)

        edges_lens = list()
        '''
        for node in G.nodes():
            if node in non_leading_nodes:
                continue
        '''
        for node in leading_nodes:

            dfs_edge = nx.dfs_edges(G, node)
            dfs_edge_list = list(dfs_edge)
            length = len(dfs_edge_list)
            edges_lens.append(length)
            if length > max_size-3:
                logger.info("node {} len {}", node, length)
                edges.append(dfs_edge_list)
                if length > max_size-2:
                    bases = list()
                    for edge in dfs_edge_list:
                        base = G[edge[0]][edge[1]]['base']
                        if isinstance(base, str):
                            bases.append(base)
                        else:
                            bases.append(base.name)
                    logger.info("bases {}", bases)

        logger.info("collected {} edges {} lens", len(edges), edges_lens)
        return edges

    def generate_graph(self):
        logger.info("generating Graph...")
        G = nx.DiGraph()
        for compound in self.cpd:
            #G.add_node((compound.rt, compound.mass))
            G.add_node(compound)

        #nx.draw(G)
        #plt.show()
        return G

    def sequences_of_path(self, G, path):
        logger.info("sequences_of_path")
        DG = nx.DiGraph()
        DG.add_path(path)
        seq = list()
        for edge in DG.edges():
            base = G[edge[0]][edge[1]]['base']
            seq.append(base)
        logger.info("sequence {}", seq)
        return seq

    def print_sequences(self, seqs):
        for idx, seq in enumerate(seqs):
            print("seq {} {}".format(idx, [item.name for item in seq]))

    def draw_scatters(self, ax):
        rt = list()
        mass = list()
        vol = list()
        df = self.df.sort_values(by='Vol', ascending=False)
        df = df.head(300)
        count = 0
        for idx, row in df.iterrows():
            if row['RT'] > 56:
                continue
            rt.append(row['RT'])
            mass.append(row['Mass'])
            vol.append(row['Vol'])
            count += 1
        #for compound in self.cpd:
            #rt.append(compound.rt)
            #mass.append(compound.mass)
            #vol.append(compound.vol)
        mp = list(range(count))
        mp = [i*100 for i in mp]
        plt.scatter(mass, rt, s=100, c=mp, cmap=cm.cool)
        plt.tick_params(labelsize=14)
        font = {'family': 'sans-serif', 'size': 14}
        ax.set_xlabel('Mass (Da)',font )
        ax.set_ylabel('Retention Time (min)', font)
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 14
        cbar.ax.set_ylabel('Intensity (counts)', size=14, rotation=270)

        #plt.show()
        #plt.savefig('compounds.png')

    def draw_paths(self, paths, ax, G):
        for path in paths:
            x = list()
            y = list()
            z = list()
            for cpd in path:
                x.append(cpd.mass)
                y.append(cpd.rt)
                z.append(cpd.vol)

                #base = nx.get_node_attributes(G, 'base')
                node = self.node(G, cpd.mass)
                if not node:
                    continue
                try:
                    node = G.nodes[node]
                    base = node['base']
                    logger.info("node base {}", base)
                    plt.annotate(s=base, xy=(cpd.mass, cpd.rt), xytext=(-10, 8), textcoords='offset points', size=13, weight='heavy')
                except KeyError as ke:
                    logger.info("KeyError {}", ke)
            #logger.info("path {}", list(zip(x, y)))
            logger.info("path {}", x)
            plt.plot(x, y, color='#708090', linewidth=1, zorder=1)
            plt.scatter(x, y, s=100, linewidth=1, c='#CF00FF', zorder=2)
            #plt.imshow(x, y, cmap=plt.cm.hot, vmin=0, vmax=1)

        return x,y

    def draw_lines(self, paths):
        x = list()
        y = list()
        for path in paths:
            for edge in path:
                start_point = edge[0]
                end_point = edge[1]
                x.append(start_point.mass)
                y.append(start_point.rt)
            plt.plot(x, y)
            #logger.info("path {}", list(zip(x,y)))

    def midline(self):
        m = np.array([d.mass for d in self.cpd])
        t = np.array([d.rt for d in self.cpd])
        midline = lowess(t, m, frac=0.25, delta=0.05)
        print(midline)
        return midline

    def create_new_points(self, G, masses):
        nodes = [self.create_new_point(G, mass) for mass in masses]
        nodes = list(set(nodes))
        return nodes

    def create_new_point(self, G, mass):
        node = self.node(G, mass)
        if node:
            return node
        row = dict()
        row['Mass'] = mass
        row['RT'] = 0.0
        row['Vol'] = 0.0
        node = self.fill_start(row)
        return node

    def find_start_point(self, G, start_mass, fill=True):
        row = self.bases.match_start(self.df, start_mass)
        node = self.node(G, row['Mass'])
        if not node and fill:
            node = self.fill_start(row)

        return node
        

    def find_start_edge(self):
        edges = list()
        i = len(self.cpd) - 1
        while i > 0:
            j = i - 1
            base = None
            while j >= 0:
                diff = self.cpd[i].mass - self.cpd[j].mass 
                if diff > self.bases.max:
                    break
                if diff < self.bases.min:
                    j -= 1
                    continue

                #diff_rt = self.cpd[i].rt - self.cpd[j].rt 
                #if diff_rt > 0:
                base = self.bases.match_start_edge(diff, self.cpd[j].mass)
                if base:
                    #break
                    edge = (self.cpd[i], self.cpd[j], base)
                    edges.append(edge)

                j -= 1

            #if base:
                #edge = (self.cpd[i], self.cpd[j], base)
                #edges.append(edge)
            i -= 1

        #logger.info("start edges {} {}", len(edges), [(edge[0].mass, edge[1].mass) for edge in edges])
        logger.info("start edges {} {}", len(edges), [edge[0].mass - edge[1].mass for edge in edges])
        return edges

    def find_terminals(self, G, nodes):
        df = self.df.sort_values(by="Mass", ascending=False)
        idxs = list()
        ret_nodes = list()
        for node in nodes:
            if node.terminal_hit:
                continue
            _df = df[(df.Mass > node.mass + self.bases.min) & (df.Mass < node.mass + self.bases.max)]
            #_df = _df[abs(_df.RT - node.rt) < node.width * 20.0]
            _idxs, base_names, base_ppms, blast_names = self.bases.match_terminal_np(node.mass, _df)
            if _idxs:
                idxs.extend(_idxs)

            node.terminal_hit = True

            vols = 0
            if len(_idxs) == 0:
                continue

            for idx in _idxs:
                row = self.df.loc[idx, :]
                _node = self.node(G, row['Mass'])
                if not _node:
                    continue
                vols += _node.vol
            aver_vol = vols / len(_idxs)

            for idx in _idxs:
                row = self.df.loc[idx, :]
                _node = self.node(G, row['Mass'])
                ret_nodes.append(_node)
                base_name = base_names.get(idx)
                blast_name = blast_names.get(idx)
                #G.add_edge(node, _node, base=base_name, ppm=base_ppms.get(idx))
                G.add_edge(node, _node, base=base_name, blast_name=blast_name, ppm=base_ppms.get(idx))
                G.nodes[_node].update({"terminal": True})
        return ret_nodes

    def find_next_multi(self, G, nodes):
        df = self.df.sort_values(by="Mass", ascending=False)
        idxs = list()
        ret_nodes = list()
        for node in nodes:
            print("processing node {}".format(node.mass))
            if node.hit:
                continue
            _df = df[(df.Mass >= node.mass + self.bases.min) & (df.Mass <= node.mass + self.bases.max)]
            #_df = _df[abs(_df.RT - node.rt) < node.width * 20.0]
            _idxs, base_names, base_ppms, blast_names = self.bases.match_bases_np(node.mass, _df)
            if _idxs:
                idxs.extend(_idxs)

            node.hit = True

        #if not idxs:
            #return ret_nodes
            vols = 0
            if len(_idxs) == 0:
                continue

            #for idx in _idxs:
                #row = self.df.loc[idx, :]
                #_node = self.node(G, row['Mass'])
                #vols += _node.vol
            #aver_vol = vols / len(_idxs)

            for idx in _idxs:
                row = self.df.loc[idx, :]
                _node = self.node(G, row['Mass'])
                if not _node:
                    print("None node found ")
                    continue
                #if _node.vol < aver_vol * 1 / 2.0:
                    #logger.info("exclude mass {} score {} ppm {:.2f}", _node.mass, _node.score, base_ppms.get(idx))
                    #continue
                #logger.info("left nodes node {}", _node.mass)
                ret_nodes.append(_node)
                base_name = base_names.get(idx)
                blast_name = blast_names.get(idx)
                #print("add_edge {} {} {} {}".format(node.mass, _node.mass, base_name, base_ppms.get(idx)))
                G.add_edge(node, _node, base=base_name, blast_name=blast_name, ppm=base_ppms.get(idx))
        return ret_nodes

    def best_choice(self, combs):
        df = pd.DataFrame(columns=['Mass', 'Vol', 'PPM', 'Name', 'Adduct', 'Index', 'TheoryMass', 'AdductValue'])
        for comb in combs:
            df = df.append({'Mass': comb[0], 'Vol': comb[1], 'PPM': comb[2], 'Name': comb[3], 'Adduct': comb[4], 'Index': comb[5], 'TheoryMass': comb[6], 'AdductValue': comb[7]}, ignore_index=True)
        df = df.sort_values(by="Vol", ascending=False)
        print(df)
        return df.iloc[0]

    def find_next(self, G, node, theory_mass):
        df = self.df.sort_values(by="Mass", ascending=False)
        df = df[(df.Mass > theory_mass + self.bases.min) & (df.Mass < theory_mass + self.bases.max)]
        #df = df[abs(df.RT - node.rt) < node.width * 10.0]
        '''
            diff_rt = abs(cpd.rt - self.cpd[j].rt)
            if diff_rt > 0 and diff_rt < cpd.width * 5.0:
        '''
        #df = df[]
        #df = df['Mass'] > self.bases.min + theory_mass
        #df = df['Mass'] < self.bases.max + theory_mass
        func = self.bases.match_base
        diff = theory_mass - df.Mass

        #matches = list()
        count = 0
        #min_ppm = 100.0
        #min_idx = -1

        max_vol = 0
        max_idx = -1
        max_base = None
        chose_ppm = 0

        combs = list()
        for idx, row in df.iterrows():
            diff = abs(theory_mass - row['Mass'])
            base, ppm, adduct = func(diff, theory_mass)
            if base:
                count += 1
                #combs.append(("{:.2f}".format(row['Mass']), "{:.2f}".format(row['Vol']), "{:.2f}".format(ppm), base.name))
                combs.append((row['Mass'], row['Vol'], ppm, base.name, adduct, idx, base.mass, self.bases.adds.adduct_value(adduct)))
                #if ppm < min_ppm:
                    #min_ppm = ppm
                    #min_idx = idx
                if row['Vol'] > max_vol:
                    max_vol = row['Vol']
                    max_idx = idx
                    max_base = base
                    chose_ppm = ppm
        #if min_idx == -1:
            #return None
        if max_idx == -1:
            return None, 0

        #logger.info("next_points {}", combs)
        chose_df = self.best_choice(combs)
        row = self.df.loc[chose_df['Index'], :]
        next_node = self.node(G, row['Mass'])
        logger.info("chose mass {} ppm {:.2f}", next_node.mass, chose_ppm)
        try:
            if chose_df.size > 0:
                G.nodes[next_node]['base'] = chose_df['Name']
            #logger.info("set node base {}", max_base)
        except KeyError as ke:
            logger.info("Wrong {}", ke)
        #nx.set_node_attributes(G, 'base', max_base.name)

        #base = G[node][path[idx+1]]['base']
        return next_node, theory_mass + chose_df['TheoryMass'] + chose_df['AdductValue']

    def node(self, G, mass):
        for node in G.nodes():
            if not node:
                continue
            if node.mass == mass:
                return node

    def fill_start(self, row):
        cpd = Compound(row)
        return cpd

    def fill_node(self, node):
        if node.special:
            return None

        row = dict()
        row['Mass'] = node.mass + self.bases.start
        row['RT'] = node.rt
        row['Vol'] = node.vol
        row['Special'] = True
        cpd = Compound(row)
        return cpd

    def find_all_leaves(self, G):
        leaves = list()
        for node in G.nodes():
            if not list(G.successors(node)):
                leaves.append(node)
        all_masses = [node.mass for node in G.nodes()]
        logger.info("G nodes {} masses {}", len(G.nodes()), all_masses)
        #leaves = [node for node in G.nodes() if not G.successors(node)]
        logger.info("G leaves {}", len(leaves))
        return leaves

    def path_avg_vol(self, path):
        if not path:
            return 0

        vol_sum = 0
        for node in path:
            vol_sum += node.vol
        return vol_sum/len(path)
    
    def path_avg_score(self, path):
        if not path:
            return 0

        score_sum = 0
        for node in path:
            score_sum += node.score
        return score_sum/len(path)

    def rank_paths_byvol(self, paths):
        data = []
        for idx, path in enumerate(paths):
            vol_sum = 0
            score_sum = 0
            for node in path:
                vol_sum += node.vol
                score_sum += node.score
            data.append([idx, vol_sum/len(path), score_sum/len(path)])

        df = pd.DataFrame(data, columns=['Index', 'Vol', 'Score'])
        df = df.sort_values(by='Vol', ascending=False)
        #df = df.sort_values(by='Score', ascending=False)
        rank_list = list(df.loc[:,'Index'])
        new_paths = list()
        for i, row in df.iterrows():
            idx = int(row['Index'])
            new_paths.append(paths[idx])
        logger.info("paths count {}", len(paths))
        return new_paths

    def rank_paths(self, paths):
        paths.sort(key=lambda path: len(path), reverse=True)
        if not paths:
            return []
        mid_len = len(paths[0]) - 3
        paths = [path for path in paths if len(path) > mid_len]
        data = []
        for idx, path in enumerate(paths):
            vol_sum = 0
            score_sum = 0
            for node in path:
                vol_sum += node.vol
                score_sum += node.score
            data.append([idx, vol_sum/len(path), score_sum/len(path)])

        df = pd.DataFrame(data, columns=['Index', 'Vol', 'Score'])
        df = df.sort_values(by='Vol', ascending=False)
        #df = df.sort_values(by='Score', ascending=False)
        rank_list = list(df.loc[:,'Index'])
        new_paths = list()
        j = 0
        for i, row in df.iterrows():
            idx = int(row['Index'])
            seq = [node.mass for node in paths[idx]]
            j += 1
            new_paths.append(paths[idx])
        logger.info("paths count {}", len(paths))
        return new_paths

    def paths_to_leaves(self, G, beg_node=None, leaves=None):
        if not leaves:
            leaves = self.find_all_leaves(G)
        start = None
        if beg_node:
            start = beg_node
        else:
            start = self.find_start_point(G)
        paths = list()
        for leave in leaves:
            try:
                paths.extend(nx.all_simple_paths(G, start, leave))
            except:
                logger.warning("no paths between {} and {}", start, leave)

        return paths

    def sequences(self, G, paths):
        seqs = list()
        #for path in paths[:200]:
        for path in paths:
            seq = list()
            for idx in range(len(path)-1):
                node1 = path[idx]
                node2 = path[idx+1]
                base = G[node1][node2]['base']
                seq.append(base)

            seqs.append(''.join(seq))
            print("\nSequence: {}\nindex\tMass\t\tRT\tVol\t\tBase".format(''.join(seq)))
            for idx, node in enumerate(path):
                print("{}\t{:.4f}\t{:.2f}\t{:.2f}\t{}".format(idx+1, node.mass, node.rt, node.vol, seq[idx-1] if idx > 0 and idx <= len(seq) else ''))
            avg_vol = self.path_avg_vol(path)
            avg_score = self.path_avg_score(path)
        seqs = list(set(seqs))
        logger.info("sequences len {} {}", len(seqs), seqs)
        #logger.info("sequences len {}", len(seqs))
        return seqs

    def output_seqs(self, file_path, seqs):
        with open(file_path, 'w+') as f:
            for idx, seq in enumerate(seqs):
                info = ">sequence {}".format(idx)
                f.write(info)
                f.write('\n')
                f.write(seq)
                f.write('\n')

    def attach_ending_edges(self, G):
        nodes = [node for node in G.nodes() if not (nx.is_isolate(G, node) or node.mass < 5000)]
        #nodes = [node for node in G.nodes() if node.mass > 6000]
        ret_nodes = self.find_terminals(G, nodes)
        logger.info("found {}/{} ending nodes", len(ret_nodes), len(nodes))
        return len(ret_nodes)

    def find_ending_nodes(self, G):
        G_terminal = nx.DiGraph()
        G_terminal.add_nodes_from(G.nodes())
        nodes = self.find_terminals(G_terminal, G_terminal.nodes())
        _G = nx.DiGraph()
        _G.add_edges_from(G_terminal.edges())
        return _G

    def walk_multi(self, nodes, G):
        """find decedent nodes, and add edges directly onto G
        """
        #start = self.find_start_point(G)
        #nodes = [start]
        while True:
            nodes = self.find_next_multi(G, nodes)
            if not nodes:
                break

        return G

    def lonely_leaves(self, G):
        leaves = list()
        for node in G.nodes():
            _nodes = list()
            sub_nodes = list(G.successors(node))
            if len(sub_nodes) <= 1:
                continue
            logger.info("find lonely_leaves for node {}", node.mass)
            for sub_node in sub_nodes:
                succs = list(G.successors(sub_node))
                if not succs:
                    _nodes.append(sub_node)
            if len(_nodes) != len(sub_nodes):
                leaves.extend(_nodes)

        return leaves

    def cut_short_branches(self, G, cpd):
        nodes = G.nodes()
        masses = [node.mass for node in nodes]
        masses.sort(reverse=True)
        ret_nodes = list()
        for mass in masses:
            _nodes = list()
            node = self.node(G, mass)
            sub_nodes = list(G.successors(node))
            if len(sub_nodes) <= 1:
                continue
            paths = cpd.paths_to_leaves(G, node)
            if not paths:
                continue
            lens = [len(path) for path in paths]
            logger.info("lengh {}", lens)
            max_len = max(lens)
            if max_len > 6:
                continue
            for path in paths:
                if len(path) < max_len:
                    _nodes = list(path[1:])
                    ret_nodes.extend(_nodes)
            G.remove_nodes_from(_nodes)
        return ret_nodes

    def walk(self, G):
        start = self.find_start_point(G)
        paths = list()
        seqs = list()
        node = start
        logger.info("wlak from first node {}", node)
        seqs.append(node)
        count = 0
        theory_mass = node.mass
        while True:
            node, theory_mass = self.find_next(G, node, theory_mass)
            '''
            if not node:
                if count > 5:
                    break
                # process Cm/Gm
                node = self.fill_node(seqs[-1])
                if not node:
                    break
            '''
            if not node:
                break
            seqs.append(node)
        paths.append(seqs)
        #logger.info("seq len {} {}", len(seqs), seqs)
        return paths

    def process(self, G):
        starts = self.find_start_edge()
        paths = list()
        for idx, start in enumerate(starts):
            seqs = list()
            seqs.append(start[0])
            seqs.append(start[1])
            node = start[1]
            '''
            while True:
                node, theory_mass = self.find_next(G, node, theory_mass)
                if not node:
                    # process Cm/Gm
                    node = self.fill_node(seqs[-1])
                    if not node:
                        break
                seqs.append(node)
            logger.info("seqs idx {} len {}", idx, len(seqs))
            '''
            paths.append(seqs)

        return paths

@click.command()
@click.option('--label', '-l', type=float, help='the label value')
@click.option('--dataset', '-d', type=str, help='the data set')
@click.option('--orientation', '-o', type=int, help='reading direction')
def main(dataset, label, orientation):
    args = sys.argv
    '''
    csv_file = None
    if len(args) == 2:
        csv_file = args[1]
    '''
    csv_file = dataset
    plt.figure(figsize=(100, 75))
    #plt.figure(figsize=(8, 6))
    ax = plt.subplot()
    if not csv_file:
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        csv_file = os.path.join(cur_dir, "statics/compounds.csv")
    file_name = csv_file[:csv_file.rfind('.')]
    graph_file = "{}.png".format(file_name)
    fasta_file = "{}.txt".format(file_name)
    if orientation == 5:
        label = 593.1743
    elif orientation == 3:
        label = 694.2397
    elif orientation == 55:
        label = 693.287
    cpd = Compounds(csv_file=csv_file, label=label)
    G = cpd.generate_graph()
    '''
    paths = cpd.walk(G)
    #seq = self.sequences_of_path(G, path)
    cpd.draw_scatters(ax)
    cpd.draw_paths(paths, ax, G)
    plt.savefig(graph_file)
    plt.show()
    '''
    start = cpd.find_start_point(G)
    G_terminal = cpd.find_ending_nodes(G)
    G = cpd.walk_multi([start], G)
    _G = nx.DiGraph()
    _G.add_edges_from(G.edges())
    #leaves = cpd.lonely_leaves(_G)
    #_G.remove_nodes_from(leaves)
    #G.remove_nodes_from(leaves)
    #short_branch_nodes = cpd.cut_short_branches(_G, cpd)
    #G.remove_nodes_from(short_branch_nodes)
    terminals = G_terminal.nodes()
    #print(G_terminal.edges())
    #_G.add_edges_from(G_terminal.edges())
    #paths = cpd.paths_to_leaves(_G, leaves=terminals)
    paths = cpd.paths_to_leaves(_G)
    paths = cpd.rank_paths(paths)
    seqs = cpd.sequences(G, paths)
    cpd.output_seqs(fasta_file, seqs)
    _G = _G.to_undirected()
    #conns = nx.connected_components(_G)
    #for conn in conns:
        #print(conn)
    labels = dict()
    for node in _G.nodes():
        labels[node] = "{} {} {:.1f}".format(node.mass, int(node.vol), node.score)
    edge_bases = nx.get_edge_attributes(G, 'base')
    edge_ppms = nx.get_edge_attributes(G, 'ppm')
    edge_labels = dict()
    for k, v in edge_ppms.items():
        edge_labels[k] = "{} {}".format(edge_bases.get(k), v)
    G_draw = _G
    pos = nx.spring_layout(G_draw)
    #nx.draw(G_draw, pos=pos, with_label=False)
    nx.draw(_G, pos=pos, with_label=True, labels=labels)
    nx.draw_networkx_edge_labels(G_draw,pos,edge_labels)
    #paths, seqs = cpd.top_n_longest_path(G, 1)
    #cpd.draw_scatters(ax)
    #cpd.draw_paths(paths, ax, G)
    plt.savefig(graph_file)
    #plt.show()
    '''
    paths = cpd.process(G)
    cpd.draw_scatters()
    cpd.draw_paths(paths)
    plt.savefig(graph_file)
    edges = cpd.generate_edges()
    G = cpd.made_edges(G, edges)
    cpd.check_graph_edges(G)
    #G = cpd.clean_invalid_edges(G)
    leads, tails = cpd.lead_tail_nodes(G)
    #paths = cpd.longest_paths(G, leads, tails)
    #path = cpd.longest_path(G)
    #paths = [path]
    paths, seqs = cpd.top_n_longest_path(G, 16)
    new_paths = cpd.conjunct_paths(paths)
    cpd.draw_paths(new_paths)
    #paths = cpd.dfs_edges(G)
    cpd.draw_scatters()
    #cpd.draw_paths(paths)
    cpd.print_sequences(seqs)
    #cpd.draw_lines(paths)
    plt.savefig(graph_file)
    '''

if __name__ == '__main__':
    main()
