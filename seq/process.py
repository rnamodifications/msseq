import os
import sys
import copy
import itertools
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.cm as cm
from statsmodels.nonparametric.smoothers_lowess import lowess
import click
from loguru import logger
from compounds import Compound, Compounds
import matplotlib.pyplot as plt

def fail_msg(msg):
    data = dict()
    data['success'] = False
    data['message'] = msg
    return data

def success(d):
    data = dict()
    data['success'] = True
    data['data'] = d
    return data

def action_oneway(dataset_path, orientation, label, write=False):
    handle_5_tail = False
    if orientation == 5:
        handle_5_tail = True

    cpd = Compounds(csv_file=dataset_path)
    G = cpd.generate_graph()

    starts = [cpd.find_start_point(G, label)]
    forward_G = nx.DiGraph()
    forward_G.add_nodes_from(G.nodes())
    forward_G.add_edges_from(G.edges())
    add_nodes(forward_G, starts)
    forward_G, _ = walk_multi(starts, forward_G, cpd, handle_5_tail)

    if write:
        # draw graph into PNG file
        file_name = dataset_path[:dataset_path.rfind('.')]
        graph_file = "{}_forward.png".format(file_name)
        draw_G(forward_G, graph_file)

    output_sequences(forward_G, starts[0])

def ranking(dataset_path, orientation, begin, end, top=100, write=False):
    print("handling dataset {} from orientation {}".format(dataset_path, orientation))
    handle_5_tail = False
    if orientation == 3:
        label = 694.2397
    elif orientation == 5:
        label = 593.1743
        handle_5_tail = True
    cpd = Compounds(csv_file=dataset_path, label=label)
    G = cpd.generate_graph()

    print("processing the forward ladders ...")
    starts = [cpd.find_start_point(G, begin)]
    print('starts {}'.format(starts))
    forward_G = copy.deepcopy(G)
    add_nodes(forward_G, starts)
    forward_G, _ = walk_multi(starts, forward_G, cpd, handle_5_tail)
    leaves = [x for x in forward_G.nodes() if forward_G.out_degree(x)==0 and forward_G.in_degree(x)==1]
    print('{} leaves {}'.format(len(leaves), [leaf.mass for leaf in leaves]))
    if end:
        end_dot = cpd.find_start_point(forward_G, end, fill=False)
        if not end_dot:
            return
        ends = [end_dot]
        paths = paths_to_leaves(forward_G, starts[0], ends)
        if not paths:
            return
        #top_path = top_path_rank_by_vol(paths)
        #ranked_paths = [top_path]
        ranked_paths = paths_rank_by_vol(paths)
        edge_paths, blast_paths = edges_for_paths(forward_G, ranked_paths)
        edge_paths_count = len(edge_paths)
        #final_seqs = [''.join(item) for item in edge_paths]
        final_seqs = [''.join(item) for item in edge_paths[:10]]
        print('top edge_path {}'.format(final_seqs))
    else:
        if not leaves:
            return
        paths = paths_to_leaves(forward_G, starts[0], leaves[-1:])
        if not paths:
            return
        ranked_paths = paths_rank_by_vol(paths)
        edge_paths, blast_paths = edges_for_paths(forward_G, ranked_paths[:top])
        edge_paths_count = len(edge_paths)
        print('got {} edge_paths'.format(edge_paths_count))
        final_seqs = [''.join(item) for item in edge_paths[:10]]
        print('top edge_paths {}'.format(final_seqs))

    if write:
        # draw graph into PNG file
        file_name = dataset_path[:dataset_path.rfind('.')]
        #graph_file = "{}_forward.png".format(file_name)
        #draw_G(forward_G, graph_file)

        # output FASTA file for sequences
        fasta_file = "{}_{}.fasta".format(file_name, int(begin))
        with open(fasta_file, 'w') as f:
            for i in range(min(top, edge_paths_count)):
                path = ranked_paths[i]
                f.write(">{} Begin: {:.4f} End: {:.4f} Count: {} Sum Vols: {:.1f}\n".format(i+1, path[0].mass, path[-1].mass, len(path)-1, sum([node.vol for node in path])))
                #edge = blast_paths[i]
                edge = edge_paths[i]
                f.write(''.join(edge))
                f.write('\n\n')

        if len(ranked_paths) > 0:
            fig = plt.figure(figsize=(16, 12))
            plt.xlabel('Mass (Da)', fontname='Arial', fontsize=15, color='black')
            plt.ylabel('Retention Time (min)', fontname='Arial', fontsize=15, color='black')
            nodes = G.nodes()
            masses = [node.mass for node in nodes]
            rts = [node.rt for node in nodes]
            plt.scatter(masses, rts, color='gray')
            for i in range(min(3, len(ranked_paths))):
                path = ranked_paths[i]
                print([i.mass for i in path])
                plt.plot([i.mass for i in path], [i.rt for i in path], 'gray')
                plt.scatter([i.mass for i in path], [i.rt for i in path], color='black')
                if i == 0:
                    edge = edge_paths[i]
                    for i in range(len(edge)):
                        plt.annotate(s=edge[i], size=15, xy=(path[i].mass, path[i].rt), xytext=(10, 10), ha='center', textcoords='offset points')
            fig.tight_layout()
            png_file = "{}_{}.png".format(file_name, int(begin))
            plt.savefig(png_file, transparent=False, dpi=300)

    body = dict()
    if final_seqs:
        body['sequences'] = final_seqs
    return success(body)

def action(dataset_path, orientation, write=False):
    print("handling dataset {} from orientation {}".format(dataset_path, orientation))
    handle_5_tail = False
    if orientation == 3:
        label = 694.2397
        #rev_starts_val = [347.0631, 323.0519, 363.058, 324.0359]
        rev_starts_val = [18.0106]
        rev_starts_level2 = [628.0932, 629.0772, 630.0612, 652.1044, 653.0884, 668.0993, 669.0833, 676.1156, 692.1105, 708.1054]
    elif orientation == 5:
        label = 593.1743
        #rev_starts_val = [267.0968, 243.0855, 283.0917, 244.0695]
        rev_starts_val = [-61.9557]
        handle_5_tail = True
    elif orientation == 50: # unlabeled 5'
        label = 18.0106
        starts_level2 = [347.0631, 323.0519, 363.058, 324.0359]
        #rev_starts_val = [267.0968, 243.0855, 283.0917, 244.0695]
        rev_starts_val = [-61.9557]
        handle_5_tail = True
    elif orientation == 55: # unlabeled 5' and labeld 5'
        label = 593.1743
        #rev_starts_val = [267.0968, 243.0855, 283.0917, 244.0695]
        rev_starts_val = [-61.9557]
        handle_5_tail = True

        #label = 18.0106
        #starts_level2 = [347.0631, 323.0519, 363.058, 324.0359]
        #rev_starts_val = [267.0968, 243.0855, 283.0917, 244.0695]
        rev_starts_val = [18.0106]
        handle_5_tail = True
    elif orientation == 53:
        label = 852.1934
        rev_starts_val = [267.0968, 243.0855, 283.0917, 244.0695]
        handle_5_tail = True
    else:
        return fail_msg("Orientation is in need to proceed.")

    cpd = Compounds(csv_file=dataset_path, label=label)
    G = cpd.generate_graph()

    print("processing the forward ladders ...")
    starts = [cpd.find_start_point(G, label)]
    forward_G = copy.deepcopy(G)
    add_nodes(forward_G, starts)
    if orientation == 50:
        forward_G, starts_level2 = add_start_edges_level2(forward_G, cpd, [label])
        starts.extend(starts_level2)
    forward_G, _ = walk_multi(starts, forward_G, cpd, handle_5_tail)
    #forward_G = break_crosstalk(forward_G)
    #all_seqs = all_sequences(forward_G, starts[0])
    #print("full sequences {}".format(len(all_seqs)))

    print("processing the backward ladders ...")
    backward_G = copy.deepcopy(G)
    rev_starts = cpd.create_new_points(backward_G, rev_starts_val)
    rev_starts_ori = rev_starts.copy()
    add_nodes(backward_G, rev_starts)
    backward_G, starts_level2_nodes = add_start_edges_level2(backward_G, cpd, rev_starts_val)
    rev_starts.extend(starts_level2_nodes)
    backward_G, _ = walk_multi(rev_starts, backward_G, cpd, not handle_5_tail)
    term_nodes = [node for node in backward_G.nodes() if backward_G.nodes[node].get('terminal')]
    #all_seqs = cpd.paths_to_leaves(backward_G, beg_node=rev_starts[0], leaves=term_nodes)
    #print("full sequences {}".format(len(all_seqs)))

    #print("probing forward_G...")
    #probing_G(forward_G)
    #print("probing backward_G...")
    #probing_G(backward_G)
    print("processing the sequences ...")
    if orientation == 3:
        stack_nodes, fragments = identify_3p_sequences_v2(forward_G, backward_G, starts)
    else:
        stack_nodes, fragments = identify_3p_sequences_v2(backward_G, forward_G, rev_starts_ori)

    final_seqs = format_output_sequences(stack_nodes, fragments)
    #final_G, edges_2b_recycle = handle_final_seqs(forward_G, cpd, stack_nodes)
    #output_sequences(final_G, starts[0])
    #mass_ladders = [node.mass_ladder for node in stack_nodes]
    #result_ladders = extract_sequences(mass_ladders, edges_2b_recycle)
    #print("{} ladders {}".format(len(result_ladders), result_ladders))
    if write:
        print("preparing the sequences output ...")
        # draw graph into PNG file
        file_name = dataset_path[:dataset_path.rfind('.')]
        graph_file = "{}_forward.png".format(file_name)
        draw_G(forward_G, graph_file)
        graph_file = "{}_backward.png".format(file_name)
        backward_G = backward_G.reverse()
        draw_G(backward_G, graph_file)
        graph_file = "{}_result.png".format(file_name)
        draw_G(final_G, graph_file)

    body = dict()
    if final_seqs:
        body['sequences'] = final_seqs
    return success(body)

def probing_G(G):
    masses = [node.mass for node in G.nodes()]
    masses.sort()
    print("masses {} edges {}".format(masses, len(G.edges())))

def add_nodes(G, points):
    print("add_nodes {}".format([p.mass for p in points]))
    G.add_nodes_from(points)

def handle_start(G, cpd, starts):
    def get_base(mass):
        # 324.0359, 323.0519, 347.0631, 363.058
        internal = {'A': 347.0631, 'C': 323.0519, 'G': 363.058, 'U': 324.0359}
        for b in internal.keys():
            if internal[b] == mass:
                return b

    orign = cpd.find_start_point(G, 0)
    for idx, start in enumerate(starts):
        G.add_edge(orign, start)
        b = get_base(start.mass)
        G[orign][start].update({'base': b, 'ppm': 0.05})

    return G

def add_start_edges_level2(G, cpd, starts):
    internal = {'A': 329.0525, 'C': 305.0413, 'G': 345.0474, 'U': 306.0253}
    nodes = list()
    for start in starts:
        node1 = cpd.node(G, start)
        for char in 'ACGU':
            mass = round(start + internal[char], 4)
            node2 = cpd.node(G, mass)
            if not node2:
                node2 = cpd.create_new_point(G, mass)
            nodes.append(node2)
            G.add_edge(node1, node2)
            G[node1][node2].update({'base': char, 'ppm': 0.05})

    nodes = list(set(nodes))
    return G, nodes

def paths_to_leaves(G, start, leaves):
    print('paths_to_leaves leaves {}'.format(leaves))
    paths = nx.all_simple_paths(G, start, leaves)
    paths = list(paths)
    if not paths:
        return []
    path_lens = [len(path) for path in paths]
    max_path_len = max(path_lens)
    long_paths = [path for path in paths if len(path) > max_path_len - 2]

    return long_paths

def paths_rank_by_vol(paths):
    paths.sort(key=lambda path: sum([node.vol for node in path])/len(path), reverse=True)
    return paths

def top_path_rank_by_vol(paths):
    path = max(paths, key=lambda path: sum([node.vol for node in path])/len(path))
    return path

def edges_for_paths(G, paths):
    edges = list()
    blast_edges = list()
    for path in map(nx.utils.pairwise, paths):
        path = list(path)
        base_path = [G[edge[0]][edge[1]]['base'] for edge in path]
        blast_path = [G[edge[0]][edge[1]]['blast_name'] for edge in path]
        edges.append(base_path)
        blast_edges.append(blast_path)

    return edges, blast_edges

def print_sequences(G, paths):
    seqs = list()
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
    seqs = list(set(seqs))
    print("\n{} sequence(s) {}".format(len(seqs), seqs))
    return seqs

def all_sequences(G, start):
    ends = [node for node in G.nodes() if len(list(G.successors(node))) == 0]
    paths = paths_to_leaves(G, start, ends)
    return paths

def output_sequences(G, start):
    paths = all_sequences(G, start)

    print_sequences(G, paths)
    return paths

class StackNode:
    def __init__(self, node, seq, mass_ladder):
        if not isinstance(seq, list):
            return None
        if not isinstance(mass_ladder, list):
            return None
        self.node = node
        self.seq = seq
        self.mass = node.mass
        self.mass_ladder = mass_ladder.copy()
        self.mass_ladder.append(node.mass)

    def __expr__(self):
        msg = "seq: {}\tmass: {}".format(self.seq, self.mass)

def identify_3p_sequences_v2(forward_G, backward_G, starts):
    """identify the sequences shared by both Graphs
    """
    backward_G = backward_G.reverse()
    #backward_G_starts = [node for node in backward_G.nodes() if len(list(backward_G.predecessors(node))) == 0 and backward_G.nodes[node].get('terminal')]
    backward_G_starts = [node for node in backward_G.nodes() if backward_G.nodes[node].get('terminal')]

    def successors_of_node(G, node):
        return list(G.successors(node))

    def edge_of_node_pair(G, node_a, node_b):
        return G.edges[node_a, node_b]

    def edges_of_node(G, node):
        if node not in G.nodes():
            return []
        successors = list(G.successors(node))
        return [G.edges[node, suc] for suc in successors]

    def edge_names_of_node(G, node):
        edges = edges_of_node(G, node)
        return [edge.get('base') for edge in edges]

    def common_edge_names(names_F, names_B):
        return [n for n in names_F if n in names_B]

    def edge_name_of_nodes(G, node_a, node_b):
        edge = G.edges[node_a, node_b]
        if edge:
            return edge.get('base')

    def nodes_with_edge_name(G, node, name):
        nodes = nodes_with_edge_names(G, node, [name])
        return nodes

    def nodes_with_edge_names(G, node, names):
        successors = successors_of_node(G, node)
        nodes = [n for n in successors if edge_name_of_nodes(G, node, n) in names]
        return nodes

    stack_nodes = list()
    fragments = list()
    final_starts = starts
    logger.info("forward_G starts {} backward_G_starts {}", [cpd.mass for cpd in final_starts], [(cpd.mass, cpd.vol) for cpd in backward_G_starts])
    logger.info("final_starts {}", starts)
    for i in range(len(backward_G_starts)):
        stack_F = list()
        stack_B = list()
        stack_F.append(StackNode(final_starts[0], [''], [0.0]))
        stack_B.append(StackNode(backward_G_starts[i], [''], [0.0]))
        logger.info("initial stacks F {} {}", stack_F, final_starts[0])
        logger.info("initial stacks F {} {}", stack_F, backward_G_starts[i])
        while True:
            top_F = stack_F.pop() if stack_F else None
            top_B = stack_B.pop() if stack_B else None
            if not top_F or not top_B:
                break
            print("pop {} from stack F seq {}".format(top_F.mass, ''.join(top_F.seq)))
            print("pop {} from stack B seq {}".format(top_B.mass, ''.join(top_B.seq)))
            if top_B.mass in [0.0, 18.0106, -61.9557]:
                # record top_F.seq
                stack_nodes.append(top_F)
                continue

            edge_names_F = edge_names_of_node(forward_G, top_F.node)
            edge_names_B = edge_names_of_node(backward_G, top_B.node)
            if not edge_names_F or not edge_names_B:
                print("no edge_names_F {} or edge_names_B {}".format(edge_names_F, edge_names_B))
                fragments.append(top_F.seq)
                continue
            comm_edge_names = common_edge_names(edge_names_F, edge_names_B)
            if not comm_edge_names:
                print("no common edge names for F {} B {}".format(edge_names_F, edge_names_B))
                fragments.append(top_F.seq)
                continue

            for comm_edge_name in comm_edge_names:
                common_nodes_F = nodes_with_edge_name(forward_G, top_F.node, comm_edge_name)
                common_nodes_B = nodes_with_edge_name(backward_G, top_B.node, comm_edge_name)
                comm_nodes_num_f = len(common_nodes_F)
                comm_nodes_num_b = len(common_nodes_B)
                if comm_nodes_num_f > 1 or comm_nodes_num_b > 1:
                    common_nodes_F.sort(key=lambda item: item.vol, reverse=True)
                    common_nodes_B.sort(key=lambda item: item.vol, reverse=True)
                for idx, node in enumerate(common_nodes_F):
                    if idx == 0:
                        mass_ladder = list()
                        mass_ladder.extend(top_F.mass_ladder)
                        #mass_ladder.append(top_F.mass)
                        seq = list()
                        seq.extend(top_F.seq)
                        seq.append(comm_edge_name)
                        stack_node_f = StackNode(node, seq, mass_ladder)
                        stack_F.append(stack_node_f)
                        print("append {} to stack F. seq {}".format(stack_node_f.mass, ''.join(stack_node_f.seq)))
                for idx, node in enumerate(common_nodes_B):
                    if idx == 0:
                        mass_ladder = list()
                        mass_ladder.extend(top_B.mass_ladder)
                        #mass_ladder.append(top_B.mass)
                        seq = list()
                        seq.extend(top_B.seq)
                        seq.append(comm_edge_name)
                        stack_node_g = StackNode(node, seq, mass_ladder)
                        stack_B.append(stack_node_g)
                        print("append {} to stack B. seq {}".format(stack_node_g.mass, ''.join(stack_node_g.seq)))

    return stack_nodes, fragments

def break_crosstalk(G):
    """compare two son-nodes, detect if they have common grandson-nodes, remove edges between weak_vol_node and common_grandson_node
    """
    edges_2b_recycle = list()
    for node in G.nodes():
        successors = list(G.successors(node))
        if len(successors) < 2:
            continue
        successors.sort(key=lambda node: node.vol, reverse=True)
        first_son = successors[0]
        grandsons_of_1stson = G.successors(first_son)
        for son in successors[1:]:
            grandsons_of_son = G.successors(son)
            common_grandsons = set(grandsons_of_1stson) & set(grandsons_of_son)
            if len(common_grandsons) >= 1:
                edges_2b_recycle.extend([(son, common_grandson) for common_grandson in common_grandsons])
                print("removing edges between {} and {}".format(son.mass, [node.mass for node in common_grandsons]))

    edges_2b_recycle = list(set(edges_2b_recycle))
    G.remove_edges_from(edges_2b_recycle)
    return G

def extract_sequences(mass_ladders, edges_2b_recycle):
    result_ladders = list()
    masses_2b_recycle = [edge[0].mass for edge in edges_2b_recycle]
    edge_masses = [(edge[0].mass, edge[1].mass) for edge in edges_2b_recycle]
    masses = set(masses_2b_recycle)
    ladders_2b_recycle = list()
    for ladder in mass_ladders:
        for idx, mass in enumerate(ladder):
            if idx > 0 and (ladder[idx-1], ladder[idx]) in edge_masses:
                    # record it
                    print("record {}-{} for ladder {}".format(ladder[idx-1], ladder[idx], ladder))
                    ladders_2b_recycle.append(ladder)
                    break
        #if set(ladder) & masses:
            #continue
        #result_ladders.append(ladder)

    print("mass_ladders {} recycle {}".format(len(mass_ladders), len(ladders_2b_recycle)))
    mass_ladders = [ladder for ladder in mass_ladders if ladder not in ladders_2b_recycle]
    return mass_ladders

def handle_final_seqs(G, cpd, stack_nodes):
    result_G = nx.DiGraph()
    for node in stack_nodes:
        for idx, mass in enumerate(node.mass_ladder):
            if idx == 0:
                continue
            pre_node = cpd.node(G, node.mass_ladder[idx-1])
            if not pre_node:
                continue
            result_G.add_edge(pre_node, cpd.node(G, mass))
            result_G[pre_node][cpd.node(G, mass)].update({'base': node.seq[idx-1], 'ppm': 0})

    """
    nodes_2b_recycle = list()
    for node in result_G.nodes():
        successors = list(result_G.successors(node))
        if len(successors) > 1:
            grandson_nodes = [result_G.successors(node) for node in successors]
            grandson_nodes = list(itertools.chain.from_iterable(grandson_nodes))
            grandson_nodes = list(set(grandson_nodes))
            if len(grandson_nodes) == 1:
                successors.sort(key=lambda node: node.vol, reverse=True)
                nodes_2b_recycle.extend(successors[1:])
                print("removed nodes {} for node {}".format([node.mass for node in successors[1:]], node.mass))

    result_G.remove_nodes_from(nodes_2b_recycle)

    #print("result_G edges {}".format([(edge[0].mass, edge[1].mass) for edge in result_G.edges()]))
    masses_2b_recycle = [node.mass for node in nodes_2b_recycle]
    return result_G, masses_2b_recycle
    """
    edges_2b_recycle = list()
    for node in result_G.nodes():
        successors = list(result_G.successors(node))
        if len(successors) < 2:
            continue
        successors.sort(key=lambda node: node.vol, reverse=True)
        first_son = successors[0]
        grandsons_of_1stson = result_G.successors(first_son)
        for son in successors[1:]:
            grandsons_of_son = result_G.successors(son)
            common_grandsons = set(grandsons_of_1stson) & set(grandsons_of_son)
            if len(common_grandsons) >= 1:
                #edges_2b_recycle.extend([(son, common_grandson) for common_grandson in common_grandsons])
                edges_2b_recycle.extend((node, son))
                print("removing edges between {} and {}".format(son.mass, [node.mass for node in common_grandsons]))
                #if son.mass in [5849.9154, 5520.8994]:
                print("node {} 1stson {} son {} grandsons {}".format(node.mass, first_son.mass, son.mass, [node.mass for node in common_grandsons]))

    edges_2b_recycle = list(set(edges_2b_recycle))
    result_G.remove_edges_from(edges_2b_recycle)
    print("edges_2b_recycle {}".format(edges_2b_recycle))
    return result_G, edges_2b_recycle

def format_output_sequences(stack_nodes, fragments):
    final_sequences = list()
    for node in stack_nodes:
        idxs = range(len(node.seq))
        print("idx {} seq {} mass_ladder {}".format(len(idxs), len(node.seq), len(node.mass_ladder)))
        base_mass_list = list(zip(idxs, node.seq, node.mass_ladder[1:]))
        print(base_mass_list)
        base_mass_list.reverse()
        seq = ''.join(node.seq)
        final_sequences.append(seq)
        print("Sequence {}".format(seq))
        print("Index\tBase\tMass")
        #for idx, item in enumerate(base_mass_list):
        for item in base_mass_list:
            if item[0] == 0:
                continue
            print("{}\t{}\t{:.4f}".format(item[0], item[1], item[2]))

        print("")

    sequences = [(len(node.seq[1:]), ''.join(node.seq)) for node in stack_nodes if len(node.seq) > 1]
    print("{} sequences {}".format(len(sequences), sequences))
    print("fragments {}".format([(len(seq), ''.join(seq)) for seq in fragments if len(seq) > 1]))
    return final_sequences

def walk_multi(starts, _G, cpd, handle_5_tail=False):
    #_G = nx.DiGraph()
    #_G.add_nodes_from(G.nodes())
    #_G.add_edges_from(G.edges())
    #_G.add_nodes_from(starts)
    _G = cpd.walk_multi(starts, _G)

    terminals_num = True
    if handle_5_tail:
        terminals_num = cpd.attach_ending_edges(_G)
        if not terminals_num:
            logger.warning("No termination node found.")
            #return _G, True
        #all_seqs = all_sequences(_G, starts[0])
        #print("full sequences {}".format(len(all_seqs)))
    _G = _refine_nodes_v3(_G)

    if not terminals_num:
        return _G, True
    _G = clear_isolates_nodes(_G)
    logger.info("G edges {}", len(_G.edges()))
    return _G, False

def _refine_nodes(G):
    """remove nearby nodes, keep the node with maximum Volume
    """
    nodes_2b_recycle = list()
    nodes = [node for node in G.nodes() if not nx.is_isolate(G, node)]
    for node in nodes:
        nbs = list(G.successors(node))
        nb_edges = [G.edges[node, nb] for nb in nbs]
        nb_bases = [item.get('base') for item in nb_edges]
        base_values = set(nb_bases)
        nb_bases_groups = [[nb for nb in nbs if G.edges[node, nb]['base'] == val] for val in base_values]
        for group in nb_bases_groups:
            if len(group) == 1:
                continue
            # get top Vol in group, and record other compounds
            group.sort(key=lambda item: item.vol, reverse=True)
            nodes_2b_recycle.extend(group[1:])

    nodes_2b_recycle = list(set(nodes_2b_recycle))
    logger.debug("Recyle {}", [item.mass for item in nodes_2b_recycle])
    G.remove_nodes_from(nodes_2b_recycle)
    return G

def probing_node(G, mass):
    print("trying to probe node for mass {}".format(mass))
    node = None
    for _node in G.nodes():
        if _node.mass == mass:
            node = _node
            break
    if not node:
        return
    successors = G.successors(node)
    predecessors = G.predecessors(node)
    print("probing node {} successors {} predecessors {}".format(node.mass, [node.mass for node in successors], [node.mass for node in predecessors]))

def _refine_nodes_v3(G):
    def rerange_terminal(G, nodes):
        terminal_nodes = [node for node in nodes if G.nodes[node].get('terminal')]
        if terminal_nodes:
            print("===>got terminal_nodes {} from ori {}".format([(node.mass, node.vol) for node in terminal_nodes], [(node.mass, node.vol) for node in nodes]))
            terminal_nodes.reverse()
            for node in terminal_nodes:
                nodes.remove(node)
            for node in terminal_nodes:
                nodes.insert(0, node)
        return nodes

    nodes = [node for node in G.nodes() if not nx.is_isolate(G, node)]

    for node in nodes:
        if not G.has_node(node):
            continue

        nbs = list(G.successors(node))
        nb_edges = [G.edges[node, nb] for nb in nbs]
        nb_bases = [item.get('base') for item in nb_edges]
        base_values = set(nb_bases)
        nb_bases_groups = [[nb for nb in nbs if G.edges[node, nb]['base'] == val] for val in base_values]
        for group in nb_bases_groups:
            if len(group) == 1:
                continue
            # get top Vol in group, and remove other compounds
            group.sort(key=lambda item: item.vol, reverse=True)
            group = rerange_terminal(G, group)
            strong_node = group[0]
            print("removing weak_nodes {} from strong node {}-{} for node {}".format([item.mass for item in group[1:]], group[0].mass, group[0].vol, node.mass))
            for weak_node in group[1:]:
                try:
                    G.remove_node(weak_node)
                except networkx.exception.NetworkXError as error:
                    print(error)
                continue

    return G

def clear_isolates_nodes(G):
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    return G

def draw_G(G, graph_file):
    print("Drawing picture {} ...".format(graph_file))
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning)
    plt.figure(figsize=(100, 75))

    pos = nx.spring_layout(G)
    _draw_nodes(G, pos)
    _draw_edges(G, pos)
    plt.savefig(graph_file)
    print("Done.")

def _draw_nodes(G, pos):
    labels = dict()
    for node in G.nodes():
        if not node:
            continue
        labels[node] = "{} {} {:.1f}".format(node.mass, int(node.vol), node.score)
    nx.draw(G, pos=pos, with_label=True, labels=labels)

def _draw_edges(G, pos):
    edge_bases = nx.get_edge_attributes(G, 'base')
    edge_ppms = nx.get_edge_attributes(G, 'ppm')
    edge_labels = dict()
    for k, v in edge_ppms.items():
        edge_labels[k] = "{} {}".format(edge_bases.get(k), v)
    nx.draw_networkx_edge_labels(G, pos, edge_labels)

@click.command()
@click.option("--dataset", "-d", help="data sourse")
@click.option("--orientation", "-o", help="reading orientation,", type=click.Choice(['3', '5', '53', '50', '55']))
@click.option("--label", "-l", type=click.FLOAT, help="tag value")
@click.option("--begin", "-b", type=float, help="start mass")
@click.option("--end", "-e", type=float, help="end mass")
@click.option("--top", "-t", default=100, type=int, help="top N, default 100")
@click.option("--write", "-w", help="data sourse")
def main(dataset, orientation, label, begin, end, top, write):
    if not dataset:
        return fail_msg("Dataset in need to proceed.")

    if orientation:
        orientation = int(orientation)

    if begin:
        ret = ranking(dataset, orientation, begin, end, top, write)
    elif label:
        logger.info("params: ori {} label {}", orientation, label)
        ret = action_oneway(dataset, orientation, label, write)
    else:
        ret = action(dataset, orientation, write)
        print(ret)

    if ret and not ret['success']:
        logger.error(ret['message'])
        print_help()

def print_help():
    with click.get_current_context() as ctx:
        click.echo(ctx.get_help())

if __name__ == '__main__':
    main()
