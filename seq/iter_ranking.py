import os
import sys
import json
import click
from process import ranking
import multiprocessing as mp

def evaluate_split_dataset(dataset, seq_mass, orientation, non_cds, tag=694.2397):
    if non_cds:
        dataset = os.path.dirname(dataset) + '.xlsx'
    if orientation == 3:
        begin = tag
        end_mass = seq_mass.get('mass')
    else:
        begin = 18.0106
        end_mass = seq_mass.get('mass') + begin - 80 - tag
    ret = ranking(dataset, orientation, begin, end_mass)
    if not ret:
        seqs = None
    else:
        seqs = ret.get('data').get('sequences')

    if seqs:
        seq_mass['output_seq'] = seqs#[0]
    else:
        seq_mass['output_seq'] = []#''

    return seq_mass

def process_split_samples(seq_mass_xlsxs, orientation, non_cds):
    with mp.Pool(mp.cpu_count()) as p:
    #with mp.Pool(4) as p:
            evals = p.starmap(evaluate_split_dataset, [(seq_mass.get('dataset'), seq_mass, orientation, non_cds) for seq_mass in seq_mass_xlsxs])

    print(evals)
    return evals

@click.command()
@click.option("--scheme_file", "-f", type=str, help="scheme file, json")
@click.option("--orientation", "-o", type=int, default=3, help="reading orientation")
@click.option("--non_cds", "-n", type=bool, default=False, help="analysis on original dataset, not the split ones")
def main(scheme_file, orientation, non_cds):
    if not scheme_file:
        print_help()
        return

    with open(scheme_file, 'r') as f:
        seq_mass_xlsxs = json.load(f)

    ret = process_split_samples(seq_mass_xlsxs, orientation, non_cds)

    ext = '_ranking.json'
    if non_cds:
        ext = '_non_cds.json'
    output_file = os.path.splitext(scheme_file)[0] + ext
    with open(output_file, 'w') as f:
        json.dump(ret, f, indent=4)

def print_help():
    with click.get_current_context() as ctx:
        click.echo(ctx.get_help())

if __name__ == '__main__':
    main()
