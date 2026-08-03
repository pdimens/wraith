import os
import sys

import click
import pysam

from wrath.utils import makewindows

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option('-l', '--auto-sv', is_flag = True, default = False, help = "automatic SV detection")
@click.option('-s', '--start', type=int, help = "start position to subset windows")
@click.option('-e', '--end', type=int,   help = "end position to subset windows")
@click.option('-w', '--window', type=int, default=50000, help = "window size to scan")
@click.option('-o', '--outdir', type=str, default="wrath_out", help = "name of output directory")
@click.option('-p', '--no-plot', is_flag = True, default = False, help = 'Skip heatmap plotting')
@click.option('-t', '--threads',type=int, default=4, help = "threads to use")
@click.option('-v', '--verbose', is_flag = True, default=False, help = "verbose output")
@click.argument('chromosome', required=True, type=str)
@click.argument('reference', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('inputs', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.help_option('--help', hidden = True)
def cli(chromosome, reference, inputs, outdir, window, threads, no_plot, verbose, start, end, auto_sv):
# inputs: list of bam files with paths of the individuals of the population/phenotype of interest
    clength = 0
    if not os.path.isfile(reference + ".fai"):
        pysam.faidx(reference)
    with open(reference + ".fai", 'r') as f:
        for i in f:
            if i.startswith(chromosome):
                clength = int(i.split()[1])
                break
    if clength == 0 :
        print(f"Unable to find chromosome {chromosome} in {reference}")
        sys.exit(1)
    if end:
        clength -= end
    starts,ends = makewindows(0 if not start else start, clength, window)
    os.makedirs(f"{outdir}/beds", exist_ok=True)
    with open(f"wrath_out/beds/{chromosome}.{window}.bed", 'w') as bed:
        for startpos,endpos in zip(starts, ends):
            bed.write(f"{chromosome}\t{startpos}\t{endpos}\n")