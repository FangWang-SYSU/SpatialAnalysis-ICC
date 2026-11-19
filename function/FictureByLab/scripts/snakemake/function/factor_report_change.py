import matplotlib
matplotlib.use("Agg")
import pandas as pd
import argparse
import sys, os, copy, gc, re, gzip, pickle, argparse, logging, warnings
import numpy as np
from jinja2 import Environment, FileSystemLoader
from importlib import resources as impresources

import ficture.scripts


def factor_report_change(_args):

    parser = argparse.ArgumentParser(prog="factor_report_change")
    parser.add_argument('--factor_info', type=str, required=True, help='')
    parser.add_argument('--ref_gene', type=str, required=True, help='')
    parser.add_argument('--outpath', type=str, required=True, help='')
    parser.add_argument('--pref', type=str, help='')
    parser.add_argument('--column_name', type=str, default='TopGene_pval', help='Column name to use (default: TopGene_pval)')
    parser.add_argument('--output_pref', type=str, default='', help='')

    args = parser.parse_args(_args)

    ejs = os.path.join(impresources.files(ficture.scripts), "factor_report.template.html")

    if not os.path.isfile(ejs):
        sys.exit(f"Template file {ejs} not found")


    pref=args.pref
    factor_info=args.factor_info
    ref_gene=args.ref_gene
    column_name=args.column_name
    outpath=args.outpath


    output_pref = args.output_pref
    if output_pref == '':
        output_pref = outpath+"/"+pref

    df = pd.read_csv(factor_info, sep="\t")

    marker = pd.read_csv(ref_gene, sep="\t")
    gene_dict = marker.groupby("cluster").apply(lambda x: dict(zip(x["gene"], x["index"]))).to_dict()

    for cluster, marker_genes in gene_dict.items():
        inter_gene_col = f"{cluster}_Genes" 

        def find_intersection(row):
            if isinstance(row[column_name], str):
                top_genes = set(row[column_name].split(", "))  
                intersection = {gene: marker_genes[gene] for gene in top_genes if gene in marker_genes} 
                sorted_intersection = sorted(intersection.items(), key=lambda x: x[1])
                return ", ".join(f"{gene}({index})" for gene, index in sorted_intersection) if sorted_intersection else " " 
            return ""

        df[inter_gene_col] = df.apply(find_intersection, axis=1)

    f = output_pref+".factor.info.tsv"
    df.to_csv(f, sep='\t', index=False, header=True, float_format="%.5f")


    f = output_pref+".factor.info.tsv"
    with open(f, 'r') as rf:
        lines = rf.readlines()
    header = lines[0].strip().split('\t')
    rows = [ list(enumerate(row.strip().split('\t') )) for row in lines[1:]]

    # Load template
    env = Environment(loader=FileSystemLoader(os.path.dirname(ejs)))
    template = env.get_template(os.path.basename(ejs))


    # Render the HTML file
    html_output = template.render(header=header, rows=rows, image_base64=None, tree_image_alt=None, tree_image_caption=None)

    f=output_pref+".factor.info.html"
    with open(f, "w") as html_file:
        html_file.write(html_output)

if __name__ == "__main__":
    import sys
    factor_report_change(sys.argv[1:])
