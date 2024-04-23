"""Module to deal with the translation of genomic coordinates.

Originating from a standard annotation file (i.e. GTF) into Hi-C fragments/bins cordinates for
integration to the GREAT framework analysis.
"""

import math
import pandas as pd


def compute_gene_coordinates(dfGenes, dfGenome, resolution):
    """Get two tab files, one with the gene names and their respective start sites and one with the 3D genome coordinates (as midpoints of fragments of fixed length i.e. resolution).
    
    - Args:
    - `starts_sites_Genes`: A gene names and start sites data frame.
    - `dfGenome3D`: The 3D genome coordinates data frame.
    - `resolution`: The fixed length of the fragments used for 3D genome coordinates.
    - Return:
    - A DataFrame with the 3D coordinates of genes.
    """
    # Check chromosome consistancy
    chromsGenome = list(set(dfGenome.loc[:, "chr"].tolist()))
    chromsGenes = list(set(dfGenes.loc[:, "chr"].tolist()))
    if chromsGenes != chromsGenome:
        raise ValueError("Chromosome inconsistency between genes file and genome file. Check the input files.")
    # Prepare the data frame lists
    names = []
    chrs = []
    starts = []
    Xs = []
    Ys = []
    Zs = []
    # Loop through all chromosomes in the files
    for ch in chromsGenes:
        dfChromGene = dfGenes[dfGenes["chr"] == ch]
        dfChromGenome = dfGenome[dfGenome["chr"] == ch]
        for gene in dfChromGene.itertuples():
            for point in dfChromGenome.itertuples():
                if abs(gene.start - point.midpoint) < resolution:
                    pointBefor = point
                    pointAfter = dfChromGenome[dfChromGenome["midpoint"] == point.midpoint + resolution]
                    break
            # A geometric approach to translate gene start bps coordinates to correct proportional point on the genome segment in 3D spce.
            xb, yb, zb = pointBefor.X, pointBefor.Y, pointBefor.Z
            xa, ya, za = pointAfter.iloc[0].X, pointAfter.iloc[0].Y, pointAfter.iloc[0].Z
            distance = math.sqrt((xa - xb)**2 + (ya - yb)**2 + (za - zb)**2)
            a_bNorm = [(xa - xb) / distance, (ya - yb) / distance, (za - zb) / distance]
            geneCoef = ((gene.start - point.midpoint) / resolution) * distance
            geneExtr = [e * geneCoef for e in a_bNorm]
            names.append(gene.Index)
            chrs.append(ch)
            starts.append(gene.start)
            Xs.append(point.X + geneExtr[0])
            Ys.append(point.Y + geneExtr[1])
            Zs.append(point.Z + geneExtr[2])
    # Generate the new data frame and reurn it
    df = pd.DataFrame({"chr": chr, "start": starts, "X": Xs, "Y": Ys, "Z": Zs})
    df.index = names
    return df
