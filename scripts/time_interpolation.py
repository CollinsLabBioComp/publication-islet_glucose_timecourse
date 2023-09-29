#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import scanpy as sc
import scvelo as scv


def get_weight(distances, index, neighbor):
    return 1 / distances[tuple([index, neighbor])]


def smooth_cell_timepoint(
    ix,
    neighbors,
    time_points,
    distances,
    self_weight
):
    weights = np.array([ get_weight(distances, ix, neighbors[ix][i]) for i in range(len(neighbors[ix])) ])
    sort_indicies = np.argsort(weights)
    
    # The first point is itself --> dist 0 = weight inf. Remove first point to find solely neighbors
    neighbor_timepoints = np.array(time_points[neighbors[ix]].values.astype(float))
    neighbor_timepoints = neighbor_timepoints[sort_indicies][:-1]

    weights = weights[sort_indicies][:-1]
    weights[np.argwhere(np.isinf(weights))] = self_weight

    normalized_weights = (weights / np.sum(weights)) * (1 - self_weight)
    weighted_times = np.multiply(normalized_weights, neighbor_timepoints)
    self_time = time_points[ix]
    self_weighted_time = self_time * self_weight

    return(np.sum(weighted_times) + self_weighted_time)


def compute_smoothed_time(
    adata,
    obs_time_key,
    inter_new_key,
    self_weight = 0.5
):
    neighbor_indices = adata.uns['neighbors']['indices']
    distances = adata.uns['neighbors']['distances']
    time_points = adata.obs[obs_time_key]

    sm_cells = np.array([
        smooth_cell_timepoint(
            i,
            neighbor_indices,
            time_points,
            distances,
            self_weight
        ) for i in range(len(adata.obs))
    ])

    adata.obs[inter_new_key] = (
        (sm_cells - np.nanmin(sm_cells)) / (np.nanmax(sm_cells) - np.nanmin(sm_cells))
    )
    return(adata)


def calculate_dimensional_data(adata, nearest_neighbors = 30):
    sc.pp.pca(adata)

    # Calculate neighbors and UMAP
    # sc.pp.neighbors(adata, n_neighbors = nearest_neighbors, random_state=42)
    scv.pp.neighbors(adata, n_neighbors=nearest_neighbors, random_state=42) # Faster than scanpy
    sc.tl.umap(adata)
    return(adata)


def add_squared_time(adata, pseudotime_key):
    adata.obs[pseudotime_key + "_sq"] = adata.obs[pseudotime_key].values ** 2
    return adata

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="Calculate interpolated time"
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='''
        AnnData input
        '''
    )

    parser.add_argument(
        '-al', '--anndata_layer',
        action='store',
        dest='ann_layer',
        default='cp10k',
        help='''
        AnnData layer to compute PCs on 
        '''
    )

    parser.add_argument(
        '-otc', '--observed_time_column',
        action='store',
        dest='otc',
        default='time_point',
        help='''
        Column in AnnData to use for observed time
        '''
    )

    parser.add_argument(
        '-itc', '--interpolated_time_column',
        action='store',
        dest='itc',
        default='smooth_point',
        help='''
        Column in AnnData to output interpolated time
        '''
    )

    parser.add_argument(
        '-sw', '--self_weight',
        action='store',
        dest='sw',
        type=float,
        default=0.5,
        help='''
        Self weight
        '''
    )

    parser.add_argument(
        '-nn', '--nearest_neighbors',
        action='store',
        dest='nn',
        type=int,
        default=30,
        help='''
        Nearest neighbor
        '''
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='out_file',
        default='output',
        help='Output file'
    )
    options = parser.parse_args()

    # Read anndata
    adata = sc.read_h5ad(options.h5)
    adata.X = adata.layers[options.ann_layer]

    # Calulate interpolated time
    nn = min(options.nn, len(adata.obs))
    adata = calculate_dimensional_data(adata, nn)
    adata = compute_smoothed_time(
        adata,
        options.otc,
        options.itc,
        options.sw
    )

    # Save
    adata.write_h5ad(options.out_file)

if __name__ == '__main__':
    main()