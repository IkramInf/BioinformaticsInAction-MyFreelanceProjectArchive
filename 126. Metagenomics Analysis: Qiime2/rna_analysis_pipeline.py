#!/usr/bin/env python
# coding: utf-8

"""
1 Processing of the raw sequence data
    1.1 Filter and trim
    1.2 Learn the error rates
    1.3 Sample inference
    1.4 Merge paired reads
    1.5 Construct sequence table and remove chimeras
    1.6 Track reads through the pipeline
    1.7 Assign taxonomy
    1.8 Renaming the samples and bringing them in the right order
2 Bacterial community composition (BCC)
    2.1 Visualise alpha-diversity
    2.2 BCC on Phylum-, Class-, Order-, Family- and Genus-level
"""

import os  # Operating system related functions
import glob  # File path pattern matching
import wget  # Download files from the internet
import argparse # Command line argument(s) parser
import qiime2  # QIIME 2 framework for microbiome analysis
import subprocess  # Run shell commands from Python
import pandas as pd  # Data manipulation and analysis library
from qiime2.plugins import demux  # Plugin for demultiplexing sequence data
from qiime2.plugins import dada2  # Plugin for denoising and quality control of sequence data
from qiime2.plugins import phylogeny  # Plugin for creating phylogenetic trees
from qiime2.plugins import empress  # Plugin for visualizing and interpreting PCoA plots
from qiime2.plugins import feature_table  # Plugin for working with feature tables (biom tables)
from qiime2.plugins import diversity  # Plugin for calculating various measures of diversity
from qiime2.plugins import feature_classifier  # Plugin for taxonomic classification of sequences
from qiime2.plugins import taxa  # Plugin for working with taxonomic annotations
from qiime2.plugins import composition  # Plugin for analyzing compositionality of microbial communities


def get_paired_end_manifest_file(data_dir):
    """
    Generate a paired-end manifest file for DNA sequencing data.

    Args:
        data_dir (str): The directory containing the sequencing data files.

    Returns:
        str: The filepath of the generated paired-end manifest file.
    """
    # Create a directory to store the manifest file
    os.makedirs("input/", exist_ok=True)
    out_filepath = "input/paired_end_manifest.csv"
    
    # Get a list of all fastq files in the provided data directory
    filenames = glob.glob(os.path.join(data_dir, "*/*.fastq*"))
    dataset = [os.path.abspath(file_path) for file_path in filenames]
    
    # Create a dictionary to store sample IDs and their respective file paths
    path_dict = {}
    for path in dataset:
        filename = os.path.basename(path)
        base_dir = os.path.dirname(path)
        sample_id = filename.strip().split("_")[0]
        path_dict.setdefault(sample_id, []).append(os.path.join(base_dir, filename))

    # Write the manifest file
    with open(out_filepath, "w") as ifile:
        ifile.write("sample-id,absolute-filepath,direction\n")
        for Id, filenames in path_dict.items():
            for filename in filenames:
                if "R1" in os.path.basename(filename).upper():
                    forward = filename
                else:
                    reverse = filename
            base_dir = os.path.dirname(forward)
            ifile.write(f"{Id},{forward},forward\n")
            ifile.write(f"{Id},{reverse},reverse\n")
        
    return out_filepath


def get_classifier(url, out_dir):
    """
    Download a classifier from a URL and save it to the specified directory.
    
    This function creates the output directory if it doesn't exist, and then
    downloads the classifier from the provided URL to that directory.
    
    Parameters:
        url (str): The URL of the classifier to be downloaded.
        out_dir (str): The directory where the downloaded classifier will be saved.
    
    Returns:
        str: The file path to the downloaded classifier.
    
    Example usage:
        classifier_path = get_classifier("https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza", 
        "classifiers/")
    """
    os.makedirs(out_dir, exist_ok=True)
    filename = os.path.basename(url)
    out_filepath = os.path.join(out_dir, filename)
    
    # Check if the file already exists before downloading
    if not os.path.exists(out_filepath):
        wget.download(url=url, out=out_dir)
        print(f"Downloaded classifier from {url} to {out_filepath}")
    else:
        print(f"Classifier already exists at {out_filepath}. Skipping download.")
    
    return out_filepath

def save_results(variable_name, output_dir):
    """
    Save QIIME 2 artifacts and visualizations from a variable containing matrices and visualizations.
    
    This function checks if the provided variable has the _fields attribute, indicating that
    it contains matrices and possibly visualizations. It then attempts to save each matrix
    or visualization in the provided output directory.
    
    Parameters:
        variable_name (object): The variable containing matrices and visualizations to be saved.
        output_dir (str): The directory where the matrices and visualizations will be saved as QIIME 2 artifacts.
    
    Returns:
        None

    Example usage:
        save_results(core_metrics, "output/diversity/")
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if the variable has the _fields attribute
    if hasattr(variable_name, '_fields'):
        # Iterate over each matrix in the _fields attribute
        for matrix_name in variable_name._fields:
            # Get the matrix attribute dynamically using getattr
            matrix = getattr(variable_name, matrix_name)
            # Check if the matrix is an Artifact instance
            if isinstance(matrix, qiime2.Artifact):
                # Save the matrix with the appropriate name
                out_filepath = os.path.join(output_dir, f"{matrix_name}.qza")
                matrix.save(out_filepath)
                print(f"Saved matrix {matrix_name} successfully at {out_filepath}.")
            # Check if the matrix is an Visualization instance
            elif isinstance(matrix, qiime2.Visualization):
                # Save the matrix with the appropriate name
                out_filepath = os.path.join(output_dir, f"{matrix_name}.qzv")
                matrix.save(out_filepath)
                print(f"Saved matrix {matrix_name} successfully at {out_filepath}.")
            else:
                print(f"Matrix {matrix_name} is not a valid Artifact or Visualization. Skipping.")
    else:
        print(f"Save failed! Provided {variable_name} objects could not be saved.")


if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="16S RNA Sequence Analysis.")
    
    # Add arguments for dataset directory and metadata filepath with short options
    parser.add_argument("-d", "--data_dir", required=True, help="Path to the directory containing fastq datasets.")
    parser.add_argument("-m", "--metadata_path", required=True, help="Path to the metadata file.")
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Import data
    paired_end_manifest = get_paired_end_manifest_file(args.data_dir)
    sequences = qiime2.Artifact.import_data('SampleData[PairedEndSequencesWithQuality]',
                                            paired_end_manifest,
                                            view_type='PairedEndFastqManifestPhred33')
    
    # Create a directory to store all output files
    os.makedirs("output/", exist_ok=True)
    sequences.save("output/sequences.qza")

    # Load Metadata
    metadata = qiime2.Metadata.load(args.metadata_path)

    classifier_name = input("Enter a Q2-Feature-Classifier Name.\n[Available Classifiers are: silva_full, silva_region, gg_full, gg_region, silva_weighted_full, gg_weighted_full, gg_weighted_region]: ").strip().lower()

    classifier_urls = {
        "silva_full": "https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza",
        "silva_region": "https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza",
        "gg_full": "https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone_full_length.nb.qza",
        "gg_region": "https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone.v4.nb.qza",
        "silva_weighted_full": "https://data.qiime2.org/2023.5/common/silva-138-99-nb-weighted-classifier.qza",
        "gg_weighted_full": "https://data.qiime2.org/2023.5/common/gg-13-8-99-nb-weighted-classifier.qza",
        "gg_weighted_region": "https://data.qiime2.org/2023.5/common/gg-13-8-99-515-806-nb-weighted-classifier.qza"}

    if classifier_name in classifier_urls:
        classifier_path = get_classifier(classifier_urls[classifier_name], "input/taxonomy_classifier/")
    else:
        print("You Entered a Wrong Classifier Name. Proceeding with Default gg_region (Greengenes 515F/806R classifier).")
        classifier_path = get_classifier(classifier_urls["gg_region"], "input/taxonomy_classifier/")
        
    # Demux summarize
    demux_summary = demux.visualizers.summarize(sequences)
    demux_summary.visualization.save("output/qualities.qzv")

    # ## Denoising and QC filtering
    # DADA2 denoise-paired
    dada2_denoised = dada2.methods.denoise_paired(demultiplexed_seqs=sequences,
                                                  trunc_len_f=285,
                                                  trunc_len_r=240,
                                                  max_ee_f=15,
                                                  max_ee_r=15,
                                                  n_threads=4)

    summary_table = feature_table.visualizers.summarize(table=dada2_denoised.table,
                                             sample_metadata=metadata)
    summary_rep_seq = feature_table.visualizers.tabulate_seqs(data=dada2_denoised.representative_sequences)

    save_results(dada2_denoised, "output/dada2/")

    table_core_features = feature_table.visualizers.core_features(table=dada2_denoised.table)
    table_core_features.visualization

    # ## Build a phylogenetic tree
    # Phylogeny align-to-tree-mafft-fasttree
    representative_sequences = dada2_denoised.representative_sequences
    phylo_tree = phylogeny.pipelines.align_to_tree_mafft_fasttree(representative_sequences)

    save_results(phylo_tree, "output/tree/")

    # Empress tree-plot
    #tree = qiime2.Artifact.load("tree/rooted_tree.qza")
    emp_tree_plot = empress.visualizers.tree_plot(phylo_tree.rooted_tree)
    emp_tree_plot.visualization.save("output/tree/rooted_tree.qzv")

    #jupyter serverextension enable --py qiime2 --sys-prefix

    # Replace 'empress-tree.qzv' with the actual path to your visualization file
    visualization = qiime2.Visualization.load('output/tree/rooted_tree.qzv')
    visualization

    alpha_rarefaction_viz = diversity.visualizers.alpha_rarefaction(
                                table=dada2_denoised.table,
                                max_depth=20,
                                phylogeny=phylo_tree.rooted_tree,
                                metadata=metadata)

    alpha_rarefaction_viz.visualization

    # ## Alpha Diversity
    # Alpha diversity
    table = dada2_denoised.table
    table_summary = feature_table.visualizers.summarize(table, sample_metadata=metadata)
    table_summary.visualization.save("output/dada2/table_summary.qzv")

    # Core metrics phylogenetic
    core_metrics = diversity.pipelines.core_metrics_phylogenetic(
        table=table,
        phylogeny=phylo_tree.rooted_tree,
        sampling_depth=16,
        metadata=metadata)

    # Make sure to replace these with actual variable names and output directory
    save_results(core_metrics, "output/diversity/")

    # Alpha group significance
    alpha_group_significance = diversity.visualizers.alpha_group_significance(
        alpha_diversity=core_metrics.shannon_vector,
        metadata=metadata)

    alpha_group_significance.visualization

    alpha_group_significance_faith = diversity.visualizers.alpha_group_significance(
        alpha_diversity=core_metrics.faith_pd_vector,
        metadata=metadata)

    #alpha_group_significance_even = diversity.visualizers.alpha_group_significance(
    #    alpha_diversity=core_metrics.evenness_vector,
    #    metadata=metadata)

    alpha_group_significance_faith.visualization

    #alpha_group_significance_even.visualization

    # Beta diversity adonis
    adonis_result = diversity.visualizers.adonis(
        distance_matrix=core_metrics.unweighted_unifrac_distance_matrix,
        metadata=metadata,
        formula="SampleName",
        n_jobs=4)

    adonis_result.visualization.save("output/diversity/permanova.qzv")

    # Taxonomy classification
    reads = dada2_denoised.representative_sequences
    # Path to the classifier
    classifier = qiime2.Artifact.load(classifier_path)
    taxa_classified = feature_classifier.methods.classify_sklearn(
        reads=reads,
        classifier=classifier,
        n_jobs=4
    )
    taxa_classified.classification.save("output/taxa.qza")

    # Taxa barplot
    taxonomy = taxa_classified.classification
    taxa_barplot = taxa.visualizers.barplot(
        table=table,
        taxonomy=taxonomy,
        metadata=metadata,
    )
    taxa_barplot.visualization.save("output/taxa_barplot.qzv")

    taxa_barplot.visualization

    # Taxa collapse
    genus_table = taxa.methods.collapse(
        table=table,
        taxonomy=taxonomy,
        level=6
    )
    genus_table.collapsed_table.save("output/genus.qza")

    com_plot_command = [
        'qiime', 'empress', 'community-plot',
        '--i-tree', 'output/tree/rooted_tree.qza',
        '--i-feature-table', 'output/dada2/table.qza',
        '--m-sample-metadata-file', 'metadata.tsv',
        '--m-feature-metadata-file', 'output/taxa.qza',
        '--o-visualization', 'output/community-tree-viz.qzv'
    ]

    subprocess.run(com_plot_command)


    # ## Differential abundance
    table_filter_samples = feature_table.methods.filter_samples(
                                table=dada2_denoised.table,
                                min_frequency=10)

    table_filter_samples.filtered_table.save("output/filtered_table.qza")

    table_filter_samples.filtered_table.view(pd.DataFrame)

    table_rel_features = feature_table.methods.relative_frequency(table=table_filter_samples.filtered_table)

    # Convert the relative frequency table to a frequency table
    #table_freq = table_rel_features.relative_frequency_table

    # Use the frequency table for core features analysis
    core_features = feature_table.visualizers.core_features(
                            table=table_filter_samples.filtered_table)

    comp_pseudo = composition.methods.add_pseudocount(table=dada2_denoised.table)

    comp_pseudo.composition_table.view(pd.DataFrame)

    comp_ancom = composition.visualizers.ancom(
                        table=comp_pseudo.composition_table,
                        metadata=metadata.get_column("SampleName"))

