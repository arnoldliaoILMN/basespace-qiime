#!/usr/bin/env python

import sys
import json
import glob
import numpy as np

import qiime2 as q2
import qiime2.plugins.demux.actions as q2_demux
import qiime2.plugins.feature_table.actions as q2_feature_table
import qiime2.plugins.quality_filter.actions as q2_quality_filter
import qiime2.plugins.deblur.actions as q2_deblur
import qiime2.plugins.dada2.actions as q2_dada2
import qiime2.plugins.phylogeny.actions as q2_phylogeny
import qiime2.plugins.feature_classifier.actions as q2_feature_classifier

from math import floor
from os.path import join
from os import makedirs
from collections import namedtuple


def main():
    with open('/data/input/AppSession.json', 'U') as fd_json:
        app = json.load(fd_json)

    # TODO: this should be a parameter
    deblur = True
    samples = []

    # TODO: confirm the size of the AMI is fixed to 12 cores, if this
    # is the case then lets always split over 11 jobs
    jobs = 11

    # TODO: Making a parser of this whole thing would be better (could also
    # include some unit testing). Also encapsulating everything in one object
    # would be better
    for item in app['Properties']['Items']:
        if item['Name'] == 'Input.Projects':
            project_id = item['Items'][0]['Id']
        # if item['Name'] == 'Input.number-of-jobs':
        #     jobs = int(item['Content'])

        if item['Name'] == 'Input.sample-id':
            for sample in item['Items']:
                s = {'id': None, 'name': None, 'path': None, 'length': None}

                s['id'] = sample['Id']
                s['name'] = sample['Name']
                s['length'] = float(sample['Read1'])

                # path to the samples
                path = join('/data/input/samples/', s['id'])

                # this assumes the per-sample directory structure is
                # small/simple otherwise this might be resource-intensive
                paths = glob.iglob(join(path, '**/*.fastq.gz'), recursive=True)
                for path in paths:
                    if '_R1_' in path or '.R1.' in path:
                        if s['path'] is None:
                            s['path'] = path
                        else:
                            raise ValueError('More than one forward reads have'
                                             ' been found for sample %s' %
                                             sample)

                samples.append(s)

    # from BaseSpace's documentation
    workspace = '/data/scratch/'

    base = join('/data/output/appresults/', project_id)
    output_dir = join(base, 'upstream-results')

    # for sanity
    makedirs(output_dir, exist_ok=True)

    # this file needs to be created based on the sequence files found in the
    # run folder, note that this manifest would only be good for BaseSpace
    # therefore we don't save this in the output folder
    manifest_fp = join(workspace, 'manifest.tsv')
    with open(manifest_fp, 'w') as f:
        f.write('sample-id\tabsolute-filepath\n')
        for sample in samples:
            f.write('{}\t{}\n'.format(sample.id, sample.path))

    # import the sequence data based on the manifest-provided information
    demux = q2.Artifact.import_data('SampleData[SequencesWithQuality]',
                                    manifest_fp,
                                    'SingleEndFastqManifestPhred33V2')

    # based on a conversation with Arnold we determined that a resonable
    # default for sequence trimming would be 70% of the sequence lengths
    trim_length = floor(np.array([s['length'] for s in samples]).mean() * 0.70)

    # create a summary visualization of the quality scores
    summary, = q2_demux.summarize(demux)
    summary.save(join(output_dir, 'sequences.qzv'))

    # quality controlling the sequences, the defaults should be good enough
    filtered, _ = q2_quality_filter.q_score(demux)
    filtered.save(join(output_dir, 'quality-controlled.sequences.qza'))

    # https://developer.basespace.illumina.com/docs/content/documentation/
    # apptools/formbuilder-overview#AdvancedToggling
    # TODO: In the form create a series of parameters that can be modified
    # depending on what method is selected
    if deblur:
        feature_table, representative_sequences, _ = \
            q2_deblur.denoise_16S(demultiplexed_seqs=filtered,
                                  sample_stats=True, trim_length=trim_length,
                                  jobs_to_start=jobs, hashed_feature_ids=True)
    else:
        feature_table, representative_sequences, _ = \
            q2_dada2.denoise_single(demultiplexed_seqs=filtered,
                                    trunc_len=trim_length, n_threads=jobs,
                                    hashed_feature_ids=True)

    # save both outputs
    feature_table.save(join(output_dir, 'feature-table.qza'))
    representative_sequences.save(join(output_dir,
                                       'representative-sequences.qza'))

    # create a summary of the feature table
    summary, = q2_feature_table.summarize(table=feature_table)
    summary.save(join(output_dir, 'feature-table-summary.qzv'))

    # build a *quick* phylogenetic tree using MAFFT. We should build this tree
    # using SEPP, however that would take more resources than what this AMI has
    _, _, _, rooted_tree = q2_phylogeny.align_to_tree_mafft_fasttree(
        sequences=representative_sequences,
        n_threads=jobs)
    rooted_tree.save(join(output_dir, 'rooted-tree.qza'))

    # taxonomic classification
    # TODO: in the form list the four options for the classifiers
    classifier = q2.Artifact.load(join(workspace,
                                       'gg-13-8-99-515-806-nb-classifier.qza'))
    classification, = q2_feature_classifier.classify_sklearn(
            reads=representative_sequences, classifier=classifier,
            reads_per_batch=1000, n_jobs=jobs)
    classification.save(join(output_dir, 'taxonomy.qza'))

    # Generated outputs are:
    #
    # For record keeping and in incase the users would like to process their
    # sequences in a different way:
    # * demultiplexed sequences
    # * sequence quality summary
    # * quality controlled sequences
    #
    # For the qiime downstream app:
    # * feature table (feature-table.qza)
    # * feature table summary per sample & feature (feature-table-summary.qza)
    # * representative sequences (representative-sequences.qza)
    # * rooted tree (rooted-tree.qza)
    # * taxonomic classifications (taxonomy.qza)
    #
    # Ask around to see if it would make sense to save as TSV, newick and FASTA
    # files

    return 0


if __name__ == '__main__':
    sys.exit(main())
