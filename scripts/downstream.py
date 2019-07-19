#!/usr/bin/env python

import sys
import json
import logging

from os.path import join
from os import makedirs

def main():

    spreadsheet_key = None
    jobs = 11

    with open('/data/input/AppSession.json', 'U') as fd_json:
        app = json.load(fd_json)

    # get command attributes, etc
    for item in app['Properties']['Items']:
        if item['Name'] == 'Input.Projects':
            project_id = item['Items'][0]['Id']
        if item['Name'] == 'Input.rarefaction-depth':
            depth = item['Content']
        if item['Name'] == 'Input.metadata-name':
            metadata_name = item['Content']

    # from BaseSpace's documentation
    input_dir = '/data/input/appresults/'

    base = join('/data/output/appresults/', project_id)
    output_dir = join(base, 'upstream-results')

    makedirs(output_dir, exist_ok=True)

    # TODO: include the path to metadata.tsv
    metadata = join(base, metadata_name)

    # qiime diversity core-metrics-phylogenetic

    # qiime diversity alpha-rarefaction

    # q2-taxa: interactive taxa barplot
    # run summarize_taxa
    #   also requires metdata

    # see https://github.com/biocore/qiime/issues/2034
    if jobs != '1':
        cmd += ' -a -O {jobs}'

    for log_file in glob(join(output_dir, 'log_*')):
        with open(log_file, 'U') as fd_log:
            print fd_log.read()

    return 0


if __name__ == '__main__':
    sys.exit(main())
