# should be updated with every new release
FROM qiime2/core

MAINTAINER Yoshiki Vazquez-Baeza <yoshiki@ucsd.edu>

# TODO: include missing sequence classifiers.
# https://data.qiime2.org/2019.4/common/gg-13-8-99-515-806-nb-classifier.qza
# https://data.qiime2.org/2019.4/common/gg-13-8-99-nb-classifier.qza
# https://data.qiime2.org/2019.4/common/silva-132-99-nb-classifier.qza
# https://data.qiime2.org/2019.4/common/silva-132-99-515-806-nb-classifier.qza

# upstream and downstream scripts to support the BaseSpace apps
ADD upstream.sh /upstream.sh
ADD downstream.sh /downstream.sh
