# Generated by packamon: do not edit by hand
# Instead of modifying this file, you can modify a template. See ?init for details.

FROM openanalytics/r-ver:4.2.2

# System libraries (incl. system requirements for R packages)
RUN apt-get update && apt-get install --no-install-recommends -y \
    libssl-dev \
    zlib1g-dev \
    pandoc \
    default-jre-headless \
    make \
    libmagic-dev \
    libnetcdf-dev \
    libxml2-dev \
    libmagick++-dev \
    libicu-dev \
    imagemagick \
    libgmp10 \
    libgmp-dev \
    libcurl4-openssl-dev \
    libpng-dev \
    librsvg2-dev \
    libglpk-dev \
    git-core \
    default-jdk-headless \
    liblzma-dev \
    libpcre3-dev \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "cat(\"local(options(repos = c(BioCsoft = 'https://bioconductor.org/packages/3.16/bioc', BioCann = 'https://bioconductor.org/packages/3.16/data/annotation', BioCexp = 'https://bioconductor.org/packages/3.16/data/experiment', BioCworkflows = 'https://bioconductor.org/packages/3.16/workflows', CRAN = 'https://cloud.r-project.org')))\n\", file = R.home('etc/Rprofile.site'), append = TRUE)"

# install dependencies
RUN R CMD javareconf
RUN R -q -e "install.packages('remotes')" && \
    R -q -e "install.packages('BiocManager'); BiocManager::install()" && \
    R -q -e "remotes::install_version('assertthat', version = '0.2.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('backports', version = '1.4.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('base64enc', version = '0.1-3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('BH', version = '1.81.0-1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('bigmemory.sri', version = '0.1.6', upgrade = FALSE)" && \
    R -q -e "install.packages('BiocGenerics')" && \
    R -q -e "remotes::install_version('BiocManager', version = '1.30.20', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('bitops', version = '1.0-7', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('brew', version = '1.0-8', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('brio', version = '1.1.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('cli', version = '3.6.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('cluster', version = '2.1.4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('codetools', version = '0.2-19', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('colorspace', version = '2.1-0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('commonmark', version = '1.9.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('cpp11', version = '0.4.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('crayon', version = '1.5.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('curl', version = '5.0.0', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('data.table', version = '1.14.8', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('DBI', version = '1.1.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('DEoptimR', version = '1.0-13', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('digest', version = '0.6.31', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('enviPat', version = '2.6', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('evaluate', version = '0.21', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('fansi', version = '1.0.4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('farver', version = '2.1.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('fastmap', version = '1.1.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('fingerprint', version = '3.5.7', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('foreign', version = '0.8-84', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('formatR', version = '1.14', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('Formula', version = '1.2-5', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('fs', version = '1.6.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('futile.options', version = '1.0.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('generics', version = '0.1.3', upgrade = FALSE)" && \
    R -q -e "install.packages('GenomeInfoDbData')" && \
    R -q -e "remotes::install_version('glue', version = '1.6.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('gtools', version = '3.9.4', upgrade = FALSE)" && \
    R -q -e "install.packages('impute')"
RUN R -q -e "remotes::install_version('isoband', version = '0.2.7', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('iterators', version = '1.0.14', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('jsonlite', version = '1.8.4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('KernSmooth', version = '2.23-21', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('labeling', version = '0.4.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('lattice', version = '0.21-8', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('lazyeval', version = '0.2.2', upgrade = FALSE)" && \
    R -q -e "install.packages('limma')" && \
    R -q -e "remotes::install_version('logger', version = '0.2.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('magrittr', version = '2.0.3', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('MALDIquant', version = '1.22.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('MASS', version = '7.3-60', upgrade = FALSE)" && \
    R -q -e "install.packages('MassSpecWavelet')" && \
    R -q -e "remotes::install_version('matrixStats', version = '0.63.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('mime', version = '0.12', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('ncdf4', version = '1.21', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('nnet', version = '7.3-19', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('OrgMassSpecR', version = '0.5-3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('pkgconfig', version = '2.0.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('png', version = '0.1-8', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('pracma', version = '2.4.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('praise', version = '1.0.0', upgrade = FALSE)" && \
    R -q -e "install.packages('preprocessCore')" && \
    R -q -e "install.packages('ProtGenerics')" && \
    R -q -e "remotes::install_version('proxy', version = '0.4-27', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('ps', version = '1.7.5', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('R.methodsS3', version = '1.8.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('R6', version = '2.5.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('RANN', version = '2.6.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rappdirs', version = '0.3.3', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('RColorBrewer', version = '1.1-3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('Rcpp', version = '1.0.10', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('remotes', version = '2.4.2', upgrade = FALSE)" && \
    R -q -e "install.packages('Rhdf5lib')" && \
    R -q -e "remotes::install_version('rJava', version = '1.0-6', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rjson', version = '0.2.21', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rlang', version = '1.1.1', upgrade = FALSE)" && \
    R -q -e "install.packages('RMassBankData')" && \
    R -q -e "remotes::install_version('rpart', version = '4.1.19', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rprojroot', version = '2.0.3', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('rstudioapi', version = '0.14', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rsvg', version = '2.4.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('snow', version = '0.4-4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('sourcetools', version = '0.1.7-1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('stringi', version = '1.7.12', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('sys', version = '3.4.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('utf8', version = '1.2.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('uuid', version = '1.1-0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('viridisLite', version = '0.4.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('withr', version = '2.5.0', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('xfun', version = '0.39', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('XML', version = '3.99-0.14', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('xml2', version = '1.3.4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('xtable', version = '1.8-4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('yaml', version = '2.3.7', upgrade = FALSE)" && \
    R -q -e "install.packages('zlibbioc')" && \
    R -q -e "install.packages('affyio')" && \
    R -q -e "remotes::install_version('askpass', version = '1.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('bigmemory', version = '4.6.1', upgrade = FALSE)" && \
    R -q -e "install.packages('Biobase')"
RUN R -q -e "remotes::install_version('cachem', version = '1.0.8', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('caTools', version = '1.18.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('checkmate', version = '2.2.0', upgrade = FALSE)" && \
    R -q -e "install.packages('ChemmineOB')" && \
    R -q -e "remotes::install_version('clue', version = '0.3-64', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('data.tree', version = '1.0.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('desc', version = '1.4.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('diffobj', version = '0.3.5', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('ellipsis', version = '0.3.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('foreach', version = '1.5.2', upgrade = FALSE)"
RUN R -q -e "install.packages('graph')" && \
    R -q -e "remotes::install_version('highr', version = '0.10', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('itertools', version = '0.1-3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('lambda.r', version = '1.2.4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('later', version = '1.3.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('lifecycle', version = '1.0.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('magick', version = '2.7.4', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('Matrix', version = '1.5-4', upgrade = FALSE)" && \
    R -q -e "install.packages('MatrixGenerics')" && \
    R -q -e "remotes::install_version('munsell', version = '0.5.0', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('nlme', version = '3.1-162', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('plyr', version = '1.8.8', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('processx', version = '3.8.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('R.oo', version = '1.25.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rcdklibs', version = '2.8', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('RCurl', version = '1.98-1.12', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rex', version = '1.2.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('robustbase', version = '0.95-1', upgrade = FALSE)" && \
    R -q -e "install.packages('S4Vectors')" && \
    R -q -e "remotes::install_version('tinytex', version = '0.45', upgrade = FALSE)"
RUN R -q -e "install.packages('affy')" && \
    R -q -e "remotes::install_version('callr', version = '3.7.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('doParallel', version = '1.0.17', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('futile.logger', version = '1.4.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('gplots', version = '3.1.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('gtable', version = '0.3.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('htmltools', version = '0.5.5', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('igraph', version = '1.4.2', upgrade = FALSE)" && \
    R -q -e "install.packages('IRanges')" && \
    R -q -e "remotes::install_version('knitr', version = '1.42', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('memoise', version = '2.0.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('mgcv', version = '1.8-42', upgrade = FALSE)" && \
    R -q -e "install.packages('MsCoreUtils')" && \
    R -q -e "install.packages('mzR')" && \
    R -q -e "remotes::install_version('openssl', version = '2.0.6', upgrade = FALSE)" && \
    R -q -e "install.packages('pcaMethods')" && \
    R -q -e "remotes::install_version('pkgload', version = '1.3.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('promises', version = '1.2.0.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('R.utils', version = '2.12.2', upgrade = FALSE)" && \
    R -q -e "install.packages('RBGL')"
RUN R -q -e "remotes::install_version('rcdk', version = '3.7.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('scales', version = '1.2.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('vctrs', version = '0.6.2', upgrade = FALSE)" && \
    R -q -e "install.packages('BiocParallel')" && \
    R -q -e "remotes::install_version('crosstalk', version = '1.2.0', upgrade = FALSE)" && \
    R -q -e "install.packages('DelayedArray')" && \
    R -q -e "remotes::install_version('fontawesome', version = '0.5.1', upgrade = FALSE)" && \
    R -q -e "install.packages('GenomeInfoDb')" && \
    R -q -e "remotes::install_version('gridExtra', version = '2.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('httpuv', version = '1.6.11', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('httr', version = '1.4.6', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('jquerylib', version = '0.1.4', upgrade = FALSE)" && \
    R -q -e "install.packages('mzID')" && \
    R -q -e "remotes::install_version('pillar', version = '1.9.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('purrr', version = '1.0.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('sass', version = '0.4.6', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('stringr', version = '1.5.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('tidyselect', version = '1.2.0', upgrade = FALSE)" && \
    R -q -e "install.packages('XVector')" && \
    R -q -e "remotes::install_version('bslib', version = '0.4.2', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('covr', version = '3.6.2', upgrade = FALSE)" && \
    R -q -e "install.packages('GenomicRanges')" && \
    R -q -e "remotes::install_version('readJDX', version = '0.6.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('roxygen2', version = '7.2.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('selectr', version = '0.4-2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('tibble', version = '3.2.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('dplyr', version = '1.1.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('ggplot2', version = '3.4.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rematch2', version = '2.1.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('rmarkdown', version = '2.21', upgrade = FALSE)"
RUN R -q -e "remotes::install_version('rvest', version = '1.0.3', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('shiny', version = '1.7.4', upgrade = FALSE)" && \
    R -q -e "install.packages('SummarizedExperiment')" && \
    R -q -e "remotes::install_version('bookdown', version = '0.34', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('htmlwidgets', version = '1.6.2', upgrade = FALSE)" && \
    R -q -e "install.packages('MsFeatures')" && \
    R -q -e "remotes::install_version('shinydashboard', version = '0.7.2', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('tidyr', version = '1.3.0', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('viridis', version = '0.6.3', upgrade = FALSE)" && \
    R -q -e "install.packages('vsn')"
RUN R -q -e "remotes::install_version('waldo', version = '0.5.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('webchem', version = '1.2.0', upgrade = FALSE)" && \
    R -q -e "install.packages('BiocStyle')" && \
    R -q -e "remotes::install_version('DT', version = '0.27', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('htmlTable', version = '2.4.1', upgrade = FALSE)" && \
    R -q -e "install.packages('MSnbase')" && \
    R -q -e "remotes::install_version('plotly', version = '4.10.1', upgrade = FALSE)" && \
    R -q -e "remotes::install_version('testthat', version = '3.1.8', upgrade = FALSE)" && \
    R -q -e "install.packages('ChemmineR')" && \
    R -q -e "remotes::install_version('Hmisc', version = '5.1-0', upgrade = FALSE)"
RUN R -q -e "install.packages('xcms')" && \
    R -q -e "install.packages('CAMERA')" && \
    R -q -e "install.packages('RMassBank')" && \
    R -q -e "remotes::install_github('WMBEdmands/compMS2Miner')" && \
    R -q -e "remotes::install_github('daqana/dqshiny')"


