appFile = [pwd, '/', 'pbc_analysis.m'];
additional_files = [pwd, '/../bayes-treesurv-gp/'];
buildResults = compiler.build.standaloneApplication(appFile, ...
    'AdditionalFiles', additional_files);