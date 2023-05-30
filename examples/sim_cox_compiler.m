appFile = [pwd, '/', 'sim_cox.m'];
additional_files = [pwd, '/../bayes-treesurv-gp/'];
buildResults = compiler.build.standaloneApplication(appFile, ...
    'AdditionalFiles', additional_files);