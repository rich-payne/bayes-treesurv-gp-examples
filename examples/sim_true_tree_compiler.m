appFile = [pwd, '/', 'sim_true_tree.m'];
additional_files = [pwd, '/../bayes-treesurv-gp/'];
buildResults = compiler.build.standaloneApplication(appFile, ...
    'AdditionalFiles', additional_files);