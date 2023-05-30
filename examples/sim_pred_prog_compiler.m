appFile = [pwd, '/', 'sim_pred_prog.m'];
additional_files = [pwd, '/../bayes-treesurv-gp/'];
buildResults = compiler.build.standaloneApplication(appFile, ...
    'AdditionalFiles', additional_files);