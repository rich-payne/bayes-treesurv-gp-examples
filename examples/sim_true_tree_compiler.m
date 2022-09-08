appFile = [pwd, '/', 'sim_true_tree.m'];
%additional_files = genpath(strcat([pwd, '/../bayes-treesurv-gp/']));
additional_files = [pwd, '/../bayes-treesurv-gp/'];
buildResults = compiler.build.standaloneApplication(appFile, ...
    'AdditionalFiles', additional_files);
% compiler.package.installer(buildResults);
%, ...
%    'InstallerName','MyMagicInstaller', ...
%    'RuntimeDelivery','installer');

addpath(genpath('../bayes-treesurv-gp/'));