function setup()
% SETUP  

old_dir = cd (fullfile ('..', '..', '..'));
run (fullfile ('solver', 'intlab', 'startintlab.m'));
run (fullfile ('solver', 'sdpt3', 'install_sdpt3.m'));
run (fullfile ('solver', 'sedumi', 'install_sedumi.m'));
addpath (fullfile (pwd (), 'solver', 'sdpa', 'mex'));
addpath (fullfile (pwd (), 'solver', 'csdp', 'matlab'));
run (fullfile ('vsdp', '2018', 'install_vsdp.m'));
cd (old_dir);

end
