#/software/bin/octave
%simulate effects based on the mat file otuput of simulations_run_analysi
arg_list = argv();
directory=arg_list{1};
mat_files=[directory '*.mat']
files=glob(mat_files);
load(mat_files)
num_effcts=100;
num_big_effcts=50;
num_small_effcts=num_effcts-num_big_effcts;
big_effcts=rnd_idx(1:num_big_effcts); 
small_effcts=rnd_idx((num_big_effcts+1):(num_big_effcts+num_small_effcts));
B=zeros(size(M,2),1)
rnd_idx=randperm(dmrs);
B(big_effcts)=normrnd(0,sqrt(0.5/num_effcts),[num_big_effcts,1]);
B(small_effcts)=	normrnd(0,sqrt(0.5/num_effcts),[num_small_effcts,1]);
save('-ascii','simulations/simulations'_B.dat'],'B');

