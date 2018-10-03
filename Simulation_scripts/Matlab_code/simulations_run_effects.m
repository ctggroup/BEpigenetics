#/software/bin/octave
arg_list = argv();
directory=arg_list{1};
rho=str2num(arg_list{2});
mat_files=[directory '*.mat']
files=glob(mat_files);
%we simulate the effects
	%we load the first file to know the dimmensions of the matrices
        eval(['load ' 'simulations/simulations_B.dat'  ]); 
B=simulations_B;
for i=1:length(files)
	sim_name=strsplit(files{i},'.ma')
	sim_name=strsplit(sim_name{1},directory)
	sim_name=[sim_name{2}];
	sim_path= arg_list{1};
	if(exist([sim_path sim_name ])!=7)
	  mkdir([sim_path sim_name ])
	end
	sim_path=[sim_path sim_name]
	eval(['load ' files{i}  ]);


	u2=zscore(M)*B;
	y=u2+normrnd(0,sqrt(0.5),[size(u2,1),1]);
	

	save('-ascii',[sim_path '/' sim_name '_R.dat'] ,'R');
	save('-ascii',[sim_path '/' sim_name '_y.dat'] ,'y');
        save('-ascii',[sim_path '/' sim_name '_M.dat'],'M');
	save('-ascii',[sim_path '/' sim_name '_B.dat'],'B');
	save('-ascii',[sim_path '/' sim_name '_O.dat'],'O');

end
