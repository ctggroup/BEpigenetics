#/software/bin/octave
arg_list = argv();
directory=arg_list{1};
rho=str2num(arg_list{2});
mat_files=[directory '*.mat']
files=glob(mat_files);

for i=1:length(files)
	sim_name=strsplit(files{i},'.ma')
	sim_name=strsplit(sim_name{1},directory)
	sim_name=[sim_name{2} '_rho_' num2str(rho) ];
	sim_path= arg_list{1};
	if(exist([sim_path sim_name ])!=7)
	  mkdir([sim_path sim_name ])
	end
	sim_path=[sim_path sim_name]
	eval(['load ' files{i}  ]);
	%num_effcts=ceil(length(dmrs)*0.10 );
	%num_big_effcts=ceil(num_effcts*0.10);
	%num_small_effcts=num_effcts-num_big_effcts;
	num_effcts=100;
        num_big_effcts=50;
        num_small_effcts=num_effcts-num_big_effcts;
	rnd_idx=randperm(dmrs);
	big_effcts=rnd_idx(1:num_big_effcts); 

	small_effcts=rnd_idx((num_big_effcts+1):(num_big_effcts+num_small_effcts));

	B=zeros(size(M,2),1)
	B(big_effcts)=normrnd(0,sqrt(0.5/num_effcts),[num_big_effcts,1]);
	B(small_effcts)=	normrnd(0,sqrt(0.5/num_effcts),[num_small_effcts,1]);
	u2=zscore(M)*B;
	y=u2+normrnd(0,sqrt(0.5),[size(u2,1),1]);
        u2ctr=bsxfun(@minus, u2, mean(u2));
	theta=acos(rho);
	[Q,Rt]=qr(u2ctr);
	P=Q'*Q;
        u1=normrnd(0,sqrt(0.1),[size(u2,1),1]);
        u1=(eye(size(y,1))-P)*u1;
	temp1=u1./sqrt(sum(u1.^2));
	temp2=u2ctr./sqrt(sum(u2ctr.^2));     
	u1=temp1+(1/tan(theta))*temp2;
	tmpR<-zscore(R)*zscore(R)'*u1/size(R,2)
	tmpM<-zscore(M)*zscore(M)'*u2/size(M,2)
	

	yu=tmpR+tmpM + normrnd(0,1-std(tmpR)-std(tmpM),[size(u2,1),1]);
	save('-ascii',[sim_path '/' sim_name '_R.dat'] ,'R');
	save('-ascii',[sim_path '/' sim_name '_y.dat'] ,'y');
        save('-ascii',[sim_path '/' sim_name '_M.dat'],'M');
	save('-ascii',[sim_path '/' sim_name '_B.dat'],'B');
	save('-ascii',[sim_path '/' sim_name '_O.dat'],'O');
	save('-ascii',[sim_path '/' sim_name '_yu.dat'],'yu');
        save('-ascii',[sim_path '/' sim_name '_u1.dat'],'u1');
	save('-ascii',[sim_path '/' sim_name '_u2.dat'],'u2');

end
