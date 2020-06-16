% function [model] = CMEC_2_ZC_2_a_running_MCM_syms(f0_user)
function [model] = CMEC_2_ZC_2_a_running_MCM_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms p_y_1_0 p_y_0_1 mu_1_y_1_0 mu_2_y_1_0 mu_3_y_1_0 mu_4_y_1_0 C_1_1_y_1_0 C_1_2_y_1_0 C_1_3_y_1_0 C_1_4_y_1_0 C_2_2_y_1_0 C_2_3_y_1_0 C_2_4_y_1_0 C_3_3_y_1_0 C_3_4_y_1_0 C_4_4_y_1_0 mu_1_y_0_1 mu_2_y_0_1 mu_3_y_0_1 mu_4_y_0_1 C_1_1_y_0_1 C_1_2_y_0_1 C_1_3_y_0_1 C_1_4_y_0_1 C_2_2_y_0_1 C_2_3_y_0_1 C_2_4_y_0_1 C_3_3_y_0_1 C_3_4_y_0_1 C_4_4_y_0_1

x = [
p_y_1_0, p_y_0_1, mu_1_y_1_0, mu_2_y_1_0, mu_3_y_1_0, mu_4_y_1_0, C_1_1_y_1_0, C_1_2_y_1_0, C_1_3_y_1_0, C_1_4_y_1_0, C_2_2_y_1_0, C_2_3_y_1_0, C_2_4_y_1_0, C_3_3_y_1_0, C_3_4_y_1_0, C_4_4_y_1_0, mu_1_y_0_1, mu_2_y_0_1, mu_3_y_0_1, mu_4_y_0_1, C_1_1_y_0_1, C_1_2_y_0_1, C_1_3_y_0_1, C_1_4_y_0_1, C_2_2_y_0_1, C_2_3_y_0_1, C_2_4_y_0_1, C_3_3_y_0_1, C_3_4_y_0_1, C_4_4_y_0_1 ...
];

% PARAMETERS

syms rm rm2 rs r0 

% KAPPA (constant parameters)

syms indmu1 indmu2 indmu3 indmu4 indmu5 indmu6 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 indC11 indC12 indC13 indC14 indC15 indC16 indC17 indC18 indC19 indC20 indC21 kmu01 kmu02 kmu03 kmu04 kmu05 kmu06 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 kC011 kC012 kC013 kC014 kC015 kC016 kC017 kC018 kC019 kC020 kC021 

syms t

p = [rm,rm2,rs,r0];

k = [indmu1,indmu2,indmu3,indmu4,indmu5,indmu6,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,indC11,indC12,indC13,indC14,indC15,indC16,indC17,indC18,indC19,indC20,indC21,kmu01,kmu02,kmu03,kmu04,kmu05,kmu06,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010,kC011,kC012,kC013,kC014,kC015,kC016,kC017,kC018,kC019,kC020,kC021];

if nargin > 0
   f0_user = varargin{1};
   if ~isnumeric(f0_user)
      p_user = setdiff(symvar(f0_user),p);
      % ADDITIONAL PARAMETERS IN INITIAL CONDITIONS
      p = [p,p_user];
   end
	fmu01 = f0_user(1); 
	fmu02 = f0_user(2); 
	fmu03 = f0_user(3); 
	fmu04 = f0_user(4); 
	fmu05 = f0_user(5); 
	fmu06 = f0_user(6); 
	fC01 = f0_user(7); 
	fC02 = f0_user(8); 
	fC03 = f0_user(9); 
	fC04 = f0_user(10); 
	fC05 = f0_user(11); 
	fC06 = f0_user(12); 
	fC07 = f0_user(13); 
	fC08 = f0_user(14); 
	fC09 = f0_user(15); 
	fC010 = f0_user(16); 
	fC011 = f0_user(17); 
	fC012 = f0_user(18); 
	fC013 = f0_user(19); 
	fC014 = f0_user(20); 
	fC015 = f0_user(21); 
	fC016 = f0_user(22); 
	fC017 = f0_user(23); 
	fC018 = f0_user(24); 
	fC019 = f0_user(25); 
	fC020 = f0_user(26); 
	fC021 = f0_user(27); 
else
	fmu01 = 1; 
	fmu02 = 0; 
	fmu03 = 0; 
	fmu04 = r0; 
	fmu05 = 0; 
	fmu06 = 0; 
	fC01 = 0; 
	fC02 = 0; 
	fC03 = 0; 
	fC04 = 0; 
	fC05 = 0; 
	fC06 = 0; 
	fC07 = 0; 
	fC08 = 0; 
	fC09 = 0; 
	fC010 = 0; 
	fC011 = 0; 
	fC012 = 0; 
	fC013 = 0; 
	fC014 = 0; 
	fC015 = 0; 
	fC016 = 0; 
	fC017 = 0; 
	fC018 = 0; 
	fC019 = 0; 
	fC020 = 0; 
	fC021 = 0; 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

f = sym(zeros(size(x)));

f(1) = -mu_4_y_1_0*p_y_1_0*rs;
f(2) = mu_4_y_1_0*p_y_1_0*rs;
f(3) = -p_y_1_0*(C_1_4_y_1_0*rs + mu_1_y_1_0*rm - mu_2_y_1_0*rm2);
f(4) = -p_y_1_0*(C_2_4_y_1_0*rs - mu_1_y_1_0*rm + 2*mu_2_y_1_0*rm2 - mu_3_y_1_0*rm2);
f(5) = -p_y_1_0*(C_3_4_y_1_0*rs - mu_2_y_1_0*rm2 + 2*mu_3_y_1_0*rm2 - mu_4_y_1_0*rm);
f(6) = -p_y_1_0*(C_4_4_y_1_0*rs - mu_3_y_1_0*rm2 + mu_4_y_1_0*rm);
f(7) = p_y_1_0*(2*C_1_2_y_1_0*rm2 - 2*C_1_1_y_1_0*rm + mu_1_y_1_0*rm + mu_2_y_1_0*rm2);
f(8) = -p_y_1_0*(C_1_2_y_1_0*rm - C_1_1_y_1_0*rm + 2*C_1_2_y_1_0*rm2 - C_1_3_y_1_0*rm2 - C_2_2_y_1_0*rm2 + mu_1_y_1_0*rm + mu_2_y_1_0*rm2);
f(9) = p_y_1_0*(C_1_2_y_1_0*rm2 - C_1_3_y_1_0*rm - 2*C_1_3_y_1_0*rm2 + C_1_4_y_1_0*rm + C_2_3_y_1_0*rm2);
f(10) = p_y_1_0*(C_1_3_y_1_0*rm2 - 2*C_1_4_y_1_0*rm + C_2_4_y_1_0*rm2);
f(11) = p_y_1_0*(2*C_1_2_y_1_0*rm - 4*C_2_2_y_1_0*rm2 + 2*C_2_3_y_1_0*rm2 + mu_1_y_1_0*rm + 2*mu_2_y_1_0*rm2 + mu_3_y_1_0*rm2);
f(12) = p_y_1_0*(C_1_3_y_1_0*rm + C_2_2_y_1_0*rm2 - 4*C_2_3_y_1_0*rm2 + C_2_4_y_1_0*rm + C_3_3_y_1_0*rm2 - mu_2_y_1_0*rm2 - mu_3_y_1_0*rm2);
f(13) = p_y_1_0*(C_1_4_y_1_0*rm + C_2_3_y_1_0*rm2 - C_2_4_y_1_0*rm - 2*C_2_4_y_1_0*rm2 + C_3_4_y_1_0*rm2);
f(14) = p_y_1_0*(2*C_2_3_y_1_0*rm2 - 4*C_3_3_y_1_0*rm2 + 2*C_3_4_y_1_0*rm + mu_2_y_1_0*rm2 + 2*mu_3_y_1_0*rm2 + mu_4_y_1_0*rm);
f(15) = -p_y_1_0*(C_3_4_y_1_0*rm - C_3_3_y_1_0*rm2 - C_2_4_y_1_0*rm2 + 2*C_3_4_y_1_0*rm2 - C_4_4_y_1_0*rm + mu_3_y_1_0*rm2 + mu_4_y_1_0*rm);
f(16) = p_y_1_0*(2*C_3_4_y_1_0*rm2 - 2*C_4_4_y_1_0*rm + mu_3_y_1_0*rm2 + mu_4_y_1_0*rm);
f(17) = C_1_4_y_1_0*p_y_1_0*rs - mu_1_y_0_1*p_y_0_1*rm - mu_1_y_0_1*mu_4_y_1_0*p_y_1_0*rs + mu_1_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
f(18) = p_y_0_1*(mu_1_y_0_1*rm - mu_2_y_0_1*rm) + C_2_4_y_1_0*p_y_1_0*rs - mu_2_y_0_1*mu_4_y_1_0*p_y_1_0*rs + mu_2_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
f(19) = p_y_0_1*(mu_2_y_0_1*rm - mu_3_y_0_1*rm) + C_3_4_y_1_0*p_y_1_0*rs - mu_3_y_0_1*mu_4_y_1_0*p_y_1_0*rs + mu_3_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
f(20) = mu_4_y_1_0^2*p_y_1_0*rs + C_4_4_y_1_0*p_y_1_0*rs + mu_3_y_0_1*p_y_0_1*rm - mu_4_y_0_1*mu_4_y_1_0*p_y_1_0*rs;
f(21) = mu_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0)^2 - C_1_4_y_1_0*p_y_1_0*rs*(2*mu_1_y_0_1 - 2*mu_1_y_1_0) - p_y_0_1*(2*C_1_1_y_0_1*rm - mu_1_y_0_1*rm) - C_1_1_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_1_1_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
f(22) = C_1_2_y_1_0*mu_4_y_1_0*p_y_1_0*rs - C_1_4_y_1_0*p_y_1_0*rs*(mu_2_y_0_1 - mu_2_y_1_0) - C_2_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0) - C_1_2_y_0_1*mu_4_y_1_0*p_y_1_0*rs - p_y_0_1*(2*C_1_2_y_0_1*rm - C_1_1_y_0_1*rm + mu_1_y_0_1*rm) + mu_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_2_y_0_1 - mu_2_y_1_0);
f(23) = p_y_0_1*(C_1_2_y_0_1*rm - 2*C_1_3_y_0_1*rm) - C_1_4_y_1_0*p_y_1_0*rs*(mu_3_y_0_1 - mu_3_y_1_0) - C_3_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0) - C_1_3_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_1_3_y_1_0*mu_4_y_1_0*p_y_1_0*rs + mu_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_3_y_0_1 - mu_3_y_1_0);
f(24) = p_y_0_1*(C_1_3_y_0_1*rm - C_1_4_y_0_1*rm) - C_1_4_y_1_0*p_y_1_0*rs*(mu_4_y_0_1 - mu_4_y_1_0) - C_4_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0) - C_1_4_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_1_4_y_1_0*mu_4_y_1_0*p_y_1_0*rs + mu_4_y_1_0*p_y_1_0*rs*(mu_1_y_0_1 - mu_1_y_1_0)*(mu_4_y_0_1 - mu_4_y_1_0);
f(25) = p_y_0_1*(2*C_1_2_y_0_1*rm - 2*C_2_2_y_0_1*rm + mu_1_y_0_1*rm + mu_2_y_0_1*rm) - C_2_4_y_1_0*p_y_1_0*rs*(2*mu_2_y_0_1 - 2*mu_2_y_1_0) + mu_4_y_1_0*p_y_1_0*rs*(mu_2_y_0_1 - mu_2_y_1_0)^2 - C_2_2_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_2_2_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
f(26) = p_y_0_1*(C_1_3_y_0_1*rm + C_2_2_y_0_1*rm - 2*C_2_3_y_0_1*rm - mu_2_y_0_1*rm) - C_2_4_y_1_0*p_y_1_0*rs*(mu_3_y_0_1 - mu_3_y_1_0) - C_3_4_y_1_0*p_y_1_0*rs*(mu_2_y_0_1 - mu_2_y_1_0) - C_2_3_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_2_3_y_1_0*mu_4_y_1_0*p_y_1_0*rs + mu_4_y_1_0*p_y_1_0*rs*(mu_2_y_0_1 - mu_2_y_1_0)*(mu_3_y_0_1 - mu_3_y_1_0);
f(27) = p_y_0_1*(C_1_4_y_0_1*rm + C_2_3_y_0_1*rm - C_2_4_y_0_1*rm) - C_2_4_y_1_0*p_y_1_0*rs*(mu_4_y_0_1 - mu_4_y_1_0) - C_4_4_y_1_0*p_y_1_0*rs*(mu_2_y_0_1 - mu_2_y_1_0) - C_2_4_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_2_4_y_1_0*mu_4_y_1_0*p_y_1_0*rs + mu_4_y_1_0*p_y_1_0*rs*(mu_2_y_0_1 - mu_2_y_1_0)*(mu_4_y_0_1 - mu_4_y_1_0);
f(28) = p_y_0_1*(2*C_2_3_y_0_1*rm - 2*C_3_3_y_0_1*rm + mu_2_y_0_1*rm + mu_3_y_0_1*rm) - C_3_4_y_1_0*p_y_1_0*rs*(2*mu_3_y_0_1 - 2*mu_3_y_1_0) + mu_4_y_1_0*p_y_1_0*rs*(mu_3_y_0_1 - mu_3_y_1_0)^2 - C_3_3_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_3_3_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
f(29) = p_y_0_1*(C_2_4_y_0_1*rm + C_3_3_y_0_1*rm - C_3_4_y_0_1*rm - mu_3_y_0_1*rm) - C_3_4_y_1_0*p_y_1_0*rs*(mu_4_y_0_1 - mu_4_y_1_0) - C_4_4_y_1_0*p_y_1_0*rs*(mu_3_y_0_1 - mu_3_y_1_0) - C_3_4_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_3_4_y_1_0*mu_4_y_1_0*p_y_1_0*rs + mu_4_y_1_0*p_y_1_0*rs*(mu_3_y_0_1 - mu_3_y_1_0)*(mu_4_y_0_1 - mu_4_y_1_0);
f(30) = p_y_0_1*(2*C_3_4_y_0_1*rm + mu_3_y_0_1*rm) - C_4_4_y_1_0*p_y_1_0*rs*(2*mu_4_y_0_1 - 2*mu_4_y_1_0) + mu_4_y_1_0*p_y_1_0*rs*(mu_4_y_0_1 - mu_4_y_1_0)^2 - C_4_4_y_0_1*mu_4_y_1_0*p_y_1_0*rs + C_4_4_y_1_0*mu_4_y_1_0*p_y_1_0*rs;
M = diag(sym([[1], [1], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_1_0], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1], [p_y_0_1]]));
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = 562949953365017/562949953421312;
x0(2) = 1/10000000000;
x0(3) = indmu3*kmu03 - fmu03*(indmu3 - 1);
x0(4) = indmu4*kmu04 - fmu04*(indmu4 - 1);
x0(5) = indmu5*kmu05 - fmu05*(indmu5 - 1);
x0(6) = indmu6*kmu06 - fmu06*(indmu6 - 1);
x0(7) = indC12*kC012 - fC012*(indC12 - 1);
x0(8) = indC13*kC013 - fC013*(indC13 - 1);
x0(9) = indC14*kC014 - fC014*(indC14 - 1);
x0(10) = indC15*kC015 - fC015*(indC15 - 1);
x0(11) = indC16*kC016 - fC016*(indC16 - 1);
x0(12) = indC17*kC017 - fC017*(indC17 - 1);
x0(13) = indC18*kC018 - fC018*(indC18 - 1);
x0(14) = indC19*kC019 - fC019*(indC19 - 1);
x0(15) = indC20*kC020 - fC020*(indC20 - 1);
x0(16) = indC21*kC021 - fC021*(indC21 - 1);
x0(17) = indmu3*kmu03 - fmu03*(indmu3 - 1);
x0(18) = indmu4*kmu04 - fmu04*(indmu4 - 1);
x0(19) = indmu5*kmu05 - fmu05*(indmu5 - 1);
x0(20) = indmu6*kmu06 - fmu06*(indmu6 - 1);
x0(21) = indC12*kC012 - fC012*(indC12 - 1);
x0(22) = indC13*kC013 - fC013*(indC13 - 1);
x0(23) = indC14*kC014 - fC014*(indC14 - 1);
x0(24) = indC15*kC015 - fC015*(indC15 - 1);
x0(25) = indC16*kC016 - fC016*(indC16 - 1);
x0(26) = indC17*kC017 - fC017*(indC17 - 1);
x0(27) = indC18*kC018 - fC018*(indC18 - 1);
x0(28) = indC19*kC019 - fC019*(indC19 - 1);
x0(29) = indC20*kC020 - fC020*(indC20 - 1);
x0(30) = indC21*kC021 - fC021*(indC21 - 1);
dx0 = sym(zeros(size(x)));
dx0(1) = -(562949953365017*rs*(indmu6*kmu06 - fmu06*(indmu6 - 1)))/562949953421312;
dx0(2) = (562949953365017*rs*(indmu6*kmu06 - fmu06*(indmu6 - 1)))/562949953421312;
dx0(3) = rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1)) - rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - rs*(indC15*kC015 - fC015*(indC15 - 1));
dx0(4) = rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - rs*(indC18*kC018 - fC018*(indC18 - 1)) - 2*rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(5) = rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1)) - rs*(indC20*kC020 - fC020*(indC20 - 1)) - 2*rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1)) + rm*(indmu6*kmu06 - fmu06*(indmu6 - 1));
dx0(6) = rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1)) - rs*(indC21*kC021 - fC021*(indC21 - 1)) - rm*(indmu6*kmu06 - fmu06*(indmu6 - 1));
dx0(7) = 2*rm2*(indC13*kC013 - fC013*(indC13 - 1)) - 2*rm*(indC12*kC012 - fC012*(indC12 - 1)) + rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1));
dx0(8) = rm*(indC12*kC012 - fC012*(indC12 - 1)) - rm*(indC13*kC013 - fC013*(indC13 - 1)) - 2*rm2*(indC13*kC013 - fC013*(indC13 - 1)) + rm2*(indC14*kC014 - fC014*(indC14 - 1)) + rm2*(indC16*kC016 - fC016*(indC16 - 1)) - rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1));
dx0(9) = rm2*(indC13*kC013 - fC013*(indC13 - 1)) - rm*(indC14*kC014 - fC014*(indC14 - 1)) - 2*rm2*(indC14*kC014 - fC014*(indC14 - 1)) + rm*(indC15*kC015 - fC015*(indC15 - 1)) + rm2*(indC17*kC017 - fC017*(indC17 - 1));
dx0(10) = rm2*(indC14*kC014 - fC014*(indC14 - 1)) - 2*rm*(indC15*kC015 - fC015*(indC15 - 1)) + rm2*(indC18*kC018 - fC018*(indC18 - 1));
dx0(11) = 2*rm*(indC13*kC013 - fC013*(indC13 - 1)) - 4*rm2*(indC16*kC016 - fC016*(indC16 - 1)) + 2*rm2*(indC17*kC017 - fC017*(indC17 - 1)) + rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + 2*rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(12) = rm*(indC14*kC014 - fC014*(indC14 - 1)) + rm2*(indC16*kC016 - fC016*(indC16 - 1)) - 4*rm2*(indC17*kC017 - fC017*(indC17 - 1)) + rm*(indC18*kC018 - fC018*(indC18 - 1)) + rm2*(indC19*kC019 - fC019*(indC19 - 1)) - rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1)) - rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(13) = rm*(indC15*kC015 - fC015*(indC15 - 1)) + rm2*(indC17*kC017 - fC017*(indC17 - 1)) - rm*(indC18*kC018 - fC018*(indC18 - 1)) - 2*rm2*(indC18*kC018 - fC018*(indC18 - 1)) + rm2*(indC20*kC020 - fC020*(indC20 - 1));
dx0(14) = 2*rm2*(indC17*kC017 - fC017*(indC17 - 1)) - 4*rm2*(indC19*kC019 - fC019*(indC19 - 1)) + 2*rm*(indC20*kC020 - fC020*(indC20 - 1)) + rm2*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + 2*rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1)) + rm*(indmu6*kmu06 - fmu06*(indmu6 - 1));
dx0(15) = rm2*(indC18*kC018 - fC018*(indC18 - 1)) + rm2*(indC19*kC019 - fC019*(indC19 - 1)) - rm*(indC20*kC020 - fC020*(indC20 - 1)) - 2*rm2*(indC20*kC020 - fC020*(indC20 - 1)) + rm*(indC21*kC021 - fC021*(indC21 - 1)) - rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1)) - rm*(indmu6*kmu06 - fmu06*(indmu6 - 1));
dx0(16) = 2*rm2*(indC20*kC020 - fC020*(indC20 - 1)) - 2*rm*(indC21*kC021 - fC021*(indC21 - 1)) + rm2*(indmu5*kmu05 - fmu05*(indmu5 - 1)) + rm*(indmu6*kmu06 - fmu06*(indmu6 - 1));
dx0(17) = (5497558138330244140625*rs*(indC15*kC015 - fC015*(indC15 - 1)))/549755813888 - rm*(indmu3*kmu03 - fmu03*(indmu3 - 1));
dx0(18) = (5497558138330244140625*rs*(indC18*kC018 - fC018*(indC18 - 1)))/549755813888 + rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - rm*(indmu4*kmu04 - fmu04*(indmu4 - 1));
dx0(19) = (5497558138330244140625*rs*(indC20*kC020 - fC020*(indC20 - 1)))/549755813888 + rm*(indmu4*kmu04 - fmu04*(indmu4 - 1)) - rm*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(20) = (5497558138330244140625*rs*(indC21*kC021 - fC021*(indC21 - 1)))/549755813888 + rm*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(21) = rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) - 2*rm*(indC12*kC012 - fC012*(indC12 - 1));
dx0(22) = rm*(indC12*kC012 - fC012*(indC12 - 1)) - 2*rm*(indC13*kC013 - fC013*(indC13 - 1)) - rm*(indmu3*kmu03 - fmu03*(indmu3 - 1));
dx0(23) = rm*(indC13*kC013 - fC013*(indC13 - 1)) - 2*rm*(indC14*kC014 - fC014*(indC14 - 1));
dx0(24) = rm*(indC14*kC014 - fC014*(indC14 - 1)) - rm*(indC15*kC015 - fC015*(indC15 - 1));
dx0(25) = 2*rm*(indC13*kC013 - fC013*(indC13 - 1)) - 2*rm*(indC16*kC016 - fC016*(indC16 - 1)) + rm*(indmu3*kmu03 - fmu03*(indmu3 - 1)) + rm*(indmu4*kmu04 - fmu04*(indmu4 - 1));
dx0(26) = rm*(indC14*kC014 - fC014*(indC14 - 1)) + rm*(indC16*kC016 - fC016*(indC16 - 1)) - 2*rm*(indC17*kC017 - fC017*(indC17 - 1)) - rm*(indmu4*kmu04 - fmu04*(indmu4 - 1));
dx0(27) = rm*(indC15*kC015 - fC015*(indC15 - 1)) + rm*(indC17*kC017 - fC017*(indC17 - 1)) - rm*(indC18*kC018 - fC018*(indC18 - 1));
dx0(28) = 2*rm*(indC17*kC017 - fC017*(indC17 - 1)) - 2*rm*(indC19*kC019 - fC019*(indC19 - 1)) + rm*(indmu4*kmu04 - fmu04*(indmu4 - 1)) + rm*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(29) = rm*(indC18*kC018 - fC018*(indC18 - 1)) + rm*(indC19*kC019 - fC019*(indC19 - 1)) - rm*(indC20*kC020 - fC020*(indC20 - 1)) - rm*(indmu5*kmu05 - fmu05*(indmu5 - 1));
dx0(30) = 2*rm*(indC20*kC020 - fC020*(indC20 - 1)) + rm*(indmu5*kmu05 - fmu05*(indmu5 - 1));

% OBSERVABLES

y = sym(zeros(41,1));

y(1) = p_y_1_0;
y(2) = p_y_0_1;
y(3) = mu_1_y_0_1*p_y_0_1 + mu_1_y_1_0*p_y_1_0;
y(4) = mu_2_y_0_1*p_y_0_1 + mu_2_y_1_0*p_y_1_0;
y(5) = mu_3_y_0_1*p_y_0_1 + mu_3_y_1_0*p_y_1_0;
y(6) = mu_4_y_0_1*p_y_0_1 + mu_4_y_1_0*p_y_1_0;
y(7) = p_y_0_1*p_y_1_0^2 + p_y_1_0*(p_y_1_0 - 1)^2;
y(8) = p_y_0_1*p_y_1_0*(p_y_0_1 + p_y_1_0 - 2);
y(9) = p_y_1_0*(p_y_1_0 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0);
y(10) = p_y_1_0*(p_y_1_0 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0);
y(11) = p_y_1_0*(p_y_1_0 - 1)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0);
y(12) = p_y_1_0*(p_y_1_0 - 1)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0);
y(13) = p_y_0_1^2*p_y_1_0 + p_y_0_1*(p_y_0_1 - 1)^2;
y(14) = p_y_0_1*(p_y_0_1 - 1)*(mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0);
y(15) = p_y_0_1*(p_y_0_1 - 1)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0);
y(16) = p_y_0_1*(p_y_0_1 - 1)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0);
y(17) = p_y_0_1*(p_y_0_1 - 1)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0) + p_y_0_1*p_y_1_0*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0);
y(18) = p_y_0_1*(C_1_1_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_1_1_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2);
y(19) = p_y_0_1*(C_1_2_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_2_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0));
y(20) = p_y_0_1*(C_1_3_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_3_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0));
y(21) = p_y_0_1*(C_1_4_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_4_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0));
y(22) = p_y_0_1*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2);
y(23) = p_y_0_1*(C_2_3_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)) + p_y_1_0*(C_2_3_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0));
y(24) = p_y_0_1*(C_2_4_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)) + p_y_1_0*(C_2_4_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0));
y(25) = p_y_0_1*(C_3_3_y_0_1 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_3_3_y_1_0 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0)^2);
y(26) = p_y_0_1*(C_3_4_y_0_1 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)) + p_y_1_0*(C_3_4_y_1_0 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0));
y(27) = p_y_0_1*(C_4_4_y_0_1 + (mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_4_4_y_1_0 + (mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0)^2);
y(28) = mu_1_y_0_1*p_y_0_1 + mu_1_y_1_0*p_y_1_0;
y(29) = mu_2_y_0_1*p_y_0_1 + mu_2_y_1_0*p_y_1_0;
y(30) = mu_3_y_0_1*p_y_0_1 + mu_3_y_1_0*p_y_1_0;
y(31) = mu_4_y_0_1*p_y_0_1 + mu_4_y_1_0*p_y_1_0;
y(32) = p_y_0_1*(C_1_1_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_1_1_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)^2);
y(33) = p_y_0_1*(C_1_2_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_2_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0));
y(34) = p_y_0_1*(C_1_3_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_3_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0));
y(35) = p_y_0_1*(C_1_4_y_0_1 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_0_1 + mu_1_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)) + p_y_1_0*(C_1_4_y_1_0 + (mu_1_y_0_1*p_y_0_1 - mu_1_y_1_0 + mu_1_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0));
y(36) = p_y_0_1*(C_2_2_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_2_2_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)^2);
y(37) = p_y_0_1*(C_2_3_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)) + p_y_1_0*(C_2_3_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)*(mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0));
y(38) = p_y_0_1*(C_2_4_y_0_1 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_0_1 + mu_2_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)) + p_y_1_0*(C_2_4_y_1_0 + (mu_2_y_0_1*p_y_0_1 - mu_2_y_1_0 + mu_2_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0));
y(39) = p_y_0_1*(C_3_3_y_0_1 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_3_3_y_1_0 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0)^2);
y(40) = p_y_0_1*(C_3_4_y_0_1 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_0_1 + mu_3_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)) + p_y_1_0*(C_3_4_y_1_0 + (mu_3_y_0_1*p_y_0_1 - mu_3_y_1_0 + mu_3_y_1_0*p_y_1_0)*(mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0));
y(41) = p_y_0_1*(C_4_4_y_0_1 + (mu_4_y_0_1*p_y_0_1 - mu_4_y_0_1 + mu_4_y_1_0*p_y_1_0)^2) + p_y_1_0*(C_4_4_y_1_0 + (mu_4_y_0_1*p_y_0_1 - mu_4_y_1_0 + mu_4_y_1_0*p_y_1_0)^2);

% SYSTEM STRUCT

model.sym.nmx = 27;
model.sym.x = x;
model.sym.u = u;
model.sym.f = f;
model.sym.M = M;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.dx0 = dx0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 0;
end