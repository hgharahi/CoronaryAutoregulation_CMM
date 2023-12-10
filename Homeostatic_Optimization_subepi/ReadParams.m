function x = ReadParams(ID, sheet)   

params = xlsread('../CMM_ParamEstimationAll/param_estimation_results_2.xlsx',sheet,[ID,'2:',ID,'23']);

x.rho_w = 1060;

x.Rh = params(1);
x.Hh = params(2);
x.Ph = params(3);

x.c1 = params(4);
x.c2 = params(5);
x.c3 = params(6);
x.c4 = params(7);
x.c5 = params(8);
x.Ghe1 = params(9);
x.Ghe2 = params(10);
x.Ghm = params(11);
x.Ghc = params(12);
x.fiber_ang = params(13);

x.Act_lvl0 = params(14);
x.lmin = params(15);
x.lmax = params(16);

x.e = params(17);
x.c = params(18);
x.m = params(19);

x.Smax = params(20);
x.beta = params(21);
x.p0 = params(22);