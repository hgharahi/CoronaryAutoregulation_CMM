function H = HtoR(R)

D = 2*R*10^3;
H = (41.1*D + 3.2)/1e6; %% From Gou and Kassab 2003

H = H/R;