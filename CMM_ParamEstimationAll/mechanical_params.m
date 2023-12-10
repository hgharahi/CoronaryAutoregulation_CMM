function x = mechanical_params(subendo, i, F, ID)

x = [subendo.c1(i), subendo.c2(i), subendo.c3(i), subendo.c4(i), subendo.c5(i), ...
    subendo.Ge1(i), subendo.Ge2(i), subendo.Gm(i), subendo.Gc(i), subendo.angle(i), ...
    F*subendo.Smax(i), subendo.l_min(i), subendo.l_max(i)];

% x = [   113.5991,   333.6904,   4.6419, 78.2082,    0.4173, ...
%         1.02,   1.02,   1.05,   .965,  1.063,  ...
%         F*subendo.Sbasal(i),  0.4,    1.7];