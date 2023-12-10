function x = mechanical_params(subendo, i, F, ID)

x = [subendo.c1(i), subendo.c2(i), subendo.c3(i), subendo.c4(i), subendo.c5(i), ...
    subendo.Ge1(i), subendo.Ge2(i), subendo.Gm(i), subendo.Gc(i), subendo.angle(i), ...
    F, subendo.l_min(i), subendo.l_max(i)];
