function [q_incr, gol] = metabolic_control(M)


M_rest = 70; % myocardial oxygen metabolism at rest. The unit is micro-l O2/min/micro-g. This will be used to non-dimensionalize the calcluations here. Strictly, based on Pradhan et al 2016.

[q_rest, gfb1_rest, gol_rest] = feedback_Pradhan2016(M_rest);

M_exercise = M;

[q_exercise, gfb1_exercise, gol_exercise] = feedback_Pradhan2016(M_exercise);

q_incr = eval(q_exercise/q_rest);

gfb = gfb1_exercise/gfb1_rest;

gol = gol_exercise/gol_rest;