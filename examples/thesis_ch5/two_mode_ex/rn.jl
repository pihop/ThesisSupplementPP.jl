using DiffEqBiological

rn_full = @reaction_network begin
  (1-z)*rm, x_1 --> x_2
  (1-z)*0.5*rm, x_2 --> x_3
  (1-z)*0.5*rm, x_3 --> x_4
  (1-z)*rm, x_4 --> x_3
  (1-z)*0.5*rm, x_3 --> x_2
  (1-z)*0.5*rm, x_2 --> x_1
  (1-z)*rs, x_4 --> x_4 + z 

  z*rm, x_1 --> x_2
  z*rm, x_2 --> x_3
  z*rm, x_3 --> x_4
end rm rs

rn_mode_1 = @reaction_network begin
  rm, x_1 --> x_2
  0.5*rm, x_2 --> x_3
  0.5*rm, x_3 --> x_4
  rm, x_4 --> x_3
  0.5*rm, x_3 --> x_2
  0.5*rm, x_2 --> x_1
end rm rs

rn_mode_2 = @reaction_network begin
  rm, x_1 --> x_2
  rm, x_2 --> x_3
  rm, x_3 --> x_4
end rm rs

p = (1.0, 0.2)
tspan = (0., 10.)
tplot = 0.0:0.01:10.0
u0 = [0., 1., 0., 0., 0.]
#u0_ = [0., 1., 0., 0.]
scaling = 100.0
u0_ssa = [0., scaling, 0., 0., 0.]
u0_i = [0., scaling, 0., 0.]
