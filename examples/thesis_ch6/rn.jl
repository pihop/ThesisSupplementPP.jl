using DiffEqBiological

@reaction_func con(r,x) = exp(-(pc*x)/N)

rn_full = @reaction_network rn begin
  (1-z)*rm*con(pc,x_1), x_1 --> x_2
  (1-z)*0.5*rm*con(pc,x_2), x_2 --> x_1
  (1-z)*0.5*rm*con(pc,x_2), x_2 --> x_3
  (1-z)*0.5*rm*con(pc,x_3), x_3 --> x_4
  (1-z)*0.5*rm*con(pc,x_3), x_3 --> x_2
  (1-z)*rm*con(pc,x_4), x_4 --> x_3
  (1-z)*rs, x_4 --> x_4 + z 

  ps*z*rm*con(pc,x_1), x_1 --> x_2
  ps*z*rm*con(pc,x_2), x_2 --> x_3
  ps*z*rm*con(pc,x_3), x_3 --> x_4
  (1-ps)*z*rm*con(pc,x_2), x_2 --> x_1
  (1-ps)*z*rm*con(pc,x_3), x_3 --> x_2 
end rm rs ps pc N 

rn_mode_1 = @reaction_network begin
  rm*con(pc,x_1), x_1 --> x_2
  0.5*rm*con(pc,x_2), x_2 --> x_1
  0.5*rm*con(pc,x_2), x_2 --> x_3
  0.5*rm*con(pc,x_3), x_3 --> x_4
  0.5*rm*con(pc,x_3), x_3 --> x_2
  rm*con(pc,x_4), x_4 --> x_3
end rm rs ps pc N

rn_mode_2 = @reaction_network begin
  ps*rm*con(pc,x_1), x_1 --> x_2
  ps*rm*con(pc,x_2), x_2 --> x_3
  ps*rm*con(pc,x_3), x_3 --> x_4
  (1-ps)*rm*con(pc,x_2), x_2 --> x_1
  (1-ps)*rm*con(pc,x_3), x_3 --> x_2 
end rm rs ps pc N 
