classdef reuse_struct < handle
   properties
      mode = 'none'; % 'periodic'
      period = 1;
      inited = 0;
      Lc = [];
      Qc = [];
      count = 0;

      RL = [];
      q = [];
      Q = [];
      grad_xy = [];
      div_xy = [];

      fine_recorder = [];

   end
end 


