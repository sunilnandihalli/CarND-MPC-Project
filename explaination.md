The model predictive controller

1. dist = v0*dt + 0.5*a0*dt*dt
   // a0 is assumed to be same as the provided throttle
   // coordinate update equations : 
   psi1 = psi0 - dist*(delta0/Lf)
   if delta0 is non-zero
       r = Lf/(delta0)
       x = x0-r*(sin(psi1)-sin(psi0))
       y = y0-r*(cos(psi0)-cos(psi1))
   else
       x = x0 + dist*cos(psi0)
       y = y0 + dist*sin(psi0)
       
2. The cross-track error was estimated as
  cte = (f(x)-y)/sqrt(1+f'(x)*f'(x))
  In order to avoid having to deal with the jump in the value of the orientation around the 0 degree region, cosine of the angle was used as a proxy for epsi (error in the orientation). This can be calculated as 
  epsi = (f'(x)*cos(psi1)-sin(psi1))/(sqrt(f'(x)*f'(x)+1))
  
3. objective-function and the constraints on various variables to be optimized.
 cte was bound between +/- 2. The objective function would approach infinity as cte goes near +/-2
 epsi was bound between +/- 0.5. This is equivalent to 60 degrees. The objective function would approach infinity near the boundaries.
 v was bounded by a maximum velocity. It objective would approach infinity when near maximum velocity. It also has a reference velocity. The objective function would grow as its value move away from the reference velocity.
 
 
The choice of Timestep Length and Elapsed duration.
  The latency played a major role in the choice of time-step. The elapsed duration was simply 30 times the timestep-length. There were no decision factors involved in deciding the number of timesteps. It was almost arbitrary.
  Since we are able change the throttle or steering angle everytime we get to make the decision, I estimated the amount of elapsed time from my previous control input to the current control-input and used that as the time-step. Because of this, the timestep used during every control-input is different. 
