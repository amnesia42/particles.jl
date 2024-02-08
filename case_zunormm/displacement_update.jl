"""
   !f(ds,s,t,i,d)

Dynamic model, computes as ds the function f in the equation ds=f(s,t)dt+g(s,t)dw 
for s at current time t for particle i and possibly using data/functions from d of type userdata.
"""
# add a function f that can apply RK4 
function f!(∂s, s, t, i, d)
   x, y, z, age = s
#    if IS_USE_RK4
#       dt = d["dt"]
#       k1 = u(x,y,z,t)
#       l1 = v(x,y,z,t)
#       m1 = w(x,y,z,t) #+ derivative_K(x, y, z, t)
#       k2 = u(x+1.0/2*dt*k1, y+1.0/2*dt*l1, z+1.0/2*dt*m1, t+1.0/2*dt)
#       l2 = v(x+1.0/2*dt*k1, y+1.0/2*dt*l1, z+1.0/2*dt*m1, t+1.0/2*dt)
#       m2 = w(x+1.0/2*dt*k1, y+1.0/2*dt*l1, z+1.0/2*dt*m1, t+1.0/2*dt) #+ derivative_K(x+1.0/2*dt*k1, y+1.0/2*dt*l1, z+1.0/2*dt*m1, t+1.0/2*dt)
#       k3 = u(x+1.0/2*dt*k2, y+1.0/2*dt*l2, z+1.0/2*dt*m2, t+1.0/2*dt)
#       l3 = v(x+1.0/2*dt*k2, y+1.0/2*dt*l2, z+1.0/2*dt*m2, t+1.0/2*dt) 
#       m3 = w(x+1.0/2*dt*k2, y+1.0/2*dt*l2, z+1.0/2*dt*m2, t+1.0/2*dt) #+ derivative_K(x+1.0/2*dt*k2, y+1.0/2*dt*l2, z+1.0/2*dt*m2, t+1.0/2*dt)
#       k4 = u(x+dt*k3,      y+dt*l3,      z+dt*m3,      t+dt)
#       l4 = v(x+dt*k3,      y+dt*l3,      z+dt*m3,      t+dt)
#       m4 = w(x+dt*k3,      y+dt*l3,      z+dt*m3,      t+dt) #+ derivative_K(x+dt*k3,      y+dt*l3,      z+dt*m3,      t+dt)
#       # dx/dt=u
#       ∂s.x = 1.0/6*(k1+2*k2+2*k3+k4)
#       # dy/dt=v
#       ∂s.y = 1.0/6*(l1+2*l2+2*l3+l4)
#       # dz/dt=0
#       ∂s.z = 1.0/6*(m1+2*m2+2*m3+m4)
#       # age=(t-t0)
#       ∂s.t = 1.0      
#    else
      # dx/dt=u
      ∂s.x = u(x, y, z, t)
      # dy/dt=v
      ∂s.y = v(x, y, z, t)
      # dz/dt=w + dk/dz
      ∂s.z = w(x, y, z, t)  #+ derivative_K(x, y, z, t) # this is quite small
      # age=(t-t0)
      ∂s.t = 1.0
   #end
end
d["f"] = f!


"""
   !g(ds,s,t,i,d)

   Dynamic model, computes as ds the function g in the equation ds=f(s,t)dt+g(s,t)dw 
   for s at current time t for particle i and possibly using data/functions from d of type userdata.
"""
function g!(∂s, s, t, i, d)
   x, y, z, age = s
#    if IS_USE_RK4
#       # currently only RK4 for background advection is used, the induced velocity by diffusion is only calculated by the most simple Euler scheme
#       ∂s.x = 0.0
#       ∂s.y = 0.0
#    # random drift velocity introduced by eddy viscosity
#       ∂s.z = sqrt(2*K(x,y,z,age))
#       ∂s.t = 0.0
#    else
      ∂s.x = 0.0
      ∂s.y = 0.0
   # random drift velocity introduced by eddy viscosity
      ∂s.z = sqrt(2*K(x,y,z,age))
      ∂s.t = 0.0
   #end
end
d["g"] = g!


"""
   !h(ds,s,t,i,d)

   Dynamic model used for implementing the M1 scheme, computes as ds the function g in the equation ds=f(s,t)dt+g(s,t)dw+h(s,t)(dw*dw-dt) 
   for s at current time t for particle i and possibly using data/functions from d of type userdata.
"""
function h!(∂s, s, t, i, d)
   x, y, z, age = s
#    if IS_USE_RK4
#       # currently only RK4 for background advection is used, the induced velocity by diffusion is only calculated by the most simple Euler scheme
#       ∂s.x = 0.0
#       ∂s.y = 0.0
#    # random drift velocity introduced by eddy viscosity
#       ∂s.z = 0.0
#       ∂s.t = 0.0
#    else
      ∂s.x = 0.0
      ∂s.y = 0.0
   # random drift velocity introduced by eddy viscosity
      ∂s.z = 0.5*K(x,y,z,age)
      ∂s.t = 0.0
   #end
end
d["h"] = h!