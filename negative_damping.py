import numpy as np, sys, os
from autograd import numpy as jnp, make_jvp
import random
from scipy.integrate import odeint as ode, solve_ivp as ode2

#commandline arguments
#usage is python negative_damping.py MODEL INTERVAL_TO_ADD_PEDESTRIANS [WHETHER_TO_GENERATE_FOR_SUPPLEMENT]
#this code generates the data for plotting the multi-panel figures (4,5,6)
#the scatter plots (fig 3) can be generated by running this code for n=1 and a sinusoidally forced bridge (see comments in the integrand function)
#and then sweeping over bridge and pedestrian frequencies
#with long interval while binary searching over the final amplitude to find the max_n_walkers value for which max_amplitude exceeds a particular value.
#sample traces are included for a single pedestrian and sinusoidally forced bridge (with slightly a different output format), which were used to generate the 
which_model=int(sys.argv[1]) #1,2,3
walker_addition_interval=int(sys.argv[2]) #try 20 (seconds)
if len(sys.argv)>3:
        use_supplemental_fig_params=int(sys.argv[3])!=0
else:
        use_supplemental_fig_params=False
        

#default params
max_n_walkers=600
walker_masses=np.random.randn(max_n_walkers)*10.0+76.9
normalized_leg_lengths=np.random.randn(max_n_walkers)*0.092+1.17 
stability_margin_hof=np.random.randn(max_n_walkers)*0.002+0.0157 
stability_margin_hof*=1.0-2.0*(np.random.rand(max_n_walkers)>0.5) #pick a random foot to be down
if which_model!=3:
        step_frequency_hz=np.abs((1.04 if use_supplemental_fig_params else 0.9)+np.random.randn(max_n_walkers)*0.0)
else:
        step_frequency_hz=(1.04 if use_supplemental_fig_params else 0.9)*np.ones(max_n_walkers)
COP_lateral_offsets=stability_margin_hof*(1.0-np.tanh(0.25*np.sqrt(9.8/normalized_leg_lengths)/step_frequency_hz))
inverted_pendulum_frequency_unadjusted_radians=np.sqrt(9.8/1.17)*(1.0 if not use_supplemental_fig_params else 1.04/0.9)+0.0*np.random.randn(max_n_walkers)
stability_margin_hof*=-1
t_next_step=(np.random.rand(max_n_walkers)*0.5/step_frequency_hz)
bridge_frequency_hz=1.03
bridge_frequency_radians=bridge_frequency_hz*2*np.pi
bridge_mass=113e3
COP_offset_fixed=.63
limit_cycle_amplitude_parameter=0.47
bridge_damping_ratio=0.04
bridge_damping_coefficient=bridge_frequency_radians*bridge_damping_ratio
bridge_frequency_squared=bridge_frequency_radians*bridge_frequency_radians
limit_cycle_damping_parameter=23.25
break_on_crit=False



t=0
def jacobian(F,N,x):
        identity=np.eye(N)
        Jcols=[]
        jvp=make_jvp(F)(x)        
        for i in range(N):
          Jcols.append(jvp(identity[i])[1])
        return np.concatenate(Jcols).reshape(identity.shape)
def integrand_model1_2(state,n_walkers,COP_lateral_offsets,sagittal_velocity_perturbed,normalized_leg_lengths,walker_masses):
    y=state[2:n_walkers+2]
    ydot=state[n_walkers+2:2*n_walkers+2]
    foot_force=9.8*(COP_lateral_offsets[:n_walkers]-y[:n_walkers])/normalized_leg_lengths[:n_walkers]
    Fdot=-9.8*ydot/normalized_leg_lengths[:n_walkers]
    xdotdot=(jnp.sum(walker_masses[:n_walkers]/bridge_mass*foot_force)-bridge_damping_coefficient*state[1]-state[0]*bridge_frequency_squared)
    xdotdotdot=(jnp.sum(walker_masses[:n_walkers]/bridge_mass*Fdot)-bridge_damping_coefficient*xdotdot-state[1]*bridge_frequency_squared)
    yddot=-foot_force-xdotdot
    eta=state[3*n_walkers+2:4*n_walkers+2]
    sagittal=state[9*n_walkers+2:10*n_walkers+2]
    deriv=jnp.concatenate([jnp.array([state[1],xdotdot]),ydot,yddot, jnp.zeros(n_walkers),9.81/normalized_leg_lengths[:n_walkers]*eta-xdotdotdot*xdotdot,y*state[1], y*state[0],eta,sagittal_velocity_perturbed,sagittal*state[1],sagittal*state[0],state[:1]**2])
    return deriv
def sgn(x):#not really limit_cycle_amplitude_parameter signum function because it has sgn(0)!=0
    return 1-2*(x<0)

def integrand_model3(state,N,omegas,Ms):
        n_walkers=N
        w0=omegas[:n_walkers]
        walker_masses=Ms[:n_walkers]
        nu=w0
        y=state[2:n_walkers+2]
        ydot=state[n_walkers+2:2*n_walkers+2]
        p=COP_offset_fixed
        r=walker_masses/(bridge_mass+jnp.sum(walker_masses))
        foot_force=limit_cycle_damping_parameter*(ydot**2+nu**2*(limit_cycle_amplitude_parameter**2-(y-p*sgn(y))**2))*ydot-w0**2*(y-p*sgn(y))
        xdotdot=1/(1-jnp.sum(r))*(jnp.sum(foot_force*r)-bridge_damping_coefficient*state[1]-state[0]*bridge_frequency_squared)
        yddot=(-foot_force-xdotdot)
        dHdydot=limit_cycle_damping_parameter*(3*ydot**2+nu*nu*(limit_cycle_amplitude_parameter**2-(y-sgn(y)*p)**2))
        dHdy=limit_cycle_damping_parameter*ydot*(-nu**2*(2*y-2*p*sgn(y)))-w0**2
        Fdot=dHdy*ydot+dHdydot*yddot
        mu=state[2*n_walkers+2:3*n_walkers+2]
        eta=state[3*n_walkers+2:4*n_walkers+2]
        hu=dHdy*eta+dHdydot*mu
        xdotdotdot=1/(1-jnp.sum(r))*(jnp.sum(Fdot*r)-bridge_damping_coefficient*xdotdot-state[1]*bridge_frequency_squared)
        deriv=jnp.concatenate([jnp.array([state[1],xdotdot]),ydot,yddot, -xdotdot*xdotdotdot-hu,mu,y*state[1],y*state[0],hu,jnp.zeros(n_walkers),dHdy,dHdydot,state[:1]**2])
#        dx=jnp.concatenate([jnp.array([state[1],xdotdot]),ydot,yddot, jnp.zeros(8*n_walkers),state[:1]**2])
        return deriv

bridge_displacement=np.zeros(2)
inverted_pendulum_frequency_radians=np.sqrt(9.81/normalized_leg_lengths)
if which_model!=3:
        COM_displacement=np.zeros(1)
        COM_velocity=COP_lateral_offsets[:1]*inverted_pendulum_frequency_radians[:1]*np.tanh(inverted_pendulum_frequency_radians[:1]*0.5/step_frequency_hz[:1])
else:
        COM_velocity=np.zeros(1)
        COM_displacement=(COP_offset_fixed-limit_cycle_amplitude_parameter)*np.ones(1)
#forward_speed=(step_frequency_hz*normalized_leg_lengths/1.34/1.352)**2*step_frequency_hz
forward_speed=0.36*np.ones(max_n_walkers)
t=0
t_prev_step=t_next_step-0.5/step_frequency_hz
tprevL=(stability_margin_hof<0)*t_prev_step+(stability_margin_hof>=0)*(t_prev_step-0.5/step_frequency_hz)
COM_displacement_prev=np.zeros(max_n_walkers)
t_prev_COM_period=np.zeros(max_n_walkers)
t_prev_prev_COM_period=np.zeros(max_n_walkers)
t_prev_prev_step=np.zeros(max_n_walkers)
t_prev_prev_prev_step=np.zeros(max_n_walkers)
for n_walkers in range(2,max_n_walkers):
          t_prev_prev_step[n_walkers-1]=t
          n_periods_COM=0
          t_prev_prev_prev_step[n_walkers-1]=t-t/0.9
          step_last=np.zeros(n_walkers)
          step_first=np.ones(n_walkers)*1.0e9
          step_interval=0.5/step_frequency_hz[:n_walkers]
          time_average_d_foot_force_d_bridge_velocity=np.zeros(n_walkers)
          COM_sin_component=np.zeros(n_walkers)
          COM_cos_component=np.zeros(n_walkers)
          sagittal_sin_component=np.zeros(n_walkers)
          sagittal_cos_component=np.zeros(n_walkers)
          mu=np.zeros(n_walkers)
          eta=np.zeros(n_walkers)
          d_foot_force_d_COM_velocity=np.zeros(n_walkers)
          d_foot_force_d_COM=np.zeros(n_walkers) 
          bridge_mean_amplitude=0
          inverted_pendulum_freq_squared=(9.81/normalized_leg_lengths[n_walkers-1])
          t_next_step[n_walkers-1]+=t #account for the fact that the current most recently added walker has his next footfall way in the past, since he was initialized at the beginning of the simulation.
          t_prev_step[n_walkers-1]+=t
          if which_model!=3:
                  COM_velocity=np.concatenate([COM_velocity,COP_lateral_offsets[n_walkers-1:n_walkers]*inverted_pendulum_freq_squared*np.tanh(inverted_pendulum_freq_squared*0.5/step_frequency_hz[n_walkers-1])])
                  COM_displacement=np.concatenate([COM_displacement,np.zeros(1)])
          else:
                  random_side=random.choice([-1,1])
                  COM_velocity=np.concatenate([COM_velocity,np.zeros(1)])
                  COM_displacement=np.concatenate([COM_displacement,random_side*(COP_offset_fixed-limit_cycle_amplitude_parameter)*np.ones(1)])
          d_foot_force_d_bridge_velocity=np.zeros(n_walkers)
          jump_term_1=np.zeros(n_walkers)
          jump_term_2=np.zeros(n_walkers)
          jump_term_3=np.zeros(n_walkers)
          jump_term_4=np.zeros(n_walkers)
          sagittal_velocity_perturbed=np.zeros(n_walkers)
          num_steps=np.zeros(n_walkers)
          t_prev_addition=t
          has_zero_cross=False
          n_footfalls=0
          tfirstp=t_prev_addition
          sagittal=np.zeros(n_walkers)
          saggital_prev=np.zeros(n_walkers)
          t_first_bridge_zero=walker_addition_interval+t_prev_addition
          t_final_bridge_zero=t_prev_addition
          bridge_max_amplitude=0
          num_zero_cross=0
          order_param_cos_component_step=0
          order_param_sin_component_step=0
          order_param_cos_component_COM=0
          order_param_sin_component_COM=0
          while t<t_prev_addition+walker_addition_interval:
            if which_model!=3:
                integrand_autonomous=lambda x: integrand_model1_2(x,n_walkers,COP_lateral_offsets[:n_walkers],sagittal_velocity_perturbed[:n_walkers],normalized_leg_lengths[:n_walkers],walker_masses[:n_walkers])
                t_next_footfall_all_walkers=max(t,np.min(t_next_step[:n_walkers]))
            else:
                integrand_autonomous=lambda x: integrand_model3(x,n_walkers,inverted_pendulum_frequency_unadjusted_radians[:n_walkers], walker_masses[:n_walkers])        
                t_next_footfall_all_walkers=min(t+0.002,t_prev_addition+walker_addition_interval)
            #Fn_=jit(Fn__)
            integrand=lambda t,x: integrand_autonomous(x)
            ts=np.linspace(0,t_next_footfall_all_walkers-t,15 if which_model==3 else 10)
            if which_model==3:
                state=np.concatenate([bridge_displacement[:2],COM_displacement[:n_walkers],COM_velocity[:n_walkers], mu,eta,  COM_sin_component,COM_cos_component,d_foot_force_d_bridge_velocity,np.zeros(n_walkers),d_foot_force_d_COM,d_foot_force_d_COM_velocity, bridge_mean_amplitude*np.ones(1)])
            else:
                state=np.concatenate([bridge_displacement[:2],COM_displacement[:n_walkers],COM_velocity[:n_walkers], mu,eta, COM_sin_component,COM_cos_component,d_foot_force_d_bridge_velocity,sagittal, sagittal_sin_component,sagittal_cos_component,bridge_mean_amplitude*np.ones(1)])
            state_prev_sample=state.copy()
            J_=lambda t,x: jacobian(integrand_autonomous,len(state),x)
            trajectory_current_sample=ode2(integrand, [0, ts[-1]],state, t_eval=ts, method='LSODA', jac=J_, rtol=1e-6, atol=1e-8).y.transpose()
            state=trajectory_current_sample[-1]
            for i,trajectory_current_subsample in enumerate(trajectory_current_sample[:]):
                if i<2:
                        continue
                tsi=ts[i]
            #    if which_model!=3:
            #            sys.stderr.write('%i %f %f %f '%(n_walkers,ts[i]+t,trajectory_current_subsample[0],trajectory_current_subsample[1])+('%f '*n_walkers)%tuple((COP_lateral_offsets[:n_walkers]-trajectory_current_subsample[2:n_walkers+2]).tolist())+'\n_walkers')
            #    else:
            #            sys.stderr.write('%i %f %f %f '%(n_walkers,ts[i]+t,trajectory_current_subsample[0],trajectory_current_subsample[1])+('%f '*n_walkers)%tuple((trajectory_current_subsample[2:n_walkers+2]).tolist())+'\n_walkers')
                if trajectory_current_sample[i-2,0]<trajectory_current_sample[i-1,0] and trajectory_current_sample[i-1,0]>trajectory_current_sample[i,0]: 
                    if has_zero_cross:
                            COM_sin_component_full_period=trajectory_current_sample[i-1,4*n_walkers+2:5*n_walkers+2]
                            COM_cos_component_full_period=trajectory_current_sample[i-1,5*n_walkers+2:6*n_walkers+2]
                            if which_model==3:
                                d_foot_force_d_COM_full_period=np.array(trajectory_current_sample[i-1,-2*n_walkers-1:-n_walkers-1])
                                d_foot_force_d_COM_velocity_full_period=np.array(trajectory_current_sample[i-1,-n_walkers-1:-1])
                            else:
                                sagittal_sin_component_full_period=trajectory_current_sample[i-1,-2*n_walkers-1:-n_walkers-1]
                                sagittal_cos_component_full_period=trajectory_current_sample[i-1,-n_walkers-1:-1]
                            mean_amplitude=trajectory_current_subsample[-1]
                    num_zero_cross+=1
                    has_zero_cross=True
                    t_final_bridge_zero=max(tsi+t,t_final_bridge_zero)
                    t_first_bridge_zero=min(tsi+t,t_first_bridge_zero)
            if which_model!=3:
                    sys.stderr.write('%i %f %f %f '%(n_walkers,t,state[0],state[1])+('%f '*n_walkers + '0.0 '*(max_n_walkers-n_walkers))%tuple((COP_lateral_offsets[:n_walkers]-state[2:n_walkers+2]).tolist())+'\n')
            else:
                    sys.stderr.write('%i %f %f %f '%(n_walkers,t,state[0],state[1])+('%f '*n_walkers+ ' 0.0 '*(max_n_walkers-n_walkers))%tuple(state[2:n_walkers+2].tolist())+'\n')
            sys.stderr.flush()
            sagittal=state[7*n_walkers+2:8*n_walkers+2]
            bridge_displacement=state[:2]
            
            if has_zero_cross:      
              COM_sin_component=state[4*n_walkers+2:5*n_walkers+2]
              COM_cos_component=state[5*n_walkers+2:6*n_walkers+2]
              d_foot_force_d_bridge_velocity=state[6*n_walkers+2:7*n_walkers+2]
              if which_model!=3:
                      sagittal_sin_component=state[-2*n_walkers-1:-n_walkers-1]
                      sagittal_cos_component=state[-n_walkers-1:-1]
              else:
                     d_foot_force_d_COM=state[-2*n_walkers-1:-n_walkers-1]
                     d_foot_force_d_COM_velocity=state[-n_walkers-1:-1]
              mu=state[2*n_walkers+2:3*n_walkers+2]
              eta=state[3*n_walkers+2:4*n_walkers+2]
              bridge_mean_amplitude=state[-1]
            k=integrand_autonomous(state)
            t=t_next_footfall_all_walkers
            bridge_displacement=state[:2]
            COM_displacement=state[2:n_walkers+2]
            COM_velocity=state[n_walkers+2:2*n_walkers+2]
            if which_model!=3:
                    indices_foot_down=np.where(t_next_step[:n_walkers]<=t+1.0e-10)[0]
            else:
                    indices_foot_down=np.where(sgn(state_prev_sample[2:n_walkers+2])!=sgn(state[2:n_walkers+2]))[0]
            time_average_d_foot_force_d_bridge_velocity[indices_foot_down]+=d_foot_force_d_bridge_velocity[indices_foot_down]
            t_prev_step[indices_foot_down]=t
            if which_model!=3:
                    for i in indices_foot_down:
                       if stability_margin_hof[i]<0:
                         tprevL[i]=t_prev_step[i]
                    u2=COM_displacement[indices_foot_down]+COM_velocity[indices_foot_down]*np.sqrt(normalized_leg_lengths[indices_foot_down]/9.8)+stability_margin_hof[indices_foot_down]
                    width=stability_margin_hof/(1-np.tanh(np.sqrt(9.81/normalized_leg_lengths)/4.0/step_frequency_hz))
                    adapt = np.maximum(-0.5,(width[indices_foot_down]**2-(COM_displacement[indices_foot_down]-u2)**2)/(4*forward_speed[indices_foot_down]**2)) if which_model==2 else 0
                    t_next_step[indices_foot_down]=t+0.5/step_frequency_hz[indices_foot_down]*(1+adapt)
            step_last[indices_foot_down]=np.maximum(step_last[indices_foot_down],t)
            step_first[indices_foot_down]=np.minimum(step_first[indices_foot_down],t)
            num_steps[indices_foot_down]+=1
            if which_model!=3:
                    saggital_prev[indices_foot_down]=sagittal[indices_foot_down]
            COM_displacement_prev[indices_foot_down]=COM_displacement[indices_foot_down]
            indices_COM_period=np.where(np.signbit(state_prev_sample[2:n_walkers+2])!=np.signbit(state[2:n_walkers+2]))[0]
            t_prev_prev_COM_period[indices_COM_period]=t_prev_COM_period[indices_COM_period]
            t_prev_COM_period[indices_COM_period]=t
            step_interval_prev=np.zeros(len(step_interval))
            step_interval_prev[:]=step_interval[:]
            if which_model==3:
                    com_offset_current=state[2:n_walkers+2]
                    com_offset_prev=state_prev_sample[2:n_walkers+2]
                    force_change_over_step_discontinuity=2*inverted_pendulum_frequency_unadjusted_radians[indices_foot_down]**2*np.array(COP_offset_fixed*sgn(com_offset_current[indices_foot_down]-com_offset_prev[indices_foot_down]))
                    sagittal_velocity_perturbed[indices_foot_down]=0.0
            else:
                    force_change_over_step_discontinuity=(u2-COP_lateral_offsets[indices_foot_down])*9.8/normalized_leg_lengths[indices_foot_down]
                    step_interval[indices_foot_down]=t_next_step[indices_foot_down]-t_prev_step[indices_foot_down]
                    sagittal_velocity_perturbed[indices_foot_down]=1-step_interval_prev[indices_foot_down]/step_interval[indices_foot_down]
                    COP_lateral_offsets[indices_foot_down]=u2
            jump_term_1[indices_foot_down]+=(force_change_over_step_discontinuity)*k[1]
            jump_term_2[indices_foot_down]+=(force_change_over_step_discontinuity)*np.array(k[2:n_walkers+2][indices_foot_down])
            jump_term_3[indices_foot_down]+=(force_change_over_step_discontinuity)*sagittal_velocity_perturbed[indices_foot_down]
            jump_term_4[indices_foot_down]+=(force_change_over_step_discontinuity)*np.array(k[2+n_walkers:2+2*n_walkers][indices_foot_down])
            t_prev_prev_prev_step[indices_foot_down]=t_prev_prev_step[indices_foot_down]
            t_prev_prev_step[indices_foot_down]=t_prev_step[indices_foot_down]
            t_prev_step[indices_foot_down]=t

            stability_margin_hof[indices_foot_down]*=-1
            com_offset_current=state[2:n_walkers+2]
            com_offset_prev=state_prev_sample[2:n_walkers+2]
            if which_model!=3:
                    phZ=(t-t_prev_step[:n_walkers])/(t_next_step[:n_walkers]-t_prev_step[:n_walkers])
                    order_param_cos_component_step+=np.sum(np.cos(2*np.pi*phZ))
                    order_param_sin_component_step+=np.sum(np.sin(2*np.pi*phZ))
            phase_COM=(t-t_prev_COM_period[:n_walkers])*2*np.pi/(t_prev_COM_period[:n_walkers]-t_prev_prev_COM_period[:n_walkers])
            phase_COM=phase_COM[np.where(t_prev_COM_period>t_prev_prev_COM_period)[0]]
            order_param_cos_component_COM+=np.sum(np.cos(2*np.pi*phase_COM))
            order_param_sin_component_COM+=np.sum(np.sin(2*np.pi*phase_COM))
            n_footfalls+=1
            n_periods_COM+=len(phase_COM)
            mu[indices_foot_down]=0
            eta[indices_foot_down]=0
            d_foot_force_d_bridge_velocity[indices_foot_down]=0
            bridge_max_amplitude=max(bridge_max_amplitude,state[0])
          time_from_first_to_last_footfall=(step_last[:n_walkers]-step_first[:n_walkers])
          mean_step_frequency=np.mean(0.5*num_steps[:n_walkers]/time_from_first_to_last_footfall)
          stdev_step_frequency=np.std(0.5*num_steps[:n_walkers]/time_from_first_to_last_footfall)
          total_integration_time=t-t_prev_addition
          time_from_first_to_last_zero_cross_bridge=t_final_bridge_zero-t_first_bridge_zero
          COM_sin_component_full_period/=(time_from_first_to_last_zero_cross_bridge*np.sqrt(mean_amplitude/time_from_first_to_last_zero_cross_bridge))
          COM_cos_component_full_period/=(time_from_first_to_last_zero_cross_bridge*np.sqrt(mean_amplitude/time_from_first_to_last_zero_cross_bridge))
          if which_model!=3:
                  sagittal_sin_component_full_period/=(time_from_first_to_last_zero_cross_bridge*np.sqrt(mean_amplitude/time_from_first_to_last_zero_cross_bridge))
                  sagittal_cos_component_full_period/=(time_from_first_to_last_zero_cross_bridge*np.sqrt(mean_amplitude/time_from_first_to_last_zero_cross_bridge))
                  d_foot_force_d_COM_full_period=-inverted_pendulum_frequency_radians[:n_walkers]**2*time_from_first_to_last_footfall
                  d_foot_force_d_COM_velocity_full_period=0
                  time_average_d_foot_force_d_bridge_velocity*=9.81/normalized_leg_lengths[:n_walkers]
          sigma_1=np.sum((time_average_d_foot_force_d_bridge_velocity+jump_term_1)/time_from_first_to_last_footfall)
          sigma_2=np.sum(-(jump_term_2+d_foot_force_d_COM_full_period)*COM_sin_component_full_period/bridge_frequency_radians/time_from_first_to_last_footfall+(jump_term_4+d_foot_force_d_COM_velocity_full_period)*COM_cos_component_full_period/time_from_first_to_last_footfall)
          sigma_3=np.sum(sagittal_sin_component_full_period*jump_term_3/time_from_first_to_last_footfall)/bridge_frequency_radians if which_model!=3 else 0
          order_footsteps=np.sqrt(order_param_cos_component_step*order_param_cos_component_step+order_param_sin_component_step*order_param_sin_component_step)/(n_walkers*n_footfalls)
          order_COM=np.sqrt(order_param_cos_component_COM*order_param_cos_component_COM+order_param_sin_component_COM*order_param_sin_component_COM)/(n_periods_COM)
          print(n_walkers,
				bridge_damping_coefficient-n_walkers*np.mean(walker_masses[:n_walkers])/bridge_mass*(sigma_1+sigma_2+sigma_3), 
				bridge_max_amplitude,
				order_footsteps,  
				order_COM, 
				mean_step_frequency, 
				stdev_step_frequency, 
				np.mean(forward_speed), 
				np.sqrt(mean_amplitude/time_from_first_to_last_zero_cross_bridge), 
				num_zero_cross/(time_from_first_to_last_zero_cross_bridge), 
				sigma_1, sigma_2, sigma_3
				)
#          if bridge_damping_coefficient-n_walkers*np.mean(walker_masses[:n_walkers])/bridge_mass*(sigma_1+sigma_2+sigma_3)<0 and break_on_crit:
          if bridge_max_amplitude>0.015 and break_on_crit:
                sys.exit(n_walkers)
          sys.stdout.flush()
sys.exit(max_n_walkers)
                  
