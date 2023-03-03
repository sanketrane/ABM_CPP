#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:13:25 2020
@author: maria
"""

import numpy as np
from joblib import Parallel, delayed
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.stats import multivariate_normal
import pandas as pd
import multiprocessing



def gamma_age(age, params):
    t0 = 5
    N0, delta, rho, r = params
    # Flat, normalised age-distribution of initial cells
    if age>=0 and age<=t0:
        return N0/t0
    else:
        return 0.

## turnover rate
def delta_age(age, params):
    N0, delta, rho, r = params
    if age>=0:
        return delta * np.exp(-r * age)
    else:
        return 0.

## division rate
def rho_age(age, params):
    N0, delta, rho, r = params
    q = 0.
    if age>=0:
        return rho * np.exp(-q * age)
    else:
        return 0.

## net growth rate
def netgrowth(age, params):
    return rho_age(age, params) - delta_age(age, params)

## net process
def netprocess(age, params):
    return rho_age(age, params) + delta_age(age, params)

# Timecourse of thymic CD4 and CD8 SP populations -- change with time
def SP_numbers(Time, Population):
    t0 = 5
    dpt0 = Time-t0 #days post t0

    ## parameters taken from spline fits
    if Population == "CD4":
        theta_0  = 4.3e5
        theta_f = 1.8e3
        n = 2.1
        X = 30
        q = 3.7
    if Population == "CD8":
        theta_0  = 9e4
        theta_f = 68
        n = 3
        X = 25
        q = 4.25

    value = theta_0 + (theta_f * dpt0**n)*(1 - ((dpt0**q)/((X**q) + (dpt0**q))))

    if Time <t0:
        return 0.
    else:
        return value

# Influx from the thymus = theta(t) = psi*SP.numbers where psi = constant
# We calculate psi from the boundary condition; zero-age cells at time t0= thymic output at t0
# We assume the age-distribution of cells at t=t0 is flat and equal to gamma(a) = 1/t0 (0<=a<=t0)
# so psi*SP.numbers(t0)  = N0 * gamma(0)

# Total influx into the naive T cell compartment from thymus (cells/day)
# can take multiple time values as argument
def theta_fun(Time, params, Population):
    t0 = 5
    psi = gamma_age(0, params)/SP_numbers(t0, Population)
    value = psi * SP_numbers(Time, Population)
    return value

# Ki67 proportion in thymic output -- changes with time -- spline (estimated from data)
def eps_func(Time, Population):
    ## spline fitted separately to the proportions of ki67+ cells within thymic SP4 cells
    ## parameters estimated from spline fit to the timecourse of ki67 proportions within the source compartment -- SP
    if Population == "CD4":
        eps_0  = 0.13732103
        eps_f = 0.03337899
        A = 2.92554110
    if Population == "CD8":
        eps_0  = 0.24510453
        eps_f = 0.01559996
        A = 14.83715328
    eps = np.exp(-eps_f * (Time + A)) + eps_0
    return eps

# k_dist_thymus is the distribution of ki67 expression in thymic output -- changes with time 
# function can return multiple values of k for one value of Time
# We found a smooth function that reproduces the timecourse of eps_func above

def k_dist_thymus(ki, Time, Population):
    n = 40
    k_bar = 1./np.exp(1)

    if ki >0 and ki<=1:
        value = (eps_func(Time, Population)/(1. - k_bar)) + (
            (((1 - eps_func(Time, Population))/(k_bar)) - (eps_func(Time, Population)/(1 - k_bar)))/(1 + (ki/k_bar)**n))
    else:
        value = 0.
    return value

# Finally, Ki67 dist of cells in periphery at t0
# option1: we make Ki67 distribution at t=t0 the same for all cell ages (0-> t0),
# and equal to the Ki67 distribution of cells leaving the thymus at t0
#k_dist_init=function(k){k_dist_thymus(k, Time=t0)}

# Option 2 - make it skewed in the other direction
def k_dist_init(k):
    if k>0 and k<=1:
        return np.exp(-k/3.5)/0.8698295
    else:
        return 0.

# The choice of this distribution matters a lot for how the initial Ki67+ fraction looks
# If this is a problem, we could try linking it to age - 
# i.e. cells born between (0,t0) entered Ki67 hi and now have Ki67=exp(-beta a)

#####################################################################################################################################
### solution of ASM PDE without K as a variable 

#--- Initial cohort ---
### age distribution of the total initial cohort at time 't' that can potentially enter division with the rate rho
def integral_fun(fun, a,b, param):
    '''
    The function integrates function fun from a to b
    using module quad from scipy.integrate library. 
    '''
    res, err = quad(fun, a, b, args=(param))
    return res

def Pooled_init_age(age, time, params):
    # Eq (17) U_i(a,t)
    # return single value for age, time
    t0 = 5
    g_func = gamma_age(t0, params)
    value = g_func * np.exp(integral_fun(netgrowth, age-(time-t0),age,params))
    return value

## total counts of initial cohort integrated across all possible cell ages
def Pooled_init_time(time, params):
    # int_t-t0^t Eq (17)
    # return vector
    t0 =5
    vec = [integral_fun(Pooled_init_age, t - t0, t, (t,params)) for t in time]
    return np.array(vec)


#--- Theta cohort ---
### age distribution of the total theta cohort at time 't' that can potentially enter division with the rate rho
def Pooled_theta_age(age, time, params, Population):
    # Eq (18) U_theta(a,t)
    # return single value for age, time
    t0 = 5
    if age <= time -t0:
        value = theta_fun(time - age, params, Population) * np.exp(integral_fun(netgrowth, 0, age, params))
    else:
        value = 0.
    return value

## total counts of theta cohort integrated across all possible cell ages
def Pooled_theta_time(time, params, Population):
    t0 =5
    vec = [integral_fun(Pooled_theta_age,0, t - t0, (t,params, Population)) for t in time]
    return np.array(vec)

def Asm_total_age(age, time, params, Population):
    # function of age and time
    # return singe value for age, time
    t0 = 5
    if age <= time - t0:
        return Pooled_theta_age(age, time, params, Population)
    else:
        return Pooled_init_age(age, time, params)

# used to fit data
def Asm_total_Time(time, params, Population):
    newtime = [0]+list(time)
    vec = [integral_fun(Asm_total_age, newtime[i], newtime[i+1], (newtime[i+1], params, Population)) for i in range(0, len(newtime)-1)]
    vec =  np.array(vec)
    return np.cumsum(vec)

#####################################################################################################################################
## tracking preexisting cells that were present at t0 -- initial cohort

#####
## UnDivided cells within initial cohort
###

### age >= Time -t0; ki67 <= exp(-beta * (Time -t0))
def u_init_nd_k_a(age, ki, Time, params, fixed_params):
    beta = fixed_params
    t0 = 5
    if age >= Time - t0 and ki<= np.exp(-beta * (Time - t0)):
        value = gamma_age(age - Time + t0, params)*k_dist_init(ki * np.exp(beta * (Time - t0)))*np.exp(beta * (Time - t0)) *np.exp(- integral_fun(netprocess, age - Time + t0, age, params))
    else:
        value = 0.0
    return value

#####
## Divided cells within initial cohort
###

### age >= Time -t0; ki67 >= exp(-beta * (Time -t0))
def u_init_div_k_a(age, ki, Time, params, fixed_params):
    beta = fixed_params
    t0 = 5
    tau = -np.log(ki)/beta    ### time since cell division

    if age >= Time - t0 and ki >= np.exp(-beta * (Time - t0)):
        value = 2 * rho_age(age - tau, params) * Pooled_init_age(age  - tau, Time - tau, params) * (1/(beta * ki)) *np.exp(- integral_fun(netprocess,age -tau, age, params))
    else:
        value = 0
    return value

# total density of the initial cohort u_init(a, k, t)
def u_init_k_a(age, ki, Time, params, fixed_params):
    beta = fixed_params
    t0 = 5
    return u_init_nd_k_a(age, ki, Time, params, beta) + u_init_div_k_a(age, ki, Time, params, beta)


#####################################################################################################################################
#####################################################################################################################################
### Cells coming from source -- cells of age 0 -- theta cohort
### age >= Time -t0; ki67 <= exp(-beta * (Time -t0))

#########
## UnDivided cells Theta cohort
#########

### age <= Time -t0; ki67 <= exp(-beta * age)
def u_theta_nd_k_a(age, ki, Time, params, fixed_params):
    beta, Population = fixed_params
    t0 = 5
    if age < Time - t0 and ki <= np.exp(-beta * age):
        return theta_fun(Time - age, params, Population) * k_dist_thymus(ki * np.exp(beta * age), Time - age, Population) *np.exp(beta*(age) - integral_fun(netprocess, 0, age, params))
    else:
        return 0

#########
## Divided cells Theta cohort
#########
# divided

### age <= Time -t0; ki67 >= exp(-beta * (Time -t0))
def u_theta_div_k_a(age, ki, Time, params, fixed_params):
    beta, Population = fixed_params
    t0 = 5
    tau = -(1/beta)*np.log(ki) # Cell with ki67=k divided this time ago

    if age < Time - t0 and ki >= np.exp(-beta * age):
        return 2 * rho_age(age - tau, params)/(beta * ki) * Pooled_theta_age(age - tau, Time - tau, params, Population) *np.exp(-integral_fun(netprocess, age-tau, age, params))
    else:
        return 0

# total density of the theta cohort u_theta(a, k, t)
def u_theta_k_a(age, ki, Time, params, fixed_params):
    beta, Population = fixed_params
    t0 = 5
    return u_theta_nd_k_a(age, ki, Time, params, fixed_params) + u_theta_div_k_a(age, ki, Time, params, fixed_params)

#####################################################################################################################################
#####################################################################################################################################
# total density of u(k, a, t)
def u_total_a_k(ki, age, Time, params, fixed_params):
    beta, Population = fixed_params
    t0 = 5

    if age < Time - t0:
        return u_theta_k_a(age, ki, Time, params, fixed_params)
    else:
        return u_init_k_a(age, ki, Time, params, beta)

## age dist u(a,t) of ki67+ cells
def Kpos_total_age(age, time, params, fixed_params):
    # integrate on ki, return single value
    t0 = 5
    k_bar = 1/np.exp(1)
    return integral_fun(u_total_a_k, k_bar, 1, (age, time, params, fixed_params))

## counts u(t) of ki67+ cells
def Kpos_total_Time(time, params, fixed_params):
    # integrate on age, return vector
    newtime = [0]+list(time)
    vec = [integral_fun(Kpos_total_age, newtime[i], newtime[i+1], (newtime[i+1], params, fixed_params)) for i in range(0,len(newtime)-1)]
    vec = np.array(vec)
    return np.cumsum(vec)

def ASM_fun(time, params, fixed_params):
    beta, Population = fixed_params
    pred_counts = Asm_total_Time(time, params, Population)
    pred_ki_counts = Kpos_total_Time(time, params, fixed_params)
    pred_ki_prop = pred_ki_counts/pred_counts
    return pred_counts, pred_ki_prop, time


# Upload the data
DD = pd.read_csv('abc_sanket/CD4_ln.csv')

# Sample 1 sets of parameters from prior
# change here all prior for theta = [psi, rho, delta] if you like
# remember to also update lines after while (7-9) and (20-22)
def sample_one_from_pi(Population):
    np.random.seed()
    # sample theta
    logN0 = np.random.normal(12,2)
    rho = np.random.normal(0.005,0.01)
    delta = np.random.normal(0.05,0.1)
    r = np.random.normal(0.01, 0.01)
    theta = [logN0,delta, rho, r]
    params = [np.exp(logN0), delta, rho, r]
    res1, res2, tres = ASM_fun(np.arange(5,300,1), params, [1/3.5, Population])
    return theta, [res1, res2, tres]

# Define the distance between simulation and data
def dist(sim, DD):
    res1, res2, tres = sim
    CC = DD['counts'].values
    KI = DD['ki67'].values
    tdata = DD['time'].values
    age_index = []
    for j in range(0, len(tdata)):
        age_index.append(np.where(tres == tdata[j])[0][0])
    resi_cc = np.log(CC) - np.log(res1[age_index])
    resi_ki67 =  np.arcsin(KI) - np.arcsin(res2[age_index])
    res = np.sum(resi_cc**2) + np.sum(resi_ki67**2)
    return res

# Sample 1 theta from prior defined above and accept if |simulation - data| < pre-defined distance
def rej_samp_popt0(i, eps, Population):
    '''
    Sample theta from pi(theta) for population t=0, particle n=i
    '''
    theta, sim = sample_one_from_pi(Population)
    dis = dist(sim, DD)
    while dis>eps:
        theta, sim = sample_one_from_pi(Population)
        dis = dist(sim, DD)
    return theta

# Sample parallel N theta from prior defined above and accept if |simulation - data| < pre-defined distance
def rej_samp_Npart_popt0(N, eps, Population, nCores):
    '''
    Sample N particles from pi(theta)
    '''
    theta = Parallel(n_jobs=nCores)(delayed(rej_samp_popt0)(i, eps, Population) for i in range(0,N)) 
    return theta

# Sample 1 theta from previous generation of posterior
def sample_one_from_theta(theta_n_particles, weights):
    j = np.random.choice(np.arange(0, len(weights),1),p=weights)
    theta = theta_n_particles[j]
    return theta

# Check if perturb theta is still within posterior bounds
# used also to calculate new weights
def joint_pi(theta,abound, bbound):
    pi = np.prod(1/(np.array(bbound)-np.array(abound)))
    for i in range(0, len(theta)):
        if theta[i]<abound[i] or theta[i]>bbound[i]:
            pi = pi*0
    return pi

# used to calculate new weights with gaussian kernel
def kernel(theta_t, theta_tm1):
    cov_theta = 2*np.cov(theta_tm1.T)
    pi_kernel = []
    for i in range(0, len(theta_tm1)):
        pi_kernel.append(multivariate_normal.pdf(theta_t, mean=theta_tm1[i], cov=cov_theta, allow_singular=True))
    return pi_kernel

def move_one_particle(i,theta_prev_pop, weights, eps, Population):
    '''
    Move particle i for population t from theta(t-1) = theta_prev_pop
    '''
    np.random.seed()
    if np.mod(i,10) == 0:
        print('Particle: ', i)
    # move particle with kernel
    # UPDATE BOUNDS FOR PARAMETERS HERE IF THEY ARE DIFFERENT
    abound = [-20, 0, 0, 0]
    bbound = [20, 1, 1,1]
    cov_theta = 2*np.cov(theta_prev_pop.T)

    dis = eps+1
    pi_newtheta = 0
    NoOfSamp = 0
    while dis>eps or pi_newtheta == 0:
        # Sample theta from previous theta with weights
        theta = sample_one_from_theta(theta_prev_pop, weights)
        # move particle with kernel
        newtheta = multivariate_normal.rvs(theta, cov_theta)
        # Check if you are in pi(theta)
        pi_newtheta = joint_pi(newtheta, abound, bbound)

        # Check the distance
        logN0, delta, rho, r = newtheta
        params = [np.exp(logN0), delta, rho, r]
        res1, res2, tres = ASM_fun(np.arange(5,300,1), params, [1/3.5, Population])
        sim =  [res1, res2, tres]
        dis = dist(sim, DD)
        NoOfSamp = NoOfSamp +1
    newweights = joint_pi(newtheta, abound, bbound)/np.sum(np.array(weights)*np.array(kernel(newtheta,theta_prev_pop)))
    return newtheta, newweights, dis, NoOfSamp

# Move N particles
def move_N_particles(N, theta_prev, weights, eps,Population, nCores):
    res = Parallel(n_jobs=nCores)(delayed(move_one_particle)(i,theta_prev, weights, eps, Population) for i in range(0,N))
    theta, weights, dis, NoOfSamp = zip(*res)
    return theta, weights, dis, NoOfSamp

def abcsmc(eps_max, eps_min, N, Population, nCores):
    '''
    eps = eps[0], ..., eps [T]  - T no of populations
    N - no of particles
    '''
    t=0
    eps = eps_max
    while eps > eps_min:
        print('Population', t, ', epsilon: ',eps)
        if t ==0:
            # Sample N partocles from pi(theta)
            thetaSS = rej_samp_Npart_popt0(N, eps, Population, nCores)
            # Assign weights 1/N
            weights = [1/N]*N
            # Theta(t=0) = [thetaSS]
            thetaS = np.array(thetaSS)
            dis = 0
            NoOfSamp = 0
        else:
            thetaSS, weightsSS, dis, NoOfSamp = move_N_particles(N, thetaS, weights, eps,Population, nCores)
            weights = np.array(weightsSS)/np.sum(np.array(weightsSS))
            thetaS = np.array(thetaSS)
            res = pd.DataFrame()
            res['N0'] = np.exp(thetaS.T[0])
            res['delta'] = thetaS.T[1]
            res['rho'] = thetaS.T[2]
            res['r'] = thetaS.T[3]
            res['weights'] = weights
            res['distance'] = dis
            res['NoOfSamp'] = NoOfSamp
            res.to_csv('abc_sanket/post_asm_eps_%s.csv'%(np.round( eps,7)))
            eps = np.median(dis)
        t =t+1
    res.to_csv('abc_sanket/post_asm_final.csv')
    return

# First always test for large max and for only 100 sample to see where you are
# you can print distance so you can see what is accepted
maxDistance = 50
minDistance = 6.35
No_of_theta_in_post = 1000
Pop = "CD4"
abcsmc(maxDistance, minDistance,No_of_theta_in_post, Pop, nCores = multiprocessing.cpu_count()-2 )

