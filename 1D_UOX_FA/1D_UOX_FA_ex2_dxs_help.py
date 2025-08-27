#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculations to match the paper
"""

frac_t = 0.004
frac_f = 0.002
frac_s = 0.0034

#%% delta_nuSigma_f

nuSigma_f1 = 0.0057671 * 2.59068
nuSigma_f2 = 0.1062200  * 2.59068

delta_nuSigma_f1 = nuSigma_f1 * frac_f
delta_nuSigma_f2 = nuSigma_f2 * frac_f

#%% delta_Sigma_t

Sigma_t1 = 0.37790 
Sigma_t2 = 0.55064

delta_Sigma_t1 = Sigma_t1 * frac_t
delta_Sigma_t2 = Sigma_t2 * frac_t
                  
#%% delta_Sigma_s12

Sigma_s12 = 0.00086471

delta_Sigma_s12 = Sigma_s12 * frac_s

#%% delta_Sigma_a

Sigma_a1 = 0.025755
Sigma_a2 = 0.15788

#  Definition of Sigma_t (Total Cross section)>
#  Sigma_t1 = Sigma_a1 + Sigma_s1->1 + Sigma_s1->2
#  Sigma_t2 = Sigma_a2 + Sigma_s2->1 + Sigma_s2->2 = Sigma_a2 + Sigma_s2->2 

Sigma_s11 = Sigma_t1 - Sigma_a1 - Sigma_s12
Sigma_s22 = Sigma_t2 - Sigma_a2 

delta_Sigma_s11 = frac_s * Sigma_s11
delta_Sigma_s22 = frac_s * Sigma_s22

delta_Sigma_a1 = delta_Sigma_t1 - delta_Sigma_s11 - delta_Sigma_s12
delta_Sigma_a2 = delta_Sigma_t2 - delta_Sigma_s22 

#%% PRINT RESULTS
print()
print('As the pertubation is a Sine, All cuatities are imaginary:')
print(f'delta_Sigma_t1: {-delta_Sigma_t1:.8e}')
print(f'delta_Sigma_t2: {-delta_Sigma_t2:.8e}')
print(f'delta_Sigma_a1: {-delta_Sigma_a1:.8e}')
print(f'delta_Sigma_a2: {-delta_Sigma_a2:.8e}')
print(f'delta_nuSigma_f1: {-delta_nuSigma_f1:.8e}')
print(f'delta_nuSigma_f2: {-delta_nuSigma_f2:.8e}')
print(f'delta_Sigma_s12: {-delta_Sigma_s12:.8e}')

#%% SIGMA R results

print()

Sigma_r1 = Sigma_a1 + Sigma_s12
Sigma_r2 = Sigma_a2 
delta_Sigma_r1 = delta_Sigma_a1 + delta_Sigma_s12
delta_Sigma_r2 = delta_Sigma_a2 

frac_r1 = delta_Sigma_r1/Sigma_r1;
frac_r2 = delta_Sigma_r2/Sigma_r2;

print(f'frac_r1: {frac_r1:.8e}')
print(f'frac_r2: {frac_r2:.8e}')

print()
print(f'delta_Sigma_r1: {-delta_Sigma_r1:.8e}')
print(f'delta_Sigma_r2: {-delta_Sigma_r2:.8e}')
#  TD
# delta_sigma_r1 1.0651547e-04
# delta_sigma_r2 8.6717600e-04
# delta_nusigma_f1 2.9904458e-05
# delta_nusigma_f2 5.5078836e-04
# delta_sigma_s12 2.9400140e-06

# FD
 # delta_sigma_r1 (0,-0.000317247)
 # delta_sigma_r2 (0,-0.000867176)
 # delta_nusigma_f1 (0,-2.99045e-05)
 # delta_nusigma_f2 (0,-0.000550788)
 # delta_sigma_s12 (0,-2.94001e-06)


