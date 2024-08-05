# A script for generating the average-case and worst-case heuristic estimates for BGV noise growth in specified circuits
# Circuits considered are the [CLP20] circuit and the "bgv deep" circuit
# For comparison with the noise observed in the implementations of these circuits in HElib 2.2.1 and in SEAL 4.0

###########
# Imports #
###########

from math import sqrt
from math import floor
from math import log
from scipy.special import erfinv


############################################################################
# Worst-case bounds after operations as presented in [Iliashenko19, CLP20] #
############################################################################

def bound_fresh(n, t):
    sigma = 3.19
    inside_sqrt = (4./3)*n + 1
    inside_sqrt *= n*sigma**2
    inside_sqrt += n/(12.)
    output = sqrt(inside_sqrt)
    output *= 6*t
    return output

def bound_add(input_bound_1, input_bound_2):
    return input_bound_1 + input_bound_2

def bound_mult(input_bound_1, input_bound_2):
    return input_bound_1 * input_bound_2

def bound_mod_switch(n, t, q, p, input_bound):
    inside_sqrt = 3*n + 2*n*n
    output_bound = sqrt(inside_sqrt)
    output_bound *= t
    output_bound += ((p/q)*input_bound)
    return output_bound


#####################################################################
# Average-case variances after operations, as presented in Figure 5 #
#####################################################################

def variance_fresh(n, t):
    sigma = 3.19
    output_variance = ((4.0/3) * n + 1) * t * t * sigma * sigma
    return output_variance

def variance_add(input_variance_1, input_variance_2):
    return input_variance_1 + input_variance_2

def variance_mult(input_variance_1, input_variance_2, n, t):
    term1 = n * input_variance_1 * input_variance_2
    term2 = input_variance_1 * n * (1./12) * (t * t -1) # component m_i uniform mod t, giving |m| \approx n*(t^2-1)/12
    term3 = input_variance_2 * n * (1./12) * (t * t -1)
    return term1 + term2 + term3

def variance_mod_switch(n, t, q, p, input_variance):
    gamma_squared_input_variance = (p/q) * (p/q) * input_variance
    output_variance = (1.0/12) * ((2.0/3) * n + 1) * (t * t - 1)
    output_variance += gamma_squared_input_variance
    return output_variance


###############################################################
# Calculate noise budget remaining after evaluating a circuit #
###############################################################

# Given variance of the noise in the output ciphertext, compute a bound on the noise, in the manner of [CCH+21]   
def alpha_bound_from_variance(variance, n):
    alpha = 0.001
    bound = sqrt(2*variance) * erfinv(pow(1-alpha,1.0/n))
    return bound

# Given bound on the noise in the output ciphertext, calculate the noise budget remaining
def get_noise_budget(bound, q):
    output = floor(log(q,2)-log(bound,2)) - 1
    return output


###########################################################################
# Estimates of worst-case and average-case heuristics for [CLP20] circuit #
###########################################################################

# Variance of the output ciphertext after evaluating the [CLP20] circuit
def variance_after_clp20_computation(n, t, q, p):
    fresh = variance_fresh(n, t)
    add = variance_add(fresh, fresh)
    mult = variance_mult(add, fresh, n, t)

    if (n > 2048):
        mod_switch = variance_mod_switch(n, t, q, p, mult)
        return fresh, add, mult, mod_switch
    else:
        return fresh, add, mult

# Top-level function for noise budget predicted for average-case approach
def average_case_clp20(n, t, q, p):
    if (n > 2048):
        fresh_ctext_variance, add_ctext_variance, mult_ctext_variance, mod_switch_variance = variance_after_clp20_computation(n, t, q, p)
    else:
        fresh_ctext_variance, add_ctext_variance, mult_ctext_variance = variance_after_clp20_computation(n, t, q, p)

    fresh_ctext_bound = alpha_bound_from_variance(fresh_ctext_variance, n)
    add_ctext_bound = alpha_bound_from_variance(add_ctext_variance, n)
    mult_ctext_bound = alpha_bound_from_variance(mult_ctext_variance, n)
    fresh_noise_budget = get_noise_budget(fresh_ctext_bound, q)
    add_noise_budget = get_noise_budget(add_ctext_bound, q)
    mult_noise_budget = get_noise_budget(mult_ctext_bound, q)

    if (n > 2048):
        mod_switch_bound = alpha_bound_from_variance(mod_switch_variance, n)
        mod_switch_noise_budget = get_noise_budget(mod_switch_bound, p)
        return fresh_noise_budget, add_noise_budget, mult_noise_budget, mod_switch_noise_budget
    else:
        return fresh_noise_budget, add_noise_budget, mult_noise_budget

# Top-level function for noise budget predicted for worst-case approach
def worst_case_clp20(n, t, q, p):
    fresh = bound_fresh(n, t)
    add = bound_add(fresh, fresh)
    mult = bound_mult(add, fresh)
    fresh_budget = get_noise_budget(fresh, q)
    add_budget = get_noise_budget(add, q)
    mult_budget = get_noise_budget(mult, q)
    if (n > 2048):
        mod_switch = bound_mod_switch(n, t, q, p, mult)
        mod_switch_budget = get_noise_budget(mod_switch, p)
        return fresh_budget, add_budget, mult_budget, mod_switch_budget 
    else:
        return fresh_budget, add_budget, mult_budget


##############################################################################
# Estimates of worst-case and average-case heuristics for "bgv deep" circuit #
##############################################################################

# Variance of the output ciphertext after evaluating the "bgv deep" circuit
def variance_after_bgv_deep(n, t, q, p):
    fresh = variance_fresh(n, t)
    mult1 = variance_mult(fresh, fresh, n, t)
    mult2 = variance_mult(mult1, mult1, n, t)
    mult3 = variance_mult(mult2, mult2, n, t)
    return fresh, mult1, mult2, mult3

# Top-level function for noise budget predicted for average-case approach
def average_case_bgv_deep(n, t, q, p):
    fresh, mult1_var, mult2_var, mult3_var = variance_after_bgv_deep(n, t, q, p)
    fresh_bound =  alpha_bound_from_variance(fresh, n)
    mult1_bound = alpha_bound_from_variance(mult1_var, n)
    mult2_bound = alpha_bound_from_variance(mult2_var, n)
    mult3_bound = alpha_bound_from_variance(mult3_var, n)
    fresh_budget = get_noise_budget(fresh_bound, q)
    mult1_budget = get_noise_budget(mult1_bound, q)
    mult2_budget = get_noise_budget(mult2_bound, q)
    mult3_budget = get_noise_budget(mult3_bound, q)
    return fresh_budget, mult1_budget, mult2_budget, mult3_budget

# Top-level function for noise budget predicted for worst-case approach
def worst_case_bgv_deep(n, t, q):
    fresh = bound_fresh(n, t)
    mult1 = bound_mult(fresh, fresh)
    mult2 = bound_mult(mult1, mult1)
    mult3 = bound_mult(mult2, mult2)
    fresh_budget = get_noise_budget(fresh, q)
    mult1_budget = get_noise_budget(mult1, q)
    mult2_budget = get_noise_budget(mult2, q)
    mult3_budget = get_noise_budget(mult3, q)
    return fresh_budget, mult1_budget, mult2_budget, mult3_budget


######################################################
# Parameters used in implementations of the circuits #
######################################################

# Common parameters
n_2048 = 2048
n_4096 = 4096
n_8192 = 8192
n_16384 = 16384
n_32768 = 32768

# HElib parameters
t_helib = 3
q_2048_helib = 18014398492704769 # log q = 54
p_2048_helib = 0 # mod switch not supported for n = 2048
q_4096_helib = 649037106476272273878613017231361 # log q = 109
p_4096_helib = 3.1517442730074012e+16 # log p = 54.807
q_8192_helib = 1.1397723799332707e+66 # log q = 219.436
p_8192_helib = 3.559126845070405e+49 # log p = 164.606
q_16384_helib = 1.3894283839645433e+132  # log q = 438.969
p_16384_helib = 4.070677950511519e+115 # log p = 384.047

# SEAL parameters
t_4096_SEAL = 1032193 # log(t) = 19.9777 
t_8192_SEAL = 1032193 
t_16384_SEAL = 786433 # log(t) = 19.5854...
t_32768_SEAL = 786433 
q_4096_SEAL = 2**72
p_4096_SEAL = 2**36
q_8192_SEAL = 2**174
p_8192_SEAL = 2**130
q_16384_SEAL = 2**389
p_16384_SEAL = 2**340
q_32768_SEAL = 2**825
p_32768_SEAL = 2**770


###########################
# Generate results tables #
###########################

# Table 1: HElib, [CLP20] circuit, worst-case
print("HElib, [CLP20] circuit, worst-case (Table 1, column [CLP20]):")
print("n: " + str(n_2048))
print(worst_case_clp20(n_2048, t_helib, q_2048_helib, p_2048_helib))
print("n: " + str(n_4096))
print(worst_case_clp20(n_4096, t_helib, q_4096_helib, p_4096_helib))
print("n: " + str(n_8192))
print(worst_case_clp20(n_8192, t_helib, q_8192_helib, p_8192_helib))
print("n: " + str(n_16384))
print(worst_case_clp20(n_16384, t_helib, q_16384_helib, p_16384_helib))
print("\n")

# Table 1: HElib, [CLP20] circuit, average-case
print("HElib, [CLP20] circuit, average-case (Table 1, column Ours):")
print("n: " + str(n_2048))
print(average_case_clp20(n_2048, t_helib, q_2048_helib, p_2048_helib))
print("n: " + str(n_4096))
print(average_case_clp20(n_4096, t_helib, q_4096_helib, p_4096_helib))
print("n: " + str(n_8192))
print(average_case_clp20(n_8192, t_helib, q_8192_helib, p_8192_helib))
print("n: " + str(n_16384))
print(average_case_clp20(n_16384, t_helib, q_16384_helib, p_16384_helib))
print("\n")

# Table 2: HElib, "bgv deep" circuit, worst-case
print("HElib, bgv deep circuit, worst-case (Table 2, column [CLP20]):")
print("n: " + str(n_4096))
print(worst_case_bgv_deep(n_4096, t_helib, q_4096_helib))
print("n: " + str(n_8192))
print(worst_case_bgv_deep(n_8192, t_helib, q_8192_helib))
print("n: " + str(n_16384))
print(worst_case_bgv_deep(n_16384, t_helib, q_16384_helib))
print("\n")

# Table 2: HElib, "bgv deep" circuit, average-case
print("HElib, bgv deep circuit, average-case (Table 2, column Ours):")
print("n: " + str(n_4096))
print(average_case_bgv_deep(n_4096, t_helib, q_4096_helib, p_4096_helib))
print("n: " + str(n_8192))
print(average_case_bgv_deep(n_8192, t_helib, q_8192_helib, p_8192_helib))
print("n: " + str(n_16384))
print(average_case_bgv_deep(n_16384, t_helib, q_16384_helib, p_16384_helib))
print("\n")

# Table 3: SEAL, [CLP20] circuit, worst-case
print("SEAL, [CLP20] circuit, worst-case (Table 3, column [CLP20]):")
print("n: " + str(n_4096))
print(worst_case_clp20(n_4096, t_4096_SEAL, q_4096_SEAL, p_4096_SEAL))
print("n: " + str(n_8192))
print(worst_case_clp20(n_8192, t_8192_SEAL, q_8192_SEAL, p_8192_SEAL))
print("n: " + str(n_16384))
print(worst_case_clp20(n_16384, t_16384_SEAL, q_16384_SEAL, p_16384_SEAL))
print("n: " + str(n_32768))
print(worst_case_clp20(n_32768, t_16384_SEAL, q_32768_SEAL, p_32768_SEAL))
print("\n")

# Table 3: SEAL, [CLP20] circuit, average-case
print("SEAL, [CLP20] circuit, average-case (Table 3, column Ours):")
print("n: " + str(n_4096))
print(average_case_clp20(n_4096, t_4096_SEAL, q_4096_SEAL, p_4096_SEAL))
print("n: " + str(n_8192))
print(average_case_clp20(n_8192, t_8192_SEAL, q_8192_SEAL, p_8192_SEAL))
print("n: " + str(n_16384))
print(average_case_clp20(n_16384, t_16384_SEAL, q_16384_SEAL, p_16384_SEAL))
print("n: " + str(n_32768))
print(average_case_clp20(n_32768, t_32768_SEAL, q_32768_SEAL, p_32768_SEAL))
print("\n")

# Table 4: SEAL, "bgv deep" circuit, worst-case
print("SEAL, bgv deep circuit, worst-case (Table 4, column [CLP20]):")
print("n: " + str(n_16384))
print(worst_case_bgv_deep(n_16384, t_16384_SEAL, q_16384_SEAL))
print("n: " + str(n_32768))
print(worst_case_bgv_deep(n_32768, t_32768_SEAL, q_32768_SEAL))
print("\n")

# Table 4: SEAL, "bgv deep" circuit, average-case
print("SEAL, bgv deep circuit, average-case (Table 4, column Ours):")
print("n: " + str(n_16384))
print(average_case_bgv_deep(n_16384, t_16384_SEAL, q_16384_SEAL, p_16384_SEAL))
print("n: " + str(n_32768))
print(average_case_bgv_deep(n_32768, t_32768_SEAL, q_32768_SEAL, p_32768_SEAL))
print("\n")