############################################################################################################################
# 10/01/2018
# Reproduction of SPN model of Mura and Nagy, 2008, with the addition of deterministic evolution for the mass.
#
# Author: Lorenzo Federico Signorini
#
# For further reference, see 'Stochastic Modelling for Systems Biology', Wilikinson, 2006, pg. 24-34 and 182-185
############################################################################################################################

from random import *
import numpy as np
from math import *
import matplotlib.pyplot as plt


#############################################
# UTILITY FUNCTIONS (all called inside 
# Gillespie algorithm) :
#############################################

def update_rates(P,p,m):
	# Calculates current value of cycb, and current rates
        # (also called in the initialization step!)
	
	# cycb:
	
	num = 2*(pow(alpha,2))*P['cycb_T']*P['CKI_T']                             # numerator
	sigma = alpha*P['cycb_T'] + alpha*P['CKI_T'] + pow(keq,-1)                # sigma factor
	den = sigma + sqrt(pow(sigma,2) - 4*pow(alpha,2)*P['cycb_T']*P['CKI_T'])  # denominator
	cycb = alpha*P['cycb_T'] - num / den

	# Rates:
	
	r= np.array([k1/alpha, k21*P['cycb_T'], k22*P['CDH1_A']*P['cycb_T']*alpha ,  k23*P['CDC20_A']*P['cycb_T']*alpha,\
		     k31*P['CDH1_I']/(J_three + P['CDH1_I']*alpha),\
		     k32*P['CDC20_A']*P['CDH1_I']*alpha/(J_three + P['CDH1_I']*alpha),\
		     k41*m*cycb*P['CDH1_A']/(J_four + P['CDH1_A']*alpha) ,(k42 * P['SK'] * P['CDH1_A']) * alpha/\
		     (J_four + P['CDH1_A']*alpha), k51/alpha, (k52*(pow(m*cycb,n_exp))/(pow(J_five,n_exp) + pow(m*cycb,n_exp)))\
		     /alpha , k6*P['CDC20_I'],k7*P['IE_A']*P['CDC20_I']*alpha/(J_seven + P['CDC20_I']*alpha) ,\
		     alpha*k8*P['CDC20_A']/(J_eight + P['CDC20_A']*alpha), k6*P['CDC20_A'], k9*m*cycb*P['IE_I'] ,\
		     k10*P['IE_A'], k11/alpha , k121*P['CKI_T'] , k122*P['SK']*P['CKI_T']*alpha , k123*m*cycb*P['CKI_T'],\
			k131/alpha , k132*P['TF_A'] , k14*P['SK'], k151*m*P['TF_I']/(J_fifteen + P['TF_I']*alpha),\
		     k152*P['SK']*P['TF_I']*alpha/(J_fifteen + P['TF_I']*alpha), k161*P['TF_A']/(J_sixteen + P['TF_A']*alpha),\
		     k162*m*cycb*P['TF_A']/(J_sixteen + P['TF_A']*alpha)])
	return r, cycb


def init_X(u, N_max, ordered_places_list):
	# Creates a u* N_max matrix to store values of tokens of all the u species, 
	# for N_max iterations of Gillespie algorithm.
	
	X = np.zeros((u,N_max)) 	# Solutions matrix: every row will be filled with a different species.
	M = [P[k] for k in ordered_places_list]
	X[:,0] = M  # Fill first column with initial conditions
	return X


def update_P(new_values,P, ordered_places_list):
	# Updates the current markings of P with the newly calculated ones.
	
	for i in range(len(new_values)): ## TODO DICTIONARIES ARE NOT SORTED DIOCANE!! quindi l'indice quando prendo le keys del dizionario chissa' se le prendo in ordine! magari si, ma magari no!
		key = ordered_places_list[i]
		P[(ordered_places_list[i])] = new_values[i]
	return P


def check_token(mass , curr_cycb , low): 
	if low == 1 and mass*curr_cycb > 0.2:
		print("\n\ntransition to HIGH\n\n")
		return("to_high")
	elif low == 0 and mass*curr_cycb < 0.1:
		print("DIVISION and comeback to low")
		return("to_low")
	else:
		return 0

	
def update_mass(m, p ,t , low, cycb,t_div, init_mass):
	C = init_mass / (1- init_mass/mass_max)
	m=(C * np.exp(mass_gr * (t-t_div) ))/(1 + (C/mass_max)*np.exp(mass_gr * (t-t_div)))
	response = check_token(m , cycb, low)
	if response == "to_high":
		low = 0
	elif response == "to_low":
		low = 1
		m = m/2
		init_mass = m
		t_div = t
		print("DIVISION OCCURRED")
	return m, low, t_div, init_mass


def fire_transition(P, Pre, iter_index, rand_index, X, S, ordered_places_list):
	# Checks that the tokens of the species involved in the Pre fase of the selected reaction
	# are > 0. If the check is passed, fires transition.
	
	t = Pre[rand_index]                                  # (rand_index+1)th row of Pre matrix
	transition_has_fired = False
	for j in range(len(t)):
		fire_transition = 0
		fire_indeed = 0
		if t[j] == 1:                                 # if there exists a Pre edge for ith species in this transition
			key=ordered_places_list[j]
			species = P[key]
			fire_transition += 1
			if species > 0:                       # if ith species is > 0
				fire_indeed += 1
	if fire_transition == fire_indeed:                    # = if all species present in t are > 0
		X[:,iter_index] = X[:,iter_index-1]+S[:,rand_index] # Calculate new markings according to selected transition
		transition_has_fired = True
	return X, transition_has_fired


def prints(x,m,t,cycb,alpha,low):
	print(low)
	print("cycb= ",cycb)
	print("m= ", m)
	print("cycb_T= ", x[0]*alpha)
	print("time : " , t)
	return


#############################################
#        PETRI NET INITIALIZATION
#############################################

def initialize_petri_net():
	# This Petri Net, N, is defined as a set such that:
        # N = [P, ordered_places_list, T, Pre, Post, p, r, m, cycb]
        # where P = #places+markings; T = transitions; Pre; Post; p=parameters( aka 'propensities'); r=rates; m=mass
        
	# Initial state:
	P = {'cycb_T':97, 'CDH1_I':418, 'CDH1_A':5.0, 'CDC20_I':24, 'CDC20_A':0, 'IE_I':384, 'IE_A':40, 'CKI_T':25, 'SK':39, 'TF_I':408, 'TF_A':15} #11 species (+ cycb) inital state ('marking') 
	ordered_places_list = ['cycb_T','CDH1_I','CDH1_A','CDC20_I','CDC20_A','IE_I','IE_A','CKI_T','SK', 'TF_I', 'TF_A']
	T = ("t1", "t21"  , "t22"  , "t23"  , "t31"  , "t32"  , "t41"  , "t42"  , "t51"  , "t52" ,
	"t6a"  , "t7"  , "t8"  , "t6b"  , "t9"  , "t10"  , "t11" , "t121" , "t122" , "t123"  , "t131" , "t132" ,
	"t14"  , "t151" , "t152","t161","t162") #27 transitions
        
	# Transition propensities:
	# First cycle:
	alpha, k1, k21, k22, k23 = 0.00236012, 0.04, 0.04, 1.0, 1.0  # cycb- mutant: k1 = 0. 03 , k22 = 0, k23 = 0.2
	# Second cycle:
	k31, k32, k41, k42, J_three, J_four, init_mass, keq = 1.0 , 10.0, 35.0, 2.0, 0.04,  0.04, 0.704045,  1000.0 
	# Third cycle:
	k51 , k52 , J_five , n_exp , k6 =  0.0050, 0.2 , 0.3 , 4.0 , 0.1 
	# Fourth cycle:
	k7, J_seven , k8 , J_eight = 1.0 , 0.0010 , 0.5 , 0.0010
	# Fifth cycle:
	k9, k10 = 0.1 , 0.02
	# Sixth cycle:
	k11 , k121 , k122 , k123 = 1.0 , 0.2 , 50.0 , 100.0  # CKI- mutant: k11 = 0
	# Seventh cycle:
	k131, k132 , k14 = 0.0 , 1.0 , 1.0   # SK- mutant (cln1delta, cln2delta, cln3delta triple mutant) : k131, 
	# Eighth cycle:
	k151 , J_fifteen , k152 , k161 , J_sixteen , k162 = 1.5 , 0.01 , 0.05, 1.0 , 0.01 , 3.0
	# Mass:
	mass_gr, mass_max = 0.01, 10.0   # cycb- mutant: 
	p = [alpha, k1, k21, k22, k23,k31, k32, k41, k42, J_three, J_four, init_mass, keq,k51 , k52 , J_five ,\
	     n_exp , k6,k7, J_seven , k8 , J_eight,k9, k10,k11 , k121 , k122 , k123, k131, k132 , k14, k151 , J_fifteen ,\
	     k152 , k161 , J_sixteen , k162, mass_gr, mass_max] # transition propensities

	# Calculate initial cycb and transition rates:
	r, cycb = update_rates(P,p,init_mass)                   # transition rates 

	# 11x27 matrices:
	Pre = np.array([[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,      #   t1
	[ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t21
	[ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t22
	[ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t23
	[ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t31
	[ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t32
	[ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t41
	[ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t42
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t51
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t52
	[ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t6a
	[ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t7
	[ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t8
	[ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t6b
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t9
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t10
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t11
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] ,                      #   t121
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] ,                      #   t122
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] ,                      #   t123
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t132
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t131
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] ,                      #   t14
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] ,                      #   t151
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] ,                      #   t152
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] ,                      #   t161
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])                      #   t162

	Post = np.array([[ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,     #   t1
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t21
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t22
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t23
	[ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t31
	[ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t32
	[ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t41
	[ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t42
	[ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t51
	[ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t52
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t6a
	[ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t7
	[ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t8
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t6b
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t9
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t10
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] ,                      #   t11
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t121
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t122
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t123
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] ,                      #   t131
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] ,                      #   t132
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,                      #   t14
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] ,                      #   t151
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0] ,                      #   t152
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0] ,                      #   t161
	[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]])                      #   t162    

	return P, T, Pre, Post, p, r, cycb, N_max, ordered_places_list, init_mass


#############################################
#          GILLESPIE ALGORITHM CYCLE
#############################################

def gillespie(P, T, Pre, Post, p, r, cycb, N_max, ordered_places_list, init_mass):
	T_tot, CYCB, M = np.zeros(N_max), np.zeros(N_max), np.zeros(N_max) # time, cycb, mass
	t = 0
	tp = 0 # printing time utility
	S = np.transpose(Post - Pre) 	# Stoichiometry matrix S (transposed for simplicity (see reference book))
	u = S.shape[0]  # n rows
	v = S.shape[1]  # n columns
	X = init_X(u, N_max, ordered_places_list)
	CYCB[0] = cycb # initial value for cycb
	M[0] = init_mass # initial mass 
	low = 1 
	t_div = 0
	for i in range(1,N_max):

		# First random number generation: define the time-step:

		r_sum=sum(r)
		rand = random()
		s = - np.log(rand)/r_sum
		t = t + s
		if t == float('nan'):
			print t
			raise ValueError("Time 'Nan', SOMETHING IS WRONG")
		T_tot[i] = t

		# Deterministic step: calculating mass:

		M[i], low, t_div, init_mass = update_mass(M[i-1], p, t, low, cycb, t_div, init_mass) #Check for division and update mass
		
		# Second random number generation: random decision on transition:

		pi = r/r_sum # converting rates to probabilities (normalization)
		if len(pi[pi<0]) != 0:
			print pi
			raise ValueError("NEGATIVE PROBABILITIES, SOMETHING IS WRONG")
		rand_index = np.random.choice(len(pi),p=pi) # Quasi random decision on transition
		X, transition_has_fired = fire_transition(P, Pre, i, rand_index, X, S, ordered_places_list)

		# Updating the Petri Net according to fired transition:

		if transition_has_fired:
			P = update_P(X[:,i],P, ordered_places_list)
			r, cycb = update_rates(P, p, M[i])
		CYCB[i] = cycb

		# Printing useful info every 1000 steps:

		if tp == 1000: 
			prints(X[:,i],M[i],t,cycb, alpha, low)
			tp = 0
		tp+=1
	return X, T_tot, CYCB, M


#############################################
#             PLOTTING FUNCTION
#############################################

def my_plot(X, P, T_tot, CYCB, M, ordered_places_list):
	for i in range(len(X[:,0])): 
		plt.plot(T_tot, X[i]*alpha , label=ordered_places_list[i])
	plt.plot(T_tot, M, label = "mass", color = "black")
	plt.legend()
	plt.show()
	return


#############################################
#             EXECUTE PROGRAM               
#############################################

if __name__ == "__main__":
	N_max = input("\n>Define the number of iterations you would like to simulate...\n-> 1244567*3 iterations are around 1400 time units\n")
	P, T, Pre, Post, p, r, cycb, N_max, ordered_places_list, init_mass = initialize_petri_net()
	X, T_tot, CYCB, M = gillespie(P, T, Pre, Post, p, r, cycb, N_max, ordered_places_list, init_mass)
	my_plot(X, P, T_tot, CYCB, M, ordered_places_list)
