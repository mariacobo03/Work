'''

@author: MARIA
'''
import math
import random

class MCMC(object):
    '''
    classdocs
    '''

    def __init__(self, values, initial_mean, sd):
        self.values = values
        self.mean = initial_mean
        self.sd = sd
        self.previous_log_likelihood = self.compute_loglikelihood_Values(self.mean)
    
    
    def log_pdf_normal(self,x, mean):
        return (-0.5*((x-mean)/self.sd)**2)-math.log(self.sd*(2*math.pi)**0.5)
    
    def compute_loglikelihood_Values(self, mean):
        log_likelihood = 0
        for v in self.values:
            likelihood = self.log_pdf_normal(v, mean)
            log_likelihood += likelihood
        return log_likelihood
    
    def next_movement(self,sd_movement):
        return random.gauss(self.mean, sd_movement)
    
    def set_mean(self,new_mean):
        self.mean = new_mean
    
    def get_mean(self):
        return self.mean
    
    def metropolis_algorithm(self, sd_movement):
        new_mean = self.next_movement(sd_movement)
        log_likelihood_new_mean = self.compute_loglikelihood_Values(new_mean)
        ratio = log_likelihood_new_mean-self.previous_log_likelihood
        shall_I_move_to_new_mean = math.log(random.uniform(0,1)) < ratio
        if(shall_I_move_to_new_mean):
            self.previous_log_likelihood = log_likelihood_new_mean
            self.mean = new_mean

        return self.mean         
    
    
def main():
    
    observed_seq = []
    
    mean = 5
    sd = 1
    
    for i in range(100):
        observed_seq.append(random.gauss(mean,sd))    
  
    mcmc_seq = []
    
    initial_mean =  random.uniform(-10,10)
    
   
    mcmc = MCMC(observed_seq, initial_mean, 1.0)
    
    for i in range(1000):
        print(mcmc.metropolis_algorithm(0.5))
    

    
    
if __name__ == "__main__":
    main () 
    
        


        
