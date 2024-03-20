from random_multinominal import RandomMultinomial
import math
# RandomMultinomial generates random numbers given a distribution


class HiddenMarkovModel(object):  # class to code methods that are used in HMM

    '''
    Create an object HiddenMarkovModel with transition states and emission probabilities for each state
    Transition is a dictionary where, for each state, we count the probability to move to another state
    Emission is a dictionary where, for each state, we store the probabilities of each category
    '''
    def __init__(self, transition_probabilities, emission_probabilities):
        # if the length of transition and emission probabilities are not the same then
        if len(transition_probabilities) != len(emission_probabilities):
            raise Exception("For each state, we must have an emission probability vector, but found " + len(transition_probabilities) + " " + len(emission_probabilities))
        self.n = len(transition_probabilities)  # take the length of transition probabilities
        self.transition_probabilities = transition_probabilities  # transition matrix
        self.emission_probabilities = emission_probabilities  # emission matrix
        self.random_transition = {}  # dictionary of random transition probabilities
        self.random_emission = {}  # dictionary of random emission probabilities

        for key in self.transition_probabilities:  # iterate over de transition probabilities
            # create a random multinomial distribution for the emission and transition probabilities
            self.random_emission[key] = RandomMultinomial(list(emission_probabilities[key].values()))
            self.random_transition[key] = RandomMultinomial(list(transition_probabilities[key].values()))

    '''
    Compute a sequence of length using a prior probabilities dictionary of the states to start the chain.
    Return the sequence and its hidden states
    '''

    def generate_sequence(self, length_sequence, prior_probabilities):
        seq = []  # list where we will store the categories
        hidden_states = []  # list where we will store the states
        # iterate over the i elements of the sequence, generating the emission and transition probabilities
        # create a random multinomial distribution of the prior probabilities
        pr = RandomMultinomial(list(prior_probabilities.values()))
        state = list(prior_probabilities.keys())[pr.sample()]  #
        for i in range(length_sequence):  # iterate over the length of the sequence
            # from the state, use its emission probability to sample one element
            category = list(self.emission_probabilities[state].keys())[self.random_emission[state].sample()]
            seq.append(category)  # add to the seq list
            hidden_states.append(state)  # add to the hidden_states list
            state = list(self.transition_probabilities[state].keys())[self.random_transition[state].sample()]
        # the putput will be a list of the sequence generated and the assigned hidden states
        output = [seq, hidden_states]
        return output  # returns the list of the seq and the corresponding states


    '''
    Compute the log probability of the sequence using forward algorithm. It uses the scale algorithm (Rabiner) to prevent underflow
    '''

    # function that takes as input an observed sequence of categories and returns the log-likelihood
    # to calculate the log-likelihood we need the transition probabilities
    def log_probability_sequence_using_scaling(self, seq, prior_probabilities):
        # scaling factor
        p_i = prior_probabilities.copy()  # replicate de prior probabilities list and leave the other untouched
        s = 0  # create count, this will be the sum of the probabilities of the events
        for i in range(len(seq)):  # iterate over the length of the sequence
            s1 = 0  # new count which is the variable to start scaling
            p_i1 = {}  # create an empty dictionary, for the probabilities
            for key2 in p_i:  # iterate over the keys in the p_i dictionary
                a = 0  # initialize the probability every time it goes through the for loop
                for key1 in p_i:  # iterate over the keys in the p_i dictionary
                    # updating the probability
                    # enter the dictionary of dictionaries
                    # multiply the transition probability by the prior probability
                    a += self.transition_probabilities[key1][key2] * p_i[key1]
                # multiply the emission probability by the one from the previous step
                a *= self.emission_probabilities[key2][seq[i]]
                p_i1[key2] = a  # update  the p_i1 dictionary with the probability of the event
                s1 += a  # add to the scaling factor the probability we just calculated
        for key2 in p_i1:  # iterate over the dictionary of probabilities
            # change value of the dictionary,
            p_i[key2] = p_i[key2]/s1  # divide the prior probability value by the scaling factor
        # we transform into logarithm to avoid the probability getting closer to zero
        s += math.log(s1)  # use log function of the math library to do the logarithm
        return s  # log probability

    '''
    Compute the log probability of the sequence using forward algorithm.
     '''

    # function that takes as input an observed sequence of categories
    # returns the log-likelihood, using the scale approach
    def log_probability_sequence_without_scaling(self, seq, prior_probabilities):
        # without scaling factor
        p_i = prior_probabilities.copy()  # replicate de prior probabilities list and leave the other one untouched
        for i in range(len(seq)):  # iterate over the length of the sequence
            p_i1 = {}  # create an empty dictionary
            for key2 in p_i:  # iterate over the keys in the p_i dictionary
                a = 0  # initialize the probability every time it goes through the for loop
                for key1 in p_i:  # iterate over the keys in the p_i dictionary
                    # updating the probability
                    # enter the dictionary of dictionaries
                    # multiply the transition probability by the prior probability
                    a += self.transition_probabilities[key1][key2] * p_i[key1]
                    # multiply the emission probability by the one from the previous step
                a *= self.emission_probabilities[key2][seq[i]]
                p_i1[key2] = a  # update  the p_i1 dictionary with the probability of the event
                p_i = p_i1  # inserting a dictionary inside another one
        s = 0  # count variable where we store the sum of probabilities
        for key2 in p_i1:  # iterate over every key
            s += p_i1[key2]  # add to the sum of probabilities, the total probability of the sequence
        # we transform into logarithm to avoid the probability getting closer to zero
        s = math.log(s)  # use log function of the math library to do the logarithm
        # very small result because there is no scaling
        return s  # return the log probability of the sequence

    '''
    given a sequence generated by a HiddenMarkovModel object, estimate the emission probabilities of the different categories for each state
    '''

    # the emission probabilities depend on the state probabilities
    def estimate_emission_probabilities(self, seq):
        emission_prob = []  # list where we will store the emission probabilities
        for i in range(len(seq)-1):  # iterate over the length of the sequence
            if i == 0:  # if the x = 0, which means that it is the first time it appears
                # we multiply the prior probability with the emission probability
                self.prob = 0.5 * emission_prob[seq[i]][seq[i]]  # create variable to store multiplication
                emission_prob.append(self.prob)  # add the probability inside the emission list
            self.s = self.transition_probabilities[seq[i]]  # s are the state probabilities
            if seq[i] == seq[i+1]:  # if the state where we are is the same as the following then
                emission_prob.append(self.s[seq[i]])  # append the state where we are to the emission_prob list
            else:  # if they are not the same
                emission_prob.append(self.s[seq[i+1]])  # append the following state (i + 1) to the list
        return emission_prob  # return the list of emission probabilities


    '''
    given a sequence generated by a HiddenMarkovModel object, estimate the transition probabilities between states
    '''

    # transition probability is probability of moving from one state to another
    def estimate_transition_probabilities(self, seq):  # seq is a list of the hidden states
        transition_prob = []  # list where we will store the transition probabilities
        for i in range(len(seq)-1):  # iterate over the length of the sequence
            if i == 0:  # if the x = 0, which means that it is the first time it appears
                transition_prob.append(0.5)  # add to the list the prior probability
            self.s = self.transition_probabilities[seq[i]]  # s are the state probabilities
            if seq[i] == seq[i+1]:  # if the state where we are is the same as the following then
                transition_prob.append(self.s[seq[i]])  # append the state where we are to the transition_prob list
            else:  # if they are not the same
                transition_prob.append(self.s[seq[i+1]])  # append the following state (i + 1)to the list
        return transition_prob  # return the transition_prob list


def main():
    # probabilities of changing from one state to another
    # these are dictionaries of dictionaries
    # fast and slow evolving 
    transition_probabilities = {"S": {"S": 0.9, "L": 0.1}, "L": {"S": 0.2, "L": 0.8}}
    # probabilities of each category depending on the state
    emission_probabilities = {"S": {"G": 0.1, "C": 0.2, "T": 0.2, "A": 0.5}, "L": {"G": 0.01, "C": 0.1, "T": 0.3, "A": 0.59}}
    # hidden markov model (hmm) is the class
    hmm = HiddenMarkovModel(transition_probabilities, emission_probabilities)
    # initial probability of starting in S or L state
    # selects the one in position 0 that is the category, not the state
    prior_probabilities = {"S": 0.5, "L": 0.5}
    # generates a sequence of length 100 with the prior probabilities
    seq = hmm.generate_sequence(100, prior_probabilities)[0]  # list of categories, with [1] would be a list of states
    print(seq)  # new random sequence generated each time we run the main function
    print(hmm.log_probability_sequence_using_scaling(seq, prior_probabilities))  # log probabilities are always negative
    print(hmm.log_probability_sequence_without_scaling(seq, prior_probabilities))  # print logarithm without scaling
    # looks at the states of the sequence
    print(hmm.estimate_transition_probabilities(hmm.generate_sequence(100, prior_probabilities)[1]))
    # emission probabilities for the whole sequence
    # print(hmm.estimate_emission_probabilities(hmm.generate_sequence(100, prior_probabilities)))


if __name__ == "__main__":
    main()
