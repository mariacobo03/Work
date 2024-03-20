'''
Class to store a symmetric distance matrix

@author: MARIA
'''

from DistanceMatrix import DistanceMatrix


class DistanceBasedAlgorithms(object):
    '''
    classdocs
    '''

    def __init__(self, d_matrix):
        '''
        Constructor
        '''
        self.d_matrix = d_matrix # distance matrix

    def UPGMA(self):
        new_matrix = self.d_matrix.copy()  # copy the distance matrix, leave the original untouched
        while new_matrix.n_rows() != 1:  # while the number of rows is not equal to 1
            # check lowest value of distance matrix
            lowest_i = 0  # the lowest i value will be 0
            lowest_j = 1  # the lowest j value will be 1
            # pick i and j with the smallest distance
            min_dist = new_matrix.get_value_i_j(lowest_i, lowest_j)
            for i in range(new_matrix.n_rows()-1):  # iterating over the number of rows - 1
                for j in range(i+1, new_matrix.n_rows()):  # one position over, iterating till the number of rows
                    # if the smallest distance is bigger than the new matrix values of i and j
                    if(min_dist > new_matrix.get_value_i_j(i, j)):
                        # create new smallest distance matrix
                        min_dist = new_matrix.get_value_i_j(i, j)
                        lowest_i = i  # lowest i will be value i
                        lowest_j = j  # lowest j will be value j
            # calculate average distance
            # iterate over the length of the number of rows
            for k in range(new_matrix.n_rows()):
                # calculate the average of the new matrix value of i and j
                avg = (new_matrix.get_value_i_j(lowest_i, k) + new_matrix.get_value_i_j(lowest_j, k))/2.0
                # compute new distance matrices
                new_matrix.add_value_i_j(lowest_i, k, avg)
            new_matrix.add_value_i_j(lowest_i, lowest_i, 0)
            # apply the newick format to return the results
            new_matrix.change_name_species_i(lowest_i, ("(" + new_matrix.get_name_of_species()[lowest_i] + ":" + str(min_dist/2.0) + ":" + new_matrix.get_name_of_species()[lowest_j] + ":" + str(min_dist/2.0) + ")"))
            # remove the species i
            new_matrix.remove_species_i(lowest_j)
        # will return the new matrix
        return new_matrix.get_name_of_species()[0]

    '''
    PSEUDOCODE:
    M → keys = names of nodes 
    while number of rows of M > 1
        pick i, j < distance in M (pick i, j with smallest distance)
        create a new label (new key) Newick (key_i:distance/2, key_j: distance/2)
        create a new matrix, nM = number of rows - 1, number of columns -1)
        for i in n row M
            for j in n now M
                if (l != i, k != j)
                    compute new distance i, l and j, k
                    nmM[i-1, j-1] = new distance
    M = nM
    '''

    # neighbor joining algorithm
    def NJ(self):
        new_matrix = self.d_matrix.copy()  # construct new matrix from distance matrix
        # this total_dist will be the sum of all the other distances
        total_dist = []  # create a vector where we will add the total distance

        for i in range(new_matrix.n_rows()):  # iterate over the rows of the matrix
            total_dist.append(sum(new_matrix.m[i]))  # sum the distances of the first row
        dnj_matrix = new_matrix.copy()  # make another copy, we can later on compare with UPGMA

        while new_matrix.n_rows() > 2:  # we while iterate while we have more than two clusters
            for i in range(new_matrix.n_rows()-1):  # iterating over the number of rows - 1
                for j in range(i+1, new_matrix.n_rows()):  # one position over, iterating till the number of rows
                    # calculate the distance from i to j
                    # compute D*ij
                    distance_ij = (dnj_matrix.n_rows()-2 * new_matrix.get_value_i_j(i, j) - total_dist[i] - total_dist[j])
                    dnj_matrix.add_value_i_j(i, j, distance_ij)  # ij and ji are symmetrical, value of Dij
            # find the smallest element D*ij of D*
            min_dist = new_matrix.get_value_i_j(0, 1)  # first value
            min_position = [0, 1]  # determinate the smallest positon
            for i in range(dnj_matrix.n_rows() - 1):  # iterating over the number of rows - 1
                for j in range(i + 1, dnj_matrix.n_rows()):  # one position over, iterating till the number of rows
                    # if value in the position we are is smaller to the minimum value
                    if min_dist > dnj_matrix.get_value_i_j(i, j):
                        min_dist = dnj_matrix.get_value_i_j(i, j)  # that value becomes the new minimum
                        min_position = [i, j]  # save the positions of the min value

            # Compute total_dist = (total_dist_i - total_dist_j) / new_matrix.n_rows() - 2.0

        return dnj_matrix


def main():
    species_names = evolution.get_list_of_species_name()
    distance_matrix = DistanceBasedAlgorithms(species_names)

    print(DistanceBasedAlgorithms(new_matrix).NJ())
    print(DistanceBasedAlgorithms(dnj_matrix).UPGMA())


if __name__ == "__main__":
    main()
