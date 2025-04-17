

from re import sub
import numpy as np
import pandas as pd
from itertools import combinations, combinations_with_replacement, count
from collections import defaultdict
import math
from scipy.stats import norm

 

def hello_world():
    print("Goodbye Cruel WOrld")



class DAFset:

    def __init__(self, df: pd.DataFrame, histones=False):
        self.raw_df = df
        self.data_df = pd.DataFrame(
            data={
                'num_blocked' : df.apply(lambda row: row.value_counts().get(-1, 0), axis=1),
                'num_bound' : df.apply(lambda row: row.value_counts().get(1, 0), axis=1)
            },
            index = df.index
        )
        self.hist_mask = self.data_df['num_blocked']>0
        self.hist_I = df.index[self.hist_mask]
        self.free_df = df[~self.hist_mask]
        self.length = df.shape[1]
        self.h = histones
    


    #### Methods

    def num_reads(self, a: np.array) -> int:
        """Counts the number of matching instances in the dataframe from a given 1-dimensional array and return an integer"""

        df = self.raw_df if self.h else self.free_df

        try:
            a = a if type(a)==np.array else np.array(a)
            assert a.shape == tuple([self.length])
        except:
            print("ERROR: Passed value could not be turned into a numpy array")
            return None
        
        n = np.count_nonzero((df.values == a).all(axis=1))

        return n

  
    def num_reads_sites(self, sites) -> int:
        """Counts the number of matching instances of bound sites in the dataframe.  Automtically converts a list of sites into the correct numpy array
        for matching with the dataframe and returns an integer"""

        df = self.raw_df if self.h else self.free_df

        a = np.zeros(self.length)
        np.put(a,[m-1 for m in sites],1)
        n = np.count_nonzero((df.values == a).all(axis=1))
        return n


    def num_contains(self, sites) -> int:
        """Similair to the num_reads family of functions.  This takes a list of sites and counts the number of instances in which the given sites
        are bound, regardless of the binding status of non-specified sites and returns a count integer"""

        df = self.raw_df if self.h else self.free_df
        
        return (df[sites].sum(axis=1).values == len(sites)).sum()


    def get_empty(self) -> int:
        """Simple function to count the number of unbound reads in the datababase, which is used as a reference state to work around an intractable parition?
        function and returns a count integer"""

        return self.num_reads(np.zeros(self.length))


    def combination_contains(self,sites):
        """similair to the num_contains function, but performs for all possible combinations of the specified list of sites and returns a dict of ints"""

        count_dict = dict()

        for i in range(len(sites)):
            site_maps = list(combinations(sites,i+1))
            for map in site_maps:
                count_dict[map] = self.num_contains(list(map))

        return count_dict


    def combination_counts(self,sites):
        """similair to the num_reads_sites function, but performs for all possible combinations of the specified list of sites and returns a dict of ints"""

        count_dict = dict()

        count_dict[tuple([0])] = self.get_empty()
        if 0 in sites:
            sites.remove(0)

        for i in range(len(sites)):
            site_maps = list(combinations(sites,i+1))
            for map in site_maps:
                bound_rep = np.zeros(self.length, dtype=int)
                np.put(bound_rep,[m-1 for m in map],1)
                count_dict[map] = self.num_reads(bound_rep)

        return count_dict


    def boltzz_factors(self, count_dict):
        '''Takes a state count dictionary and calculates the Boltzmann factor of each state relative to 
            a the unbound reference state'''


        '''Divide each state count by a reference state count (default is unbound state)
            to get the boltzmann factor for each state'''
        factor_dict = dict([(tuple([0]), 1)] + 
                           [(key, value / self.get_empty()) for key, value in count_dict.items()])
        
        return factor_dict

 
    def energy_states(self, factor_dict):
        ''' Takes in a dictionary of state Boltzmann Factors and calculates a relative energy value for
            each state and returns a dictionary of values proportionate to the true delta G'''


        factor_tups = sorted(factor_dict.keys(), key=len)
        energy_dict = dict()

        for k in factor_tups:
            # Sum up all sub states of the state
            sub_states = sum([value for key, value in energy_dict.items() if set(key).issubset(set(k)) and k != key])
            energy_dict[k] =  -1*math.log(factor_dict[k]) - sub_states 

        return energy_dict

    def p_val(self, sites):
        ''' Calculates the expected counts of a state given there is
            no additional binding interaction at the highest level and returns float'''
        
        count_dict = self.combination_counts(sites)
        factor_dict = self.boltzz_factors(count_dict)
        energy_dict = self.energy_states(factor_dict)
        del energy_dict[tuple(sites)]
        sub_energies = -1*sum(list(energy_dict.values()))
        exp_count = self.get_empty() * math.exp(sub_energies)

        c = math.comb(10,len(sites))
        n = sum(list(count_dict.values()))
        p = exp_count / n
        stdev = math.sqrt(n*p*(1-p))
        Zscore = (exp_count-count_dict[tuple(sites)])/stdev
        p_score = c*2*norm.sf(abs(Zscore))

        # left_tail = False if Zscore < 0 else True

        return {
            'p': p_score,
            'z': Zscore,
            'c': exp_count,
            'n': count_dict[tuple(sites)]
        }
    


