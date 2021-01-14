import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

###########
# CLASSES #
###########


class GaussInt():
    def __init__(self, a, b):
        '''GaussInt(a,b) => a+bi'''
        if not (isinstance(a, int) and isinstance(b, int)):
            raise Exception("GaussInt must be given two integers")
        self.a = a
        self.b = b

        self.var_isprime = None
        self.var_norm = None

    def __add__(self, other):
        '''Standard addition of Gaussian integers'''
        return GaussInt(self.a + other.a, self.b + other.b)

    def __str__(self):
        '''print(GaussInt(a,b)) prints "(a,b)"'''
        return "({},{})".format(self.a, self.b)

    def norm(self):
        '''Returns standard norm: a**2 + b**2'''
        if self.var_norm is None:
            self.var_norm = self.a**2 + self.b**2
        return self.var_norm

    def isprime(self):
        '''Tests primality of a Gaussian integer. Returns boolean'''
        if self.var_isprime is None:
            if self.a == 0 and self.b == 0:
                self.var_isprime = False
            elif self.a == 1 and self.b == 0:
                self.var_isprime = False
            elif self.a == 1 and self.b == 1:
                self.var_isprime = True
            elif self.b == 0 and (self.a % 4) == 3 and sp.isprime(self.a):
                self.var_isprime = True
            else:
                self.var_isprime = sp.isprime(self.norm())
        return self.var_isprime

#####################
# UTILITY FUNCTIONS #
#####################


def gauss_to_coord_data(ints, prime_only=False):
    '''Converts a list of GaussInt objects into two lists of x and y data for plotting'''
    x = []
    y = []
    for cur in ints:
        if not prime_only or cur.isprime():
            x.append(cur.a)
            y.append(cur.b)
    return x, y

##########################
# MATHEMATICAL FUNCTIONS #
##########################


def ti_no_waste(r, ints=None):
    '''Counts primality using a given 2D list of pre-calculated Gaussian integers.
    ints is a 2D list holding GaussInts, i.e.
    ints[4][7] is GaussInt(4,7) which already has norm, isprime calculated
    '''

    if ints is None:
        ints = [[GaussInt(0, 0), GaussInt(0, 1)], [
            GaussInt(1, 0), GaussInt(1, 1)]]

    orig_min_a = len(ints)
    orig_min_b = len(ints[0])

    sum = 0

    # counting for the ints we have already calculated
    for x in ints:
        for y in x:
            if y.norm() <= r:
                sum += y.isprime()

    # calculating and counting new ints
    for a in range(0, int(np.sqrt(r)) + 1):
        if a >= orig_min_a:
            ints.append([])
            min_b = 0
        else:
            min_b = orig_min_b
        for b in range(min_b, int(np.sqrt(r)) + 1):
            d = GaussInt(a, b)
            if (d.norm() <= r) and d.isprime():
                sum += 1
            ints[a].append(d)
    return sum, ints


def ti(r):
    '''ti(r) counts the number of Gaussian integers with norm less than r'''
    sum = 0
    for a in range(0, int(np.sqrt(r)) + 1):
        for b in range(0, int(np.sqrt(r)) + 1):
            d = GaussInt(a, b)
            if (d.norm() <= r) and d.isprime():
                sum += 1
    return sum


def plot_ti(max_r, interval=.1, filename="prime_count.png", save=False):
    '''Creates a plot of ti(r) over r'''
    t = np.arange(0, max_r, interval)
    s = list()

    _, ints = ti_no_waste(0)
    for cur_num in t:
        cur_sum, ints = ti_no_waste(cur_num, ints)
        s.append(cur_sum)

    fig, ax = plt.subplots()
    ax.plot(t, s, color='green')

    ax.set(xlabel='Radius r', ylabel='Ti(r)',
           title='Count of prime Gaussian integers with norm < r')
    ax.grid()

    if save:
        fig.savefig(filename)
    plt.show()


def plot_primes(size, filename="prime_plot.png", save=False):
    _, ints = ti_no_waste(size)
    x = []
    y = []
    for row in ints:
        cur_x, cur_y = gauss_to_coord_data(row, prime_only=True)
        x.extend(cur_x)
        y.extend(cur_y)
    plt.scatter(x, y)
    if save:
        plt.savefig(filename)
    plt.show()

#########
# DEBUG #
#########


def debug_gauss_int():
    '''Just a simple verification that things are working as expected'''
    for x in range(0, 10):
        for y in range(0, 10):
            cur = GaussInt(x, y)
            print("{}+i{}: norm: {} isprime: {}".format(x,
                                                        y, cur.norm(), cur.isprime()))


def debug_plot_ti():
    '''Compare these to Figures 38.3, 38.4, and 38.5 respectively in Mazur/Stein's book on Riemann Hypothesis'''
    plot_ti(14)
    plot_ti(100)
    plot_ti(1000)


def debug_ti_no_waste():
    '''Tests ti_no_waste up to 100 and demonstrates the proper use'''
    _, ints = ti_no_waste(1)
    for r in range(0, 100):
        sum, ints = ti_no_waste(r, ints)
        print("r: {} sum: {}".format(r, sum))


def debug_gauss_to_coord_data():
    _, ints = ti_no_waste(100)
    x = []
    y = []
    for row in ints:
        cur_x, cur_y = gauss_to_coord_data(row)
        x.extend(cur_x)
        y.extend(cur_y)
    plt.scatter(x, y)
    plt.show()


def debug_plot_primes():
    plot_primes(1000)

###########
# TESTING #
###########


debug_gauss_to_coord_data()
debug_ti_no_waste()
debug_plot_ti()
debug_plot_primes()
# debug_gauss_int()
