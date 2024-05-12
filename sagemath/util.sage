# Miscellaneous functions

from itertools import product

##### NTRU keygen stuff ##### 

def ntru_ring(n, q):
    """ Return the ring Z_q[x]/(x^n - 1)."""
    R.<z> = GF(q)[]
    RR.<x> = R.quotient(z^n - 1)
    return RR

def random_element(R, lower=-1, upper=1):
    """ 
    Random element of R with coefficients bounded by lower and upper, 
    inclusive.
    """
    n = R.degree()
    return R([randint(lower, upper) for i in range(n)])

def ntru_keygen(R, nkeys, f_coeffs=[0,1], g_coeffs=[], g=None):
    """ 
    Generate `nkeys` NTRU instances with specified support for f and g. 
    By default g has the same support as f. Optionally pass g.
    """
    # we'll assume g_coeffs = f_coeffs as the default
    if g_coeffs == []:
        g_coeffs = f_coeffs

    n = R.degree()
    if g is None:
        while True:
            g = random_element(R, g_coeffs[0], g_coeffs[-1])
            if len(g.lift().coefficients()) > 1 and g.is_unit():
                break
        
    ginv = 1/g
    
    F = []
    H = []
    while len(H) < nkeys:
        f = random_element(R, f_coeffs[0], f_coeffs[-1])
        if len(f.lift().coefficients()) < 2 or f in F:
            continue

        # NOTE: Julia code checked f coprime with all F[i]. why?

        h = f*ginv
        F.append(f)
        H.append(h)
        
    return H, F, g

##### Miscellaneous #####

def multiplication_matrix(f):
    """ Multiplication/rotation matrix of f."""
    R = f.parent()
    B = R.base_ring()
    n = R.degree()
    z = R.gen()

    temp = f
    M = zero_matrix(B, n, n)
    for i in range(n):
        M[:,i] = vector(temp.list())
        temp *= z

    return M.transpose()
