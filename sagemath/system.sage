# Generate Arora-Ge linear system.

def arora_ge_system(
    H, 
    supp_f, 
    rotations=True,
    early_abort=True):
    """ 
    Linearize the polynomial system given by 
    prod_{c in supp_f} (Constant(h_i*X) - c) for h_i in H, and return the
    corresponding matrix. Avoids manipulating any actual monomials.

    * rotations -- By default add all "rotations" of elements of H. Generally
      only meaningful in the rings Z[x]/(x^n +/- 1).

    * early_abort -- Only process enough samples such that the linearized
      system has size (n + epsilon) x n.
    """
    
    R = H[0].parent()
    n = R.degree()
    Zq = R.base_ring()
    q = Zq.order()

    monomials = all_monomials(n, supp_f)

    # defaults
    nkeys = len(H)
    nkeys_used = nkeys
    nrows = nkeys
    ncols = len(monomials)

    # if early_abort then we use the minimum number of keys such that 
    # nrows > ncols, taking into account rotations if used. 
    if early_abort and rotations:
        nkeys_used = min(nkeys, ncols//n + 1)
        nrows = nkeys_used*n
    elif early_abort:
        nkeys_used = min(nkeys, ncols + 1)
        nrows = nkeys_used
    elif rotations:
        nrows *= n

    M = zero_matrix(Zq, nrows, ncols)
    for i in range(nkeys_used):
        h = H[i]
        
        if rotations:
            A = multiplication_matrix(h)
        else:
            A = matrix(h.list())

        # m = 1 if no rotations, n if rotations
        m = A.nrows()
        for j in range(m):
            row = []
            for mon, c in monomials:
                res = 1
                for (s, e) in enumerate(mon):
                    if e != 0:
                        res *= A[j,s]^e % q

                row.append(c*res % q)
            M[m*i + j,:] = vector(row)
        
    return M


def all_monomials(n, supp_f):
    zx.<x> = PolynomialRing(ZZ, 'x')
    F = prod([x - c for c in supp_f])
    d = len(supp_f)
        
    # Compute all monomials and their coefficients. Separate monomials
    # of degree d > 1 and d == 1 so we can order them properly.
    monomials1 = []
    monomials2 = []
    for d, c1 in enumerate(F.list()):
        # If c != 0  then system contains monomials of degree d.
        if c1 != 0:
            data = monomial_data(n, d)
            temp = [(mon, c1*c2) for mon, c2 in data]
            if d != 1:
                monomials1 += temp
            else:
                monomials2 += temp


    # put degree > 1 monomials in grlex order
    monomials1.sort(key=lambda x: x[0])
    monomials1.reverse()

    # put degree 1 monomials in grlex order
    monomials2.sort(key=lambda x: x[0])
    monomials2.reverse()

    return monomials1 + monomials2


def monomial_data(n, d):
    """ 
    Return exponent vectors of all monomials of degree d in n variables
    as well as all monomial coefficients. 
    """
    monomials = WeightedIntegerVectors(d, [1]*n)

    multinomial_coeffs = [falling_factorial(d, d-i) for i in range(1,d+1)]
    monomials = [(mon, multinomial_coeffs[max(mon)-1]) for mon in monomials]

    return monomials
    
