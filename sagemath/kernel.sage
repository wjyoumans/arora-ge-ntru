# Kernel reduction algorithm

def reduce_kernel(ker, n, supp_f=[0,1]):
    """
    Reduce kernel dimension by computing the kernels of submatrices.

    Returns True if the reduced kernel has rank 1 or False otherwise, as well
    as the last n columns of the reduced kernel matrix.

    If the rank is 1, then this is a vector of the coefficients of a rotation
    of g where the highest order terms come first.
    """
    
    if ker.nrows() > n:
        print(f"FAILURE: Linear system underdetermined (too few keys).")
        return False, ker[:,-n:]
    
    d = len(supp_f)
    runs = [binomial(n+d-1-k, d-1) for k in range(1,n+1)]        
    assert(sum(runs) == ker.ncols()-n)
    
    offset = 0
    zero_coeffs = ker[:,:runs[0]]
    offset += runs[0]
    m = ker.nrows()

    print(f"Computing kernel of {m} x {runs[0]} submatrix...", end="\n")
    k = zero_coeffs.kernel().matrix()
    r = k.nrows()
    print(f"Initial kernel rank: {r}.\n")

    if r == 0:
        print(f"FAILURE: Denominator has hamming weight n = {n}.")
        return False, ker[:,-n:]
    elif r > n:
        print(f"FAILURE: Linear system underdetermined (too few keys).")
        return False, ker[:,-n:]
    elif r == 1:
        print(f"SUCCESS: Denominator has hamming weight n-1 = {n-1}.")
        temp = k*ker
        return True, temp[:,-n:]

    for i in range(1, len(runs)):
        temp = zero_coeffs.augment(ker[:,offset+1:offset+runs[i]])

        print(f"Computing kernel of {m} x {temp.ncols()} submatrix...", end="\n")
        k = temp.kernel().matrix()
        r = k.nrows()

        if r > 0:
          print(f"Reduced kernel rank: {r}.\n")
          zero_coeffs = temp
          if r == 1:
            break
        else:
            print("Kernel rank unchanged (Dim(S) = 0).\n")
            
        offset += runs[i]

    reduced = k*ker
    
    if r != 1:
        print("FAILURE: Need more samples.")
        return False, reduced 

    print("SUCCESS.")
    return True, reduced[:,-n:]
