# Solve NTRU or NTRU with multiple keys.

from timeit import default_timer as timer

def solve_arora_ge(H, supp_f, rotations=True, early_abort=True):
    """
    Try to recover the secret g from multiple instances h_i = f_i/g, where
    f_i have support supp_f.
    """
    R = H[0].parent()
    n = R.degree()

    print("Building linear system... ", end="", flush=True)
    t = timer()
    M = arora_ge_system(H, supp_f, rotations, early_abort)
    t = timer() - t
    print(f"{t} seconds.")
    print(f"System dimensions: {M.nrows()} x {M.ncols()}.\n")
        
    print("Computing right kernel... ", end="", flush=True)
    t = timer()
    ker = M.right_kernel_matrix()
    rk = ker.nrows()
    t = timer() - t
    print(f"{t} seconds.")
    print(f"Kernel dimensions: {ker.nrows()} x {ker.ncols()}\n")
    
    print("Doing kernel reduction...", flush=True)
    t = timer()
    b, gg = reduce_kernel(ker, n, supp_f)
    t = timer() - t
    print(f"Kernel reduction time: {t} seconds.\n")

    if b:
        return True, R([x for x in gg[0]][::-1])

    return False, R(0)
