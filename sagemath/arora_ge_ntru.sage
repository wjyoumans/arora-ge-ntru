# Load all files, and define script behavior.

# attach causes an error when called as a script so we use load
load("util.sage")
load("system.sage")
load("kernel.sage")
load("solve.sage")

if __name__ == '__main__' and '__file__' in globals():
    import argparse
    parser = argparse.ArgumentParser("arora_ge_ntru.sage")
    parser.add_argument("n", help="degree of the underlying ring")
    parser.add_argument("q", help="prime modulus")
    parser.add_argument("m", help="number of keys/samples to generate")
    parser.add_argument("-c", "--coeffs",
        help="number of coefficients i.e. 2 for binary, 3 for ternary",
        default=2
    )
    parser.add_argument("-r", "--ring",
        help="set to 1 for NTRU: x^n - 1, 2 for NTRU2: x^n + 1",
        default=1
    )
    args = parser.parse_args()

    n = int(args.n)
    q = int(args.q)
    m = int(args.m)
    c = int(args.coeffs)
    r = int(args.ring)

    if c == 2:
        supp_f = [0,1]
    elif c == 3:
        supp_f = [-1,0,1]
    else:
        print("c must be 2 or 3.")
        exit()
    
    R = ntru_ring(n, q)
    H, F, g = ntru_keygen(R, m, supp_f)
    b, gg = solve_arora_ge(H, supp_f, True, False)
    
    if b:
        print(f"Recovered gg with gg = g*{gg/g}.")
    else:
        print("FAILURE.")

    
