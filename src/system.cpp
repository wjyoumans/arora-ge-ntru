
#include <cassert>
#include <iostream>
#include <algorithm>
#include <sstream>

#include <flint.h>
#include <nmod.h>
#include <nmod_poly.h>
#include <nmod_mat.h>

#include "keygen.hpp"
#include "logging.hpp"
#include "system.hpp"
#include "extras.hpp"

// Binomial coefficients n + k - i choose k for i = 1 to n.
std::vector<ulong> binomials(ulong n, ulong k) {
  std::vector<ulong> res(n);
  ulong temp = n + k;
  for (uint i = 1; i <= n; i++) {
    res[i-1] = bin_uiui(temp - i, k);
  }
  return res;
}

ulong num_variables(int n, int mode) {
  ulong res = bin_uiui((ulong)(n + mode - 1), (ulong)mode) + n;
  return res;
}

// Set mat to multiplication matrix of h in Z_q[x]/(mod)
void multiplication_matrix(nmod_mat_t mat, nmod_poly_t h, const NTRUKeyGen& keygen) {
  int n = keygen.degree();
  assert(nmod_mat_ncols(mat) == n);
  assert(nmod_mat_nrows(mat) == n);

  nmod_poly_t x;
  nmod_poly_init_mod(x, h->mod);
  nmod_poly_set_coeff_ui(x, 1, 1);
  
  ulong c;
  for (int i = 0; i < n; i++) {
    for (int j =0; j < n; j++) {
      c = nmod_poly_get_coeff_ui(h, j);
      nmod_mat_set_entry(mat, j, i, c);
    }
    nmod_poly_mul(h, h, x);
    nmod_poly_rem(h, h, keygen.modulus);
  }
  nmod_poly_clear(x);
}

void arora_ge_system(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& keygen) {
  set_log_level(keygen.log_level());

  int r = keygen.ring();
  if (r == 1) {
    arora_ge_system_ntru(res, H_mat, keygen);
  } else if (r == 2) {
    arora_ge_system_ntru2(res, H_mat, keygen);
  } else {
    arora_ge_system_generic(res, H_mat, keygen);
  }

}

int index(std::vector<int> comb) {
  int d = comb.size();
  int u = 0;
  for (int i = 1; i < d; i++) {
    if (comb[i-1] == comb[i]) {
      u++;
    } 
  }
  return u;
}

void arora_ge_system_ntru(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& ctx) {
  int n = ctx.degree();
  int d = ctx.coeffs();
  int q = ctx.q();
  int nkeys = nmod_mat_nrows(H_mat);
  nmod_t q_nmod = ctx.q_nmod();
  
  std::vector<int> monomial_coeffs(d);
  for (int i = 1; i <= d; i++) {
    monomial_coeffs[i-1] = factorial(d, i);    
  }
  
  int m = n + d - 1;
  std::vector<bool> perm(m);
  std::fill(perm.begin(), perm.begin() + d, true);

  // Iterate over unique monomials.
  int i, j, k, x, idx, n_unique;
  std::vector<int> comb(d);
  std::vector<std::vector<int>> combs;
  std::vector<int> unique;

  // precompute combinations:
  do {
    idx = 0;
    n_unique = 0;
    for (i = 0; i < m; ++i) {
      if (perm[i]) {
        comb[idx] = i - idx;
        if (idx > 0 && comb[idx-1] == comb[idx]) {
          n_unique++;
        }
        idx++;
      }
    }
    combs.push_back(comb);
    unique.push_back(n_unique);
  } while (std::prev_permutation(perm.begin(), perm.end()));

  nmod_mat_t mult, window;  
  nmod_mat_init(mult, n, n, q);
  for (i = 0; i < nkeys; i++) {
    // make multiplication matrix
    for (j = 0; j < n; j++) {
      x = nmod_mat_get_entry(H_mat, i, j);
      for (k = 0; k < n; k++) {
        nmod_mat_set_entry(mult, k, (k-j+n) % n, x);
      }
    }

    // set first block to negative of multiplication matrix
    nmod_mat_window_init(window, res, n*i, 0, n*i+n, n);
    nmod_mat_neg(window, mult);
    nmod_mat_window_clear(window);

    for (j = 0; j < n; j++) {
      nmod_mat_window_init(window, mult, j, 0, j+1, n);
      for (size_t k = 0; k < combs.size(); k++) {
        comb = combs[k];
        n_unique = unique[k];        
      
        x = power_product(window, comb);
        x = nmod_mul(monomial_coeffs[n_unique], x, q_nmod);
        nmod_mat_set_entry(res, n*i + j, n + k, x);
      }
      nmod_mat_window_clear(window);
    }
  }
  nmod_mat_clear(mult);
}

void arora_ge_system_ntru2(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& ctx) {
  int n = ctx.degree();
  int d = ctx.coeffs();
  int q = ctx.q();
  int nkeys = nmod_mat_nrows(H_mat);
  nmod_t q_nmod = ctx.q_nmod();
  
  std::vector<int> monomial_coeffs(d);
  for (int i = 1; i <= d; i++) {
    monomial_coeffs[i-1] = factorial(d, i);    
  }
  
  int m = n + d - 1;
  std::vector<bool> perm(m);
  std::fill(perm.begin(), perm.begin() + d, true);

  // Iterate over unique monomials.
  int i, j, k, x, y, idx, n_unique;
  std::vector<int> comb(d);
  std::vector<std::vector<int>> combs;
  std::vector<int> unique;

  // precompute combinations:
  do {
    idx = 0;
    n_unique = 0;
    for (i = 0; i < m; ++i) {
      if (perm[i]) {
        comb[idx] = i - idx;
        if (idx > 0 && comb[idx-1] == comb[idx]) {
          n_unique++;
        }
        idx++;
      }
    }
    combs.push_back(comb);
    unique.push_back(n_unique);
  } while (std::prev_permutation(perm.begin(), perm.end()));

  nmod_mat_t mult, window;  
  nmod_mat_init(mult, n, n, q);
  for (i = 0; i < nkeys; i++) {
    // make multiplication matrix
    for (j = 0; j < n; j++) {
      x = nmod_mat_get_entry(H_mat, i, j);
      for (k = 0; k < j; k++) {
        y = nmod_neg(x, q_nmod);
        nmod_mat_set_entry(mult, k, (k-j+n) % n, y);
      }      
      for (k = j; k < n; k++) {
        nmod_mat_set_entry(mult, k, (k-j+n) % n, x);
      }
    }

    // set first block to negative of multiplication matrix
    nmod_mat_window_init(window, res, n*i, 0, n*i+n, n);
    nmod_mat_neg(window, mult);
    nmod_mat_window_clear(window);

    for (j = 0; j < n; j++) {
      nmod_mat_window_init(window, mult, j, 0, j+1, n);
      for (size_t k = 0; k < combs.size(); k++) {
        comb = combs[k];
        n_unique = unique[k];        
      
        x = power_product(window, comb);
        x = nmod_mul(monomial_coeffs[n_unique], x, q_nmod);
        nmod_mat_set_entry(res, n*i + j, n + k, x);
      }
      nmod_mat_window_clear(window);
    }
  }
  nmod_mat_clear(mult);
}

void arora_ge_system_generic(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& ctx) {
  //throw "arora_ge_system_generic not yet implemented";
  int n = ctx.degree();
  int d = ctx.coeffs();
  int q = ctx.q();
  int nkeys = nmod_mat_nrows(H_mat);
  nmod_t q_nmod = ctx.q_nmod();
  std::vector<int> monomial_coeffs(d);
  for (int i = 1; i <= d; i++) {
    monomial_coeffs[i-1] = factorial(d, i);
  }
  
  int m = n + d - 1;
  std::vector<bool> perm(m);
  std::fill(perm.begin(), perm.begin() + d, true);

  // Iterate over unique monomials.
  int i, j, k, x, idx, n_unique;
  std::vector<int> comb(d);
  std::vector<std::vector<int>> combs;
  std::vector<int> unique;

  // precompute combinations:
  do {
    idx = 0;
    n_unique = 0;
    for (i = 0; i < m; ++i) {
      if (perm[i]) {
        comb[idx] = i - idx;
        if (idx > 0 && comb[idx-1] == comb[idx]) {
          n_unique++;
        }
        idx++;
      }
    }
    combs.push_back(comb);
    unique.push_back(n_unique);
  } while (std::prev_permutation(perm.begin(), perm.end()));

  nmod_mat_t mult, window;
  nmod_mat_init(mult, n, n, q);
  nmod_poly_t hi, hij, xj, x1;
  nmod_poly_init_mod(hi, q_nmod);
  nmod_poly_init_mod(hij, q_nmod);
  nmod_poly_init_mod(xj, q_nmod);
  nmod_poly_init_mod(x1, q_nmod);
  // for each i
  for (i = 0; i < nkeys; i++) {
    // make multiplication matrix for each h_i
    for (j = 0; j < n; j++) {
      x = nmod_mat_get_entry(H_mat, i, j);
      nmod_poly_set_coeff_ui(hi, j, x);
    }
    nmod_poly_set_coeff_ui(x1, 1, 1);
    nmod_poly_zero (xj);
    nmod_poly_set_coeff_ui(xj, 0, 1);
    for (j = 0; j < n; j++) {
      nmod_poly_mul(hij, hi, xj);
      nmod_poly_rem(hij, hij, ctx.modulus);
      for (k = 0; k < n; k++) {
        x = nmod_poly_get_coeff_ui (hij, k);
        nmod_mat_set_entry (mult, k, j, x);
      }
      nmod_poly_mul(xj, xj, x1);
    }
    // set first block to negative of multiplication matrix
    nmod_mat_window_init(window, res, n*i, 0, n*i+n, n);
    nmod_mat_neg(window, mult);
    nmod_mat_window_clear(window);
    for (j = 0; j < n; j++) {
      nmod_mat_window_init(window, mult, j, 0, j+1, n);
      for (size_t k = 0; k < combs.size(); k++) {
        comb = combs[k];
        n_unique = unique[k];
        x = power_product(window, comb);
        x = nmod_mul(monomial_coeffs[n_unique], x, q_nmod);
        nmod_mat_set_entry(res, n*i + j, n + k, x);
      }
      nmod_mat_window_clear(window);
    }
  }
  nmod_mat_clear(mult);
  nmod_poly_clear(hi);
  nmod_poly_clear(hij);
  nmod_poly_clear(xj);
  nmod_poly_clear(x1);
}
