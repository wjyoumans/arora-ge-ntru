#include <cassert>
#include <flint.h>
#include <nmod.h>
#include <nmod_poly.h>
#include <nmod_mat.h>

#include "keygen.hpp"
#include "logging.hpp"


const char *var = "x";

void NTRUKeyGen::init(int degree, int q, int coeffs, int ring, int seed, int log_level) {
  
  this->degree_ = degree;
  this->q_ = q;
  this->coeffs_ = coeffs;
  this->ring_ = ring;
  this->seed_ = seed;
  this->log_level_ = log_level;

  if (seed == -1) {
    seed = time(NULL);
  }

  flint_randinit(this->state);
  flint_randseed(state, seed, seed);

  nmod_init(&this->q_nmod_, q);
  
  nmod_poly_init_mod(this->modulus, this->q_nmod());
  nmod_poly_init_mod(this->den_, this->q_nmod());
  nmod_poly_set_coeff_ui(this->modulus, degree, 1);
  if (ring == 1) {
    nmod_poly_set_coeff_ui(this->modulus, 0, q-1); 
  } else if(ring == 2) {    
    nmod_poly_set_coeff_ui(this->modulus, 0, 1);
  } else if(ring == 3) {
    nmod_poly_set_coeff_ui(this->modulus, 0, q-1);
    nmod_poly_set_coeff_ui(this->modulus, 1, q-1);
  } else if(ring == 4) {
    assert(this->degree_ % 2 == 0);
    nmod_poly_set_coeff_ui(this->modulus, 0, 1);
    nmod_poly_set_coeff_ui(this->modulus, this->degree_/2, q-1);
  } else {
    throw std::invalid_argument("ring must be 1 (NTRU: x^n - 1) or 2 (NTRU2: x^n + 1) or 3 (NTRUPrime: x^p - x - 1) or or 4 (NTTRU: x^n - x^(n/2) + 1).");
  }
}

NTRUKeyGen::NTRUKeyGen(int degree, int q, int coeffs, int ring, int seed, int log_level) {
  this->init(degree, q, coeffs, ring, seed, log_level);
}

NTRUKeyGen::NTRUKeyGen(int degree, int q, int coeffs, int ring) {
  this->init(degree, q, coeffs, ring, 0, 0);
}

NTRUKeyGen::~NTRUKeyGen() {
  nmod_poly_clear(this->den_);
  nmod_poly_clear(this->modulus);
  flint_randclear(this->state);
}

void NTRUKeyGen::rand_poly(nmod_poly_t f) {
  int n = this->degree();
  int c = this->coeffs();
  flint_rand_t &state = this->state;
  
  nmod_poly_fit_length(f, n);

  do {
    // c == 3 then ternary
    if (c == 3) {
      int q = nmod_poly_modulus(f);
      for (int i = 0; i < n; i++) {
        if (n_randint(state, 2)) {      
          nmod_poly_set_coeff_ui(f, i, q - n_randint(state, 2));
        } else {
          nmod_poly_set_coeff_ui(f, i, n_randint(state, 2));
        }
      }
    // c == 2 then binary
    } else if (c == 2) {
      for (int i = 0; i < n; i++) {
        nmod_poly_set_coeff_ui(f, i, n_randint(state, 2));
      }
    } else {
      throw std::invalid_argument("coeffs must be 2 or 3.");
    }
  } while ( nmod_poly_is_zero(f) );
}

bool NTRUKeyGen::generate(nmod_mat_t H_mat, int nkeys) {
  set_log_level(this->log_level());

  int n = this->degree();
  nmod_t q_nmod = this->q_nmod();

  nmod_poly_t temp;
  nmod_poly_init_mod(temp, q_nmod);

  // find invertible denominator g  
  this->rand_poly(this->den_);
  nmod_poly_gcd(temp, this->den_, this->modulus);
  
  while (!nmod_poly_is_one(temp)) {
    this->rand_poly(this->den_);
    nmod_poly_gcd(temp, this->den_, this->modulus);
  }
  //std::cout << "# g = ";
  //nmod_poly_print_pretty(this->den_, var);
  //std::cout << '\n';
  nmod_poly_t g_inv;
  nmod_poly_init_mod(g_inv, q_nmod);
  nmod_poly_invmod(g_inv, this->den_, this->modulus);
  //nmod_poly_gcdinv(res, g_inv, keygen.g, keygen.modulus);

  // check g_inv
  nmod_poly_mul(temp, this->den_, g_inv);
  nmod_poly_rem(temp, temp, this->modulus);
  if (!nmod_poly_is_one(temp)) {
    debug("Something went wrong in invmod!\n");
    //nmod_poly_print_pretty(temp, var);
    //std::cout << '\n';
    return 0;
  }
  // generate private and public keys
  int j;
  for (int i = 0; i < nkeys; i++) {    
    this->rand_poly(temp);
    nmod_poly_mul(temp, temp, g_inv);
    nmod_poly_rem(temp, temp, this->modulus);
    for (j = 0; j < n; j++) {
      nmod_mat_set_entry(H_mat, i, j, nmod_poly_get_coeff_ui(temp, j));
    }
  }
  
  nmod_poly_clear(g_inv);
  nmod_poly_clear(temp);
  return 1;
}
