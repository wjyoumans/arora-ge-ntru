
#pragma once

#include <vector>

#include <flint.h>
#include <nmod.h>
#include <nmod_poly.h>
#include <nmod_mat.h>


// coeffs == 3 then ternary, 2 for binary
// ring = 1 for x^n - 1, ring = 2 for x^n + 1
class NTRUKeyGen {
  int degree_;
  int q_;
  int coeffs_;
  int ring_;
  int seed_;
  int log_level_;

  nmod_t q_nmod_;  
  nmod_poly_t den_;

  void init(int degree, int q, int coeffs, int ring, int seed, int log_level);
  
  public:
    NTRUKeyGen(int degree, int q, int coeffs, int ring, int seed, int log_level);
    NTRUKeyGen(int degree, int q, int coeffs, int ring);
    ~NTRUKeyGen();

    // ring modulus
    nmod_poly_t modulus;
    flint_rand_t state;

    // Accessors
    int degree() const { return degree_; }
    int q() const { return q_; }
    int coeffs() const { return coeffs_; }
    int ring() const {return ring_; }
    int seed() const { return seed_; }
    int log_level() const { return log_level_; }
    nmod_t q_nmod() const { return q_nmod_; }

    void denominator(nmod_poly_t g) { nmod_poly_set(g, this->den_); }
    //void modulus(nmod_poly_t mod) { nmod_poly_set(mod, this->mod_); }
  
    void rand_poly(nmod_poly_t f);
    bool generate(nmod_mat_t H_mat, int nkeys);
};