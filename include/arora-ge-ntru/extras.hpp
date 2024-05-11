
#pragma once

#include <vector>
#include <iostream>
#include <nmod_mat.h>

void nmod_mat_from_nmod_poly(nmod_mat_t mat, nmod_poly_t poly);
void nmod_poly_from_nmod_mat(nmod_poly_t poly, nmod_mat_t mat);

int power_product(nmod_mat_t row, std::vector<int> const &idx_vec);

int factorial(int n, int k);

ulong bin_uiui(ulong n, ulong k);

void nmod_mat_init_from_stream(nmod_mat_t mat, int q, std::istream& is);

void nmod_mat_to_stream(nmod_mat_t mat, std::ostream& os);

