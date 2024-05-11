
#pragma once

#include <vector>

#include <flint.h>
#include <nmod_mat.h>

#include "keygen.hpp"

ulong num_variables(int n, int mode);

std::vector<ulong> binomials(ulong n, ulong k);

void arora_ge_system(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& keygen);

void arora_ge_system_ntru_new(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& keygen);

void arora_ge_system_ntru(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& keygen);

void arora_ge_system_ntru2(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& keygen);

void arora_ge_system_generic(nmod_mat_t res, nmod_mat_t H_mat, const NTRUKeyGen& keygen);
