
#pragma once

#include <nmod_mat.h>
#include "keygen.hpp"

int arora_ge_recover(nmod_mat_t den, nmod_mat_t system, NTRUKeyGen& ctx);

void arora_ge_recover_nullonly(nmod_mat_t ker, nmod_mat_t system);
