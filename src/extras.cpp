
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <flint.h>
#include <nmod.h>
#include <nmod_mat.h>
#include <nmod_poly.h>

#include "extras.hpp"

// convert 1xn matrix to polynomial. Assumes mat and poly have
// correct size
void nmod_poly_from_nmod_mat(nmod_poly_t poly, nmod_mat_t mat) {
  int n = nmod_mat_ncols(mat);
  for (int i = 0; i < n; i++) {
    nmod_poly_set_coeff_ui(poly, i, nmod_mat_get_entry(mat, 0, i));
  }
}

// convert polynomial to 1xn matrix. Assumes mat and poly have
// correct size
void nmod_mat_from_nmod_poly(nmod_mat_t mat, nmod_poly_t poly) {
  int n = nmod_mat_ncols(mat);
  for (int i = 0; i < n; i++) {
    nmod_mat_set_entry(mat, 0, i, nmod_poly_get_coeff_ui(poly, i));
  }
}

int power_product(nmod_mat_t row, std::vector<int> const &idx_vec) {
  int res = 1;
  int x;
  nmod_t ctx = row->mod;
  
  for (size_t i = 0; i < idx_vec.size(); i++) {
    x = nmod_mat_get_entry(row, 0, idx_vec[i]);
    res = nmod_mul(res, x, ctx);
  }
  return res;
}

// Only for small values.
int factorial(int n, int k) {
  int res = 1;
  for (int i = k + 1; i <= n; i++) {
    res *= i;
  }
  return res;
}

// Uses mpz_bin_uiui, valgrind complains memory possibly lost with
// fmpz_bin_uiui
ulong bin_uiui(ulong n, ulong k) {
  mpz_t bin_mpz;
  mpz_init(bin_mpz);

  if (n >= (UWORD(1) << 32)) {
    throw std::invalid_argument("Value n too large.");
  }

  if (k >= (UWORD(1) << 32)) {
    throw std::invalid_argument("Value k too large.");
  }

  mpz_bin_uiui(bin_mpz, n, k);
  ulong res = mpz_get_ui(bin_mpz);

  mpz_clear(bin_mpz);
  return res;
}

std::vector<std::string> parse_line(std::string line) {
  assert(line[line.length() - 1] == ']');

  while (line[line.length() - 1] == ']' || line[line.length() - 1] == ' ') {
    line.pop_back();
  }

  std::istringstream iss(line);
  std::string f;
  std::vector<std::string> v;

  while (!iss.eof()) {
    iss >> f;
    v.push_back(f);
  }
  return v;
}

void nmod_mat_init_from_stream(nmod_mat_t mat, int q, std::istream& is) {
    char c;
    std::string line;
    int nrows, ncols;

    std::vector<std::vector<std::string>> data;

    is >> c;
    if (c != '[') {
        throw "Invalid matrix";
        //std::cout << "error? 1\n";
    }
    while (!is.eof()) {
        is >> c;
        if (is.eof()) {
            break;
        }
        if (c == '[') {
            std::getline(is, line);
            auto row = parse_line(line);
            data.push_back(row);
        } else if (c == ']') {
            break;
        } else {
            throw "Invalid matrix";
            //std::cout << "error? 1\n";
            //break;
        }
    }

    nrows = data.size();
    ncols = data[0].size();
    nmod_mat_init(mat, nrows, ncols, q);
    //assert(nrows == nmod_mat_nrows(mat));
    //assert(ncols == nmod_mat_ncols(mat));

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            nmod_mat_set_entry(mat, i, j, stoi(data[i][j]));
        }
    }
}

void nmod_mat_to_stream(nmod_mat_t mat, std::ostream& os) {
    int nrows = nmod_mat_nrows(mat);
    int ncols = nmod_mat_ncols(mat);

    int i, j;
    os << "[";
    for (i = 0; i < nrows; i++) {
        os << "[";
        for (j = 0; j < ncols; j++) {
            os << std::to_string(nmod_mat_get_entry(mat, i, j));
            if (j < ncols - 1) {
                os << " ";
            } else {
                os << "]" << '\n';
            }
        }
    }
    os << "]";
}