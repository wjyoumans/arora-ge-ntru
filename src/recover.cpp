
#include <iostream>
#include <fstream>
#include <chrono>
#include <unistd.h>

#include <flint.h>
#include <nmod.h>
#include <nmod_poly.h>
#include <nmod_mat.h>

#include "system.hpp"
#include "keygen.hpp"
#include "logging.hpp"

using namespace std;
using namespace std::chrono;


int arora_ge_recover(nmod_mat_t den, nmod_mat_t system, NTRUKeyGen& ctx) {
  set_log_level(ctx.log_level());

  //int n = 31;
  int n = ctx.degree();
  int q = ctx.q();
  int d = ctx.coeffs();
  int ncols = nmod_mat_ncols(system);
  int status = 0;

  std::vector<ulong> bins = binomials(n, d-1);
    
  nmod_mat_t initial_kernel, window;
  nmod_mat_init(initial_kernel, ncols, ncols, q);
  int rank = nmod_mat_nullspace(initial_kernel, system);
      
  debug("Initial kernel rank: ", rank, "\n");

  bool terminate = false;
  if (rank == 1) {
    debug("SUCCESS: Kernel rank is 1.\n");
    terminate = true;
    nmod_mat_window_init(window, initial_kernel, 0, 0, n, 1);
    nmod_mat_transpose(den, window);
    nmod_mat_window_clear(window);
  } else if (rank > n) {
    debug("FAILURE: Linear system underdetermined!\n");
    terminate = true;
    status = 1;
  } else if (rank != n) {
    debug("Something went wrong. Kernel rank is ", rank, ".\n");
    terminate = true;
    status = 1;
  }
  
  if (terminate) {
    nmod_mat_clear(initial_kernel);
    return status;
  }
  
  // kernel rank is n.
  nmod_mat_t kernel, temp;
  nmod_mat_window_init(window, initial_kernel, 0, 0, ncols, n);
  nmod_mat_init_set(kernel, window);
  nmod_mat_window_clear(window);
  nmod_mat_clear(initial_kernel);

  int offset = n;
  nmod_mat_t res;
  nmod_mat_window_init(window, kernel, offset, 0, offset+bins[0], n);
  nmod_mat_init(res, n, ncols, q);
  
  rank = nmod_mat_nullspace(res, window);
  debug("New kernel rank: ", rank, "\n");
  int hw = n - rank;
  debug("Denominator hamming weight: ", hw, "\n");

  if (rank == 0) {
    debug("FAILURE: Denominator has hamming weight n.\n");
    terminate = true;
    status = 1;
  } else if (rank == 1) {
    debug("SUCCESS: Denominator has hamming weight n - 1.\n");
    terminate = true;

    nmod_mat_window_clear(window); 
    nmod_mat_window_init(window, res, 0, 0, n, 1);
  
    nmod_mat_init(temp, ncols, 1, q);
    nmod_mat_mul(temp, kernel, window);
    
    nmod_mat_window_clear(window);
    nmod_mat_window_init(window, temp, 0, 0, n, 1);
    nmod_mat_transpose(den, window);
    nmod_mat_clear(temp);
  }
  
  if (terminate) {
    nmod_mat_clear(res);
    nmod_mat_clear(kernel);
    nmod_mat_window_clear(window);
    return status;
  }

  nmod_mat_t submat;
  nmod_mat_init_set(submat, window);
  nmod_mat_init(temp, 0, 0, q);

  offset += bins[0];
  for (int i = 1; i < n; i++) {
    nmod_mat_window_clear(window);
    nmod_mat_window_init(window, kernel, offset, 0, offset+bins[i], n);
    
    nmod_mat_clear(temp);
    //nmod_mat_init(temp, offset+bins[i]-n, n, q);
    nmod_mat_init(temp, nmod_mat_nrows(submat) + nmod_mat_nrows(window), n, q);
    nmod_mat_concat_vertical(temp, submat, window);

    rank = nmod_mat_nullspace(res, temp);
    debug("New kernel rank: ", rank, "\n");

    if (rank > 0) {
      nmod_mat_clear(submat);
      nmod_mat_init_set(submat, temp);
      if (rank == 1) {
        break;
      }
    }
    offset += bins[i];
  }

  if (rank != 1) {
    debug("FAILURE: Reason unknown.\n");
    status = 1;
  }

  nmod_mat_window_clear(window); 
  nmod_mat_window_init(window, res, 0, 0, n, 1);
  
  nmod_mat_clear(temp);
  nmod_mat_init(temp, ncols, 1, q);
  nmod_mat_mul(temp, kernel, window);

  nmod_mat_window_clear(window);
  nmod_mat_window_init(window, temp, 0, 0, n, 1);
  nmod_mat_transpose(den, window);

  nmod_mat_clear(res);
  nmod_mat_clear(kernel);
  nmod_mat_clear(temp);
  nmod_mat_clear(submat);
  nmod_mat_window_clear(window);
  
  return status;
}


int arora_ge_recover_nullonly(nmod_mat_t ker, nmod_mat_t system) {
  auto t0 = high_resolution_clock::now();  
  nmod_mat_nullspace(ker, system);
  auto t1 = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(t1-t0);
  //nmod_mat_print(ker);
  
  int tSize = 0, resident = 0, share = 0;
  ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
  double rss = resident * page_size_kb;

  std::cout << "# mem: " << rss/1000.0 << " M," << " time: " << duration.count()/1000000.0 << " s" << endl;

  return 0;
}