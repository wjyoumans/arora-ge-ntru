#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

#include <flint.h>
#include <nmod_mat.h>
#include <nmod_poly.h>
#include <system_error>
#include <chrono>
#include <unistd.h>

#include "argparse/argparse.hpp"

#include "arora-ge-ntru/keygen.hpp"
#include "arora-ge-ntru/system.hpp"
#include "arora-ge-ntru/recover.hpp"
#include "arora-ge-ntru/extras.hpp"
#include "arora-ge-ntru/logging.hpp"

using namespace std::chrono;

void keygen(argparse::ArgumentParser& program, NTRUKeyGen& ctx) {
  int n = ctx.degree();
  int q = ctx.q();
  
  int nkeys = program.get<int>("--keys");

  std::string out_fn = std::string();
  if (auto out = program.present("--pk_output")) {
    out_fn = *out;
  }

  std::string sk_out_fn = out_fn + ".sk";
  if (auto out = program.present("--sk_output")) {
    sk_out_fn = *out;
  }

  debug("Generating ", nkeys, " keys.\n");
  
  // generate keys
  nmod_mat_t H_mat;
  nmod_mat_init(H_mat, nkeys, n, q);
  ctx.generate(H_mat, nkeys);

  // Save keys to file if output file specified, otherwise print to stdout
  nmod_poly_t den_poly;
  nmod_poly_init(den_poly, q);
  ctx.denominator(den_poly);

  nmod_mat_t den;
  nmod_mat_init(den, 1, n, q);
  nmod_mat_from_nmod_poly(den, den_poly);
  
  if (out_fn.empty()) {
    std::stringstream ss;
    
    nmod_mat_to_stream(H_mat, ss);
    ss << '\n';
    nmod_mat_to_stream(den, ss);
    std::cout << ss.str() << '\n';
  } else {
    std::ofstream file1, file2;
    file1.open(out_fn);
    nmod_mat_to_stream(H_mat, file1);
    
    file2.open(sk_out_fn);
    nmod_mat_to_stream(den, file2);
  }

  nmod_poly_clear(den_poly);
  nmod_mat_clear(den);
  nmod_mat_clear(H_mat);
}

void system(argparse::ArgumentParser& program, NTRUKeyGen& ctx) {
  int n = ctx.degree();
  int q = ctx.q();
  int c = ctx.coeffs();

  std::string in_fn = program.get("--input");
  std::string out_fn = std::string();
  if (auto out = program.present("-o")) {
    out_fn = *out;
  }

  debug("Reading key file.\n");
  nmod_mat_t H_mat;
  std::ifstream file;
  file.open(in_fn);
  nmod_mat_init_from_stream(H_mat, q, file);
  
  int nkeys = nmod_mat_nrows(H_mat);
  ulong nvars = num_variables(n, c);

  int nrows = n*nkeys;
  debug("Building ", nrows, " x ", nvars,  " system.\n");
  nmod_mat_t res;
  nmod_mat_init(res, nrows, nvars, q);
  arora_ge_system(res, H_mat, ctx);
  nmod_mat_clear(H_mat);
  
  // Save system to file if output file specified, otherwise print to stdout
  if (out_fn.empty()) {
    std::stringstream ss;
    nmod_mat_to_stream(res, ss);
    std::cout << ss.str() << '\n';
  } else {
    std::ofstream file;
    file.open(out_fn);
    nmod_mat_to_stream(res, file);
  }
  nmod_mat_clear(res);
}

void recover(argparse::ArgumentParser& program, NTRUKeyGen& ctx) {
  int n = ctx.degree();
  int q = ctx.q();

  std::string in_fn = program.get("--input");
  std::string out_fn = std::string();
  if (auto out = program.present("-o")) {
    out_fn = *out;
  }

  debug("Reading linear system file.\n");
  nmod_mat_t system;
  std::ifstream file;
  file.open(in_fn);
  nmod_mat_init_from_stream(system, q, file);

  if (program["--nullonly"] == true && !out_fn.empty()) {
    int ncols = nmod_mat_ncols(system);
    nmod_mat_t ker;
    nmod_mat_init(ker, ncols, n, q);
    
    debug("Computing nullspace only.\n");
    arora_ge_recover_nullonly(ker, system);
    std::ofstream file;
    file.open(out_fn);
    nmod_mat_to_stream(ker, file);
    nmod_mat_clear(ker);
  }
  else {
    nmod_mat_t den;
    nmod_mat_init(den, 1, n, q);
    
    debug("Attempting full key recovery.\n");
    int ret = arora_ge_recover(den, system, ctx);
    if (ret == 0) {
      debug("Saving key.\n");
      std::ofstream file;
      file.open(out_fn);
      nmod_mat_to_stream(den, file);
    }
    nmod_mat_clear(den);
  }

  nmod_mat_clear(system);
}

void verify(argparse::ArgumentParser& program, NTRUKeyGen& ctx) {
  int q = ctx.q();
  debug("Verifying result.\n");
  
  std::string in_fn1 = program.get("--sk_input1");
  std::string in_fn2 = program.get("--sk_input2");

  nmod_mat_t sk1, sk2;
  nmod_poly_t sk1_poly, sk2_poly;
  std::ifstream file1, file2;
  
  file1.open(in_fn1);
  nmod_mat_init_from_stream(sk1, q, file1);
  nmod_poly_init_mod(sk1_poly, ctx.q_nmod());
  nmod_poly_from_nmod_mat(sk1_poly, sk1);
  
  file2.open(in_fn2);
  nmod_mat_init_from_stream(sk2, q, file2);
  nmod_poly_init_mod(sk2_poly, ctx.q_nmod());
  nmod_poly_from_nmod_mat(sk2_poly, sk2);

  nmod_poly_t sk2_inv, temp;
  nmod_poly_init_mod(sk2_inv, ctx.q_nmod());

  nmod_poly_invmod(sk2_inv, sk2_poly, ctx.modulus);

  nmod_poly_mul(sk1_poly, sk1_poly, sk2_inv);
  nmod_poly_rem(sk1_poly, sk1_poly, ctx.modulus);

  const char* X = "x";
  
  debug("sk1/sk2 = ");
  if (log_level > 0)
    nmod_poly_print_pretty(sk1_poly, X);

  int deg = nmod_poly_degree(sk1_poly);
  nmod_poly_init_mod(temp, ctx.q_nmod());
  nmod_poly_set_coeff_ui(temp, deg, 1);

  bool success = false;
  if (nmod_poly_equal(sk1_poly, temp) == 1) {
    success = true;
  } else {
    nmod_poly_neg(temp, temp);
    if (nmod_poly_equal(sk1_poly, temp) == 1) {
      success = true;
    }
  }
  
  nmod_poly_clear(temp);
  nmod_poly_clear(sk1_poly);
  nmod_poly_clear(sk2_poly);
  nmod_poly_clear(sk2_inv);
  
  nmod_mat_clear(sk1);
  nmod_mat_clear(sk2);
  
  assert(success);
}


void all (argparse::ArgumentParser& program, NTRUKeyGen& ctx) {
  int n = ctx.degree();
  int q = ctx.q();
  int c = ctx.coeffs();
  int nkeys = program.get<int>("--keys");

  // generate keys and output poly
  nmod_mat_t H_mat, den, den_found;
  nmod_poly_t den_poly, den_found_poly;
  std::stringstream ss;
  nmod_mat_init(H_mat, nkeys, n, q);
  nmod_poly_init(den_poly, q);
  nmod_mat_init(den, 1, n, q);
  nmod_mat_init(den_found, 1, n, q);
  
  ctx.generate(H_mat, nkeys);
  ctx.denominator(den_poly);
  nmod_mat_from_nmod_poly(den, den_poly);
  nmod_mat_to_stream(H_mat, ss);
  ss << '\n';
  nmod_mat_to_stream(den, ss);
  std::cout << ss.str() << '\n';

  // build system
  ulong nvars = num_variables(n, c);
  nmod_mat_t system;
  nmod_mat_init(system, n*nkeys, nvars, q);
  arora_ge_system(system, H_mat, ctx);

  if (0) {
    nmod_mat_to_stream(system, ss);
    std::cout << ss.str() << '\n';
  }
  
  // solve linear system
  auto t0 = high_resolution_clock::now();  
  arora_ge_recover(den_found, system, ctx);
  auto t1 = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(t1-t0);  

  // copied from https://stackoverflow.com/questions/669438
  int tSize = 0, resident = 0, share = 0;
  std::ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();
  long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
  double rss = resident * page_size_kb;
  //double shared_mem = share * page_size_kb;
  nmod_mat_clear(system);
  
  if (0) {
    nmod_mat_to_stream(den_found, ss);
    std::cout << ss.str() << '\n';
  }

  // double check solution
  nmod_poly_init_mod(den_found_poly, ctx.q_nmod());
  nmod_poly_from_nmod_mat(den_found_poly, den_found);

  nmod_poly_t den_inv, temp;
  nmod_poly_init_mod(den_inv, ctx.q_nmod());
  nmod_poly_invmod(den_inv, den_poly, ctx.modulus);
  nmod_poly_mul(den_found_poly, den_found_poly, den_inv);
  nmod_poly_rem(den_found_poly, den_found_poly, ctx.modulus);

  const char* X = "x";
  nmod_poly_print_pretty(den_found_poly, X);
  std::cout << '\n';

  int deg = nmod_poly_degree(den_found_poly);
  nmod_poly_init_mod(temp, ctx.q_nmod());
  nmod_poly_set_coeff_ui(temp, deg, 1);
  std::cout << "# time: " << duration.count()/1000000.0 << ", mem: "
    "# total: " << rss/1000.0 << " M, " << std::endl;
  
  if (nmod_poly_equal(den_found_poly, temp) == 1) {
    std::cout << "# Success" << std::endl;
  }
  else {
    nmod_poly_neg(temp, temp);
     if (nmod_poly_equal(den_found_poly, temp) == 1) {
       std::cout << "# Success" << std::endl;       
     }
     else {
       std::cout << "# Fail" << std::endl;       
     }
  }

  nmod_poly_clear(temp);
  nmod_poly_clear(den_inv);
  nmod_poly_clear(den_poly);
  nmod_poly_clear(den_found_poly);  
  nmod_mat_clear(den);
  nmod_mat_clear(den_found);  
  nmod_mat_clear(H_mat);
}

// See https://github.com/p-ranav/argparse for argparse help
int main(int argc, char** argv) {
  argparse::ArgumentParser program("arora-ge-ntru");
  program.add_description("Arora-Ge algorithm for NTRU with multiple keys.");
  program.add_argument("n")
    .help("degree of the underlying ring")
    .scan<'i', int>();
  program.add_argument("q")
    .help("prime modulus")
    .scan<'i', int>();
  program.add_argument("-c", "--coeffs")
    .default_value(2)
    .help("number of coefficients. 2 for binary, 3 for ternary")
    .scan<'i', int>();
  program.add_argument("-s", "--seed")
    .default_value(-1)
    .help("optionally fix seed. If seed is -1 then use a random seed.")
    .scan<'i', int>();
  program.add_argument("-r", "--ring")
    .help("use 1 for NTRU: x^n - 1, 2 for NTRU2: x^n + 1, 3 for NTRUPrime: x^n - x - 1 or 4 for NTTRU: x^n - x^(n/2) + 1.")
    .default_value(1)
    .scan<'i', int>();
  program.add_argument("--verbose")
    .help("increase output verbosity")
    .flag();
  
  argparse::ArgumentParser keygen_cmd("keygen");
  keygen_cmd.add_description("Generate NTRU keys with a shared denominator.");
  
  keygen_cmd.add_argument("-k", "--keys")
    .default_value(1)
    .help("number of keys to generate")
    .scan<'i', int>();
  keygen_cmd.add_argument("--pk_output")
    .help("optional output file for public keys");
  keygen_cmd.add_argument("--sk_output")
    .help("optional output file for shared denominator");
  keygen_cmd.add_argument("-s", "--seed")
    .default_value(1)
    .help("optional seed")
    .scan<'i', int>();
  program.add_subparser(keygen_cmd);

  argparse::ArgumentParser system_cmd("system");
  system_cmd.add_description("Create linearized system.");

  system_cmd.add_argument("-i", "--input")
    .required()
    .help("input file of keys (output of keygen subcommand)");
  system_cmd.add_argument("-o", "--output")
    .help("optional output file");  
  program.add_subparser(system_cmd);
  
  argparse::ArgumentParser recover_cmd("recover");
  recover_cmd.add_description("Recover key from linearized system.");
  
  recover_cmd.add_argument("-i", "--input")
    .required()
    .help("input file of linearized system (output of system subcommand)");
  recover_cmd.add_argument("-o", "--output")
    .help("optional output file");
  recover_cmd.add_argument("--nullonly")
    .help("flag -- only output nullspace and then stop")
    .flag();
  program.add_subparser(recover_cmd);

  argparse::ArgumentParser verify_cmd("verify");
  verify_cmd.add_description("Verify the files contain secret keys which are rotations of each other.");
  verify_cmd.add_argument("--sk_input1")
    .required()
    .help("first input file with (a rotation of) the private key");
  verify_cmd.add_argument("--sk_input2")
    .required()
    .help("second input file with (a rotation of) the private key");
  program.add_subparser(verify_cmd);
  
  argparse::ArgumentParser all_cmd("all");
  all_cmd.add_description("All-in-one command.");  
  all_cmd.add_argument("-k", "--keys")
    .default_value(1)
    .help("number of keys to generate")
    .scan<'i', int>();
  program.add_subparser(all_cmd);

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::exception& err) {
    if (program.is_subcommand_used("keygen")) {
      std::cerr << keygen_cmd;
      std::exit(1);
    } else if (program.is_subcommand_used("system")) {
      std::cerr << system_cmd;
      std::exit(1);
    } else if (program.is_subcommand_used("recover")) {
      std::cerr << recover_cmd;
      std::exit(1);
    } else if (program.is_subcommand_used("verify")) {
      std::cerr << verify_cmd;
      std::exit(1);
    } else if (program.is_subcommand_used("all")) {
      std::cerr << all_cmd;
      std::exit(1);
    } else {
      std::cerr << program;
      std::exit(1);
    }
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  int n = program.get<int>("n");
  int q = program.get<int>("q");
  int c = program.get<int>("--coeffs");
  int s = program.get<int>("--seed");
  int r = program.get<int>("--ring");

  assert(n > 1);
  assert(q > 2);
  assert(c == 2 || c == 3);
  assert(r == 1 || r == 2 || r == 3 || r == 4);

  int level = 0;
  if (program["--verbose"] == true) {
    level = 1;
    set_log_level(level);
  }
  
  NTRUKeyGen ctx(n, q, c, r, s, level);

  debug("Parameters:", 
    "\n  n = ", n,
    "\n  q = ", q,
    "\n  coeffs = ", c,
    "\n  ring = ", r,
    "\n  seed = ", s,
    "\n"
  );
  
  if (program.is_subcommand_used("keygen")) {
    keygen(keygen_cmd, ctx);
  } else if (program.is_subcommand_used("system")) {
    system(system_cmd, ctx);
  } else if (program.is_subcommand_used("recover")) {
    recover(recover_cmd, ctx);
  } else if (program.is_subcommand_used("verify")) {
    verify(verify_cmd, ctx);
  } else if (program.is_subcommand_used("all")) {
    all(all_cmd, ctx);
  } else {
    std::cerr << program;
    std::exit(1);    
  }
  
  return 0;
}
