
#pragma once

#include <iostream>
#include <cstdarg>

static int log_level = 0;

inline void set_log_level(int level) {
  log_level = level;
}

template<typename ...T>
inline void debug(T&... msg) {
  if (log_level > 0)
    (std::cout << ... << msg);
}