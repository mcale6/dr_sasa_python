#include "common.hpp"
#include <stdexcept>

class SASAError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};
