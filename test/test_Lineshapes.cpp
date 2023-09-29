#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE "Lineshapes"

#include <boost/test/unit_test.hpp>
#include "AmpGen/Lineshapes.h"
#include "AmpGen/CompiledExpression.h"

#include "AmpGen/simd/avx2d_types.h"
using namespace AmpGen::AVX2d;

namespace tt = boost::test_tools;
using namespace AmpGen;

BOOST_AUTO_TEST_CASE(test_BW)
{
  DebugSymbols db;

  double s = 1.0;
  auto expression = Lineshape::BW().get(Parameter("x1[0]", s, true), 0.1, 0.1, "rho(770)0", 1, "", &db);
  auto compiled_expression = make_expression_with_db<std::complex<double>>(expression, "expression", db);
  auto z = compiled_expression(&s);

  BOOST_TEST(std::real(z) == -0.27051973, tt::tolerance(1e-5));
  BOOST_TEST(std::imag(z) == 0.26100285, tt::tolerance(1e-5));
  BOOST_TEST(std::real(expression()) == -0.27051973, tt::tolerance(1e-5));
  BOOST_TEST(std::imag(expression()) == 0.26100285, tt::tolerance(1e-5));

  compiled_expression.debug(&s);
}
