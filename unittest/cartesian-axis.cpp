//
// Copyright (c) 2017 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#include "pinocchio/spatial/cartesian-axis.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

using namespace se3;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(basics)
{
  typedef Eigen::Vector3d Vector3;
  const Vector3 v = Vector3::Random();
  BOOST_CHECK(Vector3::UnitX().isApprox((Vector3)UnitXd()));
  BOOST_CHECK(Vector3::UnitX().cross(v).isApprox(UnitXd::cross(v)));
  
  BOOST_CHECK(Vector3::UnitY().isApprox((Vector3)UnitYd()));
  BOOST_CHECK(Vector3::UnitY().cross(v).isApprox(UnitYd::cross(v)));
  
  BOOST_CHECK(Vector3::UnitZ().isApprox((Vector3)UnitZd()));
  BOOST_CHECK(Vector3::UnitZ().cross(v).isApprox(UnitZd::cross(v)));
}

BOOST_AUTO_TEST_SUITE_END()
