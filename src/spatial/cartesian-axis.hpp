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

#ifndef __se3_spatial_cartesian_axis_hpp__
#define __se3_spatial_cartesian_axis_hpp__

#include <Eigen/Dense>
#include "pinocchio/spatial/fwd.hpp"

namespace se3
{
  ///
  /// \brief Base class for Cartesian Axis. A Cartesian Axis is a canonic Base vector.
  ///
  template<typename Derived>
  struct CartesianAxisBase
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::Vector Vector;
    enum
    {
      axis = traits<Derived>::axis,
      Options = traits<Derived>::Options,
      dim = traits<Derived>::dim
    };
    
    operator typename Vector::BasisReturnType () const
    { return Vector::Unit(axis); }
    
    operator Vector () const { return (Vector)(Vector::Unit(axis)); }
  };
  
  
  template<int _axis, int _dim, typename _Scalar, int _Options>
  struct traits< CartesianAxisTpl<_axis,_dim,_Scalar,_Options> >
  {
    enum
    {
      axis = _axis,
      Options = _Options,
      dim = _dim
    };
    
    typedef _Scalar Scalar;
    typedef Eigen::Matrix<Scalar,dim,1,Options> Vector;
  };
  
  ///
  /// \brief Generic templated version of a Cartesian Axis.
  ///
  /// \tparam _axis Index of the canonic axis.
  /// \tparam _dim  Dimension of the space.
  /// \tparam _Scalar Scalar type of the vector space.
  /// \tparam _Options Eigen option for matrices.
  ///
  template<int _axis, int _dim, typename _Scalar, int _Options>
  struct CartesianAxisTpl : CartesianAxisBase< CartesianAxisTpl<_axis,_dim,_Scalar,_Options> >
  {
    typedef CartesianAxisBase< CartesianAxisTpl<_axis,_dim,_Scalar,_Options> > Base;
  };
  
  ///
  /// \brief Special treatment for the 3d case.
  ///
  template<typename Derived>
  struct CartesianAxis3Base;
  
  template<typename Derived>
  struct traits< CartesianAxis3Base<Derived> > : traits<Derived>
  {
    typedef typename traits<Derived>::Scalar Scalar;
    typedef typename traits<Derived>::Vector Vector;
    
  };
  
  template<typename Derived>
  struct CartesianAxis3Base : CartesianAxisBase< Derived >
  {
    typedef CartesianAxisBase< CartesianAxis3Base<Derived> > Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Vector Vector;
    
    template<typename VectorDerived, typename ReturnType>
    static void cross(const Eigen::MatrixBase<VectorDerived> & v,
                      Eigen::MatrixBase<ReturnType> const & res)
    { Derived::cross(v,res); }
    
    template<typename VectorDerived>
    static Vector cross(const Eigen::MatrixBase<VectorDerived> & v)
    { Vector res; Derived::cross(v,res); return res; }
  };
  
  template<typename _Scalar, int _Options>
  struct CartesianAxisTpl<0,3,_Scalar,_Options> : CartesianAxis3Base< CartesianAxisTpl<0,3,_Scalar,_Options> >
  {
    typedef CartesianAxis3Base< CartesianAxisTpl<0,3,_Scalar,_Options> > Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Vector Vector;
    
    using Base::cross;
   
    template<typename VectorDerived, typename ReturnType>
    static void cross(const Eigen::MatrixBase<VectorDerived> & v,
                      Eigen::MatrixBase<ReturnType> const & res)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VectorDerived,3);
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(ReturnType,3);
      
      Eigen::MatrixBase<ReturnType> & res_ = const_cast<Eigen::MatrixBase<ReturnType> &>(res);
      res_[0] = 0.; res_[1] = -v[2]; res_[2] = v[1];
    }
  };
  
  template<typename _Scalar, int _Options>
  struct CartesianAxisTpl<1,3,_Scalar,_Options> : CartesianAxis3Base< CartesianAxisTpl<1,3,_Scalar,_Options> >
  {
    typedef CartesianAxis3Base< CartesianAxisTpl<1,3,_Scalar,_Options> > Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Vector Vector;
    
    using Base::cross;
    
    template<typename VectorDerived, typename ReturnType>
    static void cross(const Eigen::MatrixBase<VectorDerived> & v,
                      Eigen::MatrixBase<ReturnType> const & res)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VectorDerived,3);
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(ReturnType,3);
      
      Eigen::MatrixBase<ReturnType> & res_ = const_cast<Eigen::MatrixBase<ReturnType> &>(res);
      res_[0] = v[2]; res_[1] = 0.; res_[2] = -v[0];
    }
  };
  
  template<typename _Scalar, int _Options>
  struct CartesianAxisTpl<2,3,_Scalar,_Options> : CartesianAxis3Base< CartesianAxisTpl<2,3,_Scalar,_Options> >
  {
    typedef CartesianAxis3Base< CartesianAxisTpl<2,3,_Scalar,_Options> > Base;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Vector Vector;
    
    using Base::cross;
    
    template<typename VectorDerived, typename ReturnType>
    static void cross(const Eigen::MatrixBase<VectorDerived> & v,
                      Eigen::MatrixBase<ReturnType> const & res)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(VectorDerived,3);
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(ReturnType,3);
      
      Eigen::MatrixBase<ReturnType> & res_ = const_cast<Eigen::MatrixBase<ReturnType> &>(res);
      res_[0] = -v[1]; res_[1] = v[0]; res_[2] = 0.;
    }
  };
  
} // namespace se3

#endif // ifndef __se3_spatial_cartesian_axis_hpp_
