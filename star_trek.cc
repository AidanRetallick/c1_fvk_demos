//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
#include <fenv.h>

//Generic routines
#include "generic.h"

// The equations
#include "c1_foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"
#include "src/generic/elements.h"
#include "src/generic/my_geom_object.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

//                      OUTLINE OF PROBLEM CONSTRUCTION
// The basic constuction is much the same as the usual order of things in a
// problem. Underneath is the order of actions (with stars next to non actions
// that are unique to these types of problems).
// 1.  Setup mesh parameters
// 2.  Build the mesh

// hierher: combine 3 and 4?

// 3.* Upgrade Elements
//     We upgrade edge elements on relevant boundaries to be curved C1 elements.
//     This involves working out which edge is to be upgraded and then passing
//     information about the global curve and start and end points of the
//     element edge on that curve to the element.
// 4.* Rotate edge degrees of freedom.
//     We rotate the Hermite dofs that lie on the edge into the normal -
//     tangential basis so that we can set physical boundary conditions like
//     clamping or resting conditions.
// 5.  Complete problem Setup and set Boundary conditions.

//                            REQUIRED DEFINITIONS
// Per Curve Section we will need:
// 1.  A parametric function defining the curve section.
// 2.  The tangential derivative of the parametric function defining
//     the curve section.


// hierher: we should probably insist on this being done in any case. Is it
//          obvious to the user that they have a fifth order boundary
//          representation?

// 3.* (For 5 order boundary representation) The second tangential derivative
//     of the parametric function defining the curve section.
// 4.  A unit normal and tangent to each curve section and corresponding
//     derivatives, to allow the rotation of boundary coordinates.

// hierher: doesn't really make sense and how is it convenient? It's
//          not used here. I think

// It also convenient to define:
// 1.  An inverse function (x,y) -> s (the arc coordinate) to help in setting
//     up the nodal positions in terms of this parametric coordinate.




//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////



namespace oomph
{
  //========= start_of_duplicate_node_constraint_element ==================
  /// Non-geometric element used to constrain dofs between duplicated
  /// vertices where the Hemite data at each node is different but must
  /// be compatible.
  //=======================================================================
  class ConstrainDuplicateHermiteFvKNodeElement : public virtual GeneralisedElement
  {
    
  public:
    /// Construcor. Needs the two node pointers so that we can retrieve the
    /// boundary data at solve time
    ConstrainDuplicateHermiteFvKNodeElement(Node* left_node_pt,
				       Node* right_node_pt)
      : Left_node_pt(left_node_pt),
	Right_node_pt(right_node_pt)
    {
      // Add internal data which stores the eight Lagrange multipliers
      add_internal_data(new Data(8));
      
      // Add each node as external data
      add_external_data(Left_node_pt);
      add_external_data(Right_node_pt);
    }

    /// Destructor
    ~ConstrainDuplicateHermiteFvKNodeElement()
    {
      // Must remove Lagrange multiplier data?
    }

    /// Add the contribution to the residuals from the Lagrange multiplier
    /// constraining equations
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      fill_in_generic_residual_contribution_foeppl_von_karman(
	residuals,
	GeneralisedElement::Dummy_matrix,
	0);
    }

    /// Add the contribution to the Jacobian from the Lagrange multiplier
    /// constraining equations
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
					 DenseMatrix<double> &jacobian)
    {
      fill_in_generic_residual_contribution_foeppl_von_karman(
	residuals,
	jacobian,
	1);
    }
    
  private:
    /// Add the contribution to the residuals (and jacobain if flag is 1) from
    /// the Lagrange multiplier constraining equations
    void fill_in_generic_residual_contribution_foeppl_von_karman(
      Vector<double> &residuals,
      DenseMatrix<double> &jacobian,
      const unsigned &flag)
    {
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Calculate Jacobian and Hessian of coordinate transform between
      // each boundary coordinate

      //----------------------------------------------------------------------
      // We need the parametrisations either side of the vertex which define
      // the coordinates each node uses for its Hermite dofs.
      Vector<double> left_dxids(2,0.0);
      Vector<double> left_d2xids2(2,0.0);
      Vector<double> right_dxids(2,0.0);
      Vector<double> right_d2xids2(2,0.0);
      Left_boundary_pt->dposition(Left_node_coord, left_dxids);
      Left_boundary_pt->d2position(Left_node_coord, left_d2xids2);
      Right_boundary_pt->dposition(Right_node_coord, right_dxids);
      Right_boundary_pt->d2position(Right_node_coord, right_d2xids2);

      // Get the speed of each parametrisation
      double left_mag =
	sqrt(left_dxids[0] * left_dxids[0] + left_dxids[1] * left_dxids[1]);
      double right_mag =
	sqrt(right_dxids[0] * right_dxids[0] + right_dxids[1] * right_dxids[1]);

      //----------------------------------------------------------------------
      // Normalise dxids to find the tangent vectors, normal vectors and their
      // derivatives either side of the vertex
      Vector<double> left_ti(2,0.0);
      Vector<double> left_ni(2,0.0);
      Vector<double> left_dtids(2,0.0);
      Vector<double> left_dnids(2,0.0);
      Vector<double> right_ti(2,0.0);
      Vector<double> right_ni(2,0.0);
      Vector<double> right_dtids(2,0.0);
      Vector<double> right_dnids(2,0.0);
      for(unsigned alpha=0; alpha<2; alpha++)
      {
	// Fill in the tangents either side of the vertex
	left_ti[alpha]  =  left_dxids[alpha] / left_mag;
	right_ti[alpha] = right_dxids[alpha] / right_mag;
	// Fill in the derivatives of the (normalised) tangents either side of
	// the vertex
	left_dtids[alpha] = left_d2xids2[alpha] / std::pow(left_mag, 2)
	  - (left_dxids[0]*left_d2xids2[0]+left_dxids[1]*left_d2xids2[1])
	  * left_dxids[alpha] / std::pow(left_mag,4);
	right_dtids[alpha] = right_d2xids2[alpha] / std::pow(right_mag, 2)
	  - (right_dxids[0]*right_d2xids2[0]+right_dxids[1]*right_d2xids2[1])
	  * right_dxids[alpha] / std::pow(right_mag,4);
	// Use these to fill out the corresponding vectors for the normal
	// direction (nx,ny) = (ty,-tx)
	left_ni[alpha] = pow(-1,alpha)*left_ti[(alpha+1)%2];
	right_ni[alpha] = pow(-1,alpha)*right_ti[(alpha+1)%2];
	left_dnids[alpha] = pow(-1,alpha)*left_dtids[(alpha+1)%2];
	right_dnids[alpha] = pow(-1,alpha)*right_dtids[(alpha+1)%2];
      }

      //----------------------------------------------------------------------
      // We need to fill out the Jacobians and Hessians of the boundary
      // coordinates either side of the vertex
      DenseMatrix<double> left_jac(2,2,0.0);
      DenseMatrix<double> right_jac(2,2,0.0);
      Vector<DenseMatrix<double>> left_hess(2,DenseMatrix<double>(2,2,0.0));
      Vector<DenseMatrix<double>> right_hess(2,DenseMatrix<double>(2,2,0.0));
      for(unsigned alpha=0; alpha<2; alpha++)
      {
	// Fill in Jacobians {{nx,tx},{ny,ty}}
	left_jac(alpha,0) = left_ni[alpha];
	left_jac(alpha,1) = left_ti[alpha];
	right_jac(alpha,0) = right_ni[alpha];
	right_jac(alpha,1) = right_ti[alpha];
	// Fill in Hessians
	// left_hess[alpha](0,0) = 0.0;
	left_hess[alpha](0,1) = left_dnids[alpha];
	left_hess[alpha](1,0) = left_dnids[alpha];
	left_hess[alpha](1,1) = left_dtids[alpha];
	// right_hess[alpha](0,0) = 0.0;
	right_hess[alpha](0,1) = right_dnids[alpha];
	right_hess[alpha](1,0) = right_dnids[alpha];
	right_hess[alpha](1,1) = right_dtids[alpha];
      }

      //----------------------------------------------------------------------
      // We need the inverse Jacobian and Hessian for the left parametrisation
      DenseMatrix<double> left_jac_inv(2,2,0.0);
      Vector<DenseMatrix<double>> left_hess_inv(2,DenseMatrix<double>(2,2,0.0));
      left_jac_inv(0,0) = left_jac(1,1);
      left_jac_inv(0,1) =-left_jac(0,1);
      left_jac_inv(1,0) =-left_jac(1,0);
      left_jac_inv(1,1) = left_jac(0,0);
      // Fill out inverse of Hessian
      // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
	for (unsigned beta = 0; beta < 2; beta++)
	{
	  for (unsigned gamma = 0; gamma < 2; gamma++)
	  {
	    for (unsigned alpha2 = 0; alpha2 < 2; alpha2++)
	    {
	      for (unsigned beta2 = 0; beta2 < 2; beta2++)
	      {
		for (unsigned gamma2 = 0; gamma2 < 2; gamma2++)
		{
		  left_hess_inv[alpha](beta, gamma) -=
		    left_jac_inv(alpha, alpha2)
		    * left_hess[alpha2](beta2, gamma2)
		    * left_jac_inv(beta2, beta)
		    * left_jac_inv(gamma2, gamma);
		}
	      }
	    }
	  }
	}
      }

      //----------------------------------------------------------------------
      // Use these to calculate the Jacobian of the left->right transform
      //     J = J_{left}^{-1}J_{right}
      // and the Hessian of the left->right transform
      //     H = H_{left}^{-1}J_{right}J_{right} + J_{left}^{-1}H_{right}
      DenseMatrix<double> jac_of_transform(2,2,0.0);
      Vector<DenseMatrix<double>> hess_of_transform(2,DenseMatrix<double>(2,2,0.0));
      for(unsigned alpha = 0; alpha < 2; alpha++)
      {
	for(unsigned beta = 0; beta < 2; beta++)
	{
	  for(unsigned gamma = 0; gamma < 2; gamma++)
	  {
	    // Add contribution to J
	    jac_of_transform(alpha,beta) +=
	      left_jac_inv(alpha,gamma) * right_jac(gamma,beta);
	    for(unsigned mu = 0; mu < 2; mu++)
	    {
	      // Add second term contribution to H
	      hess_of_transform[alpha](beta,gamma) +=
		left_jac_inv(alpha,mu) * right_hess[mu](beta,gamma);
	      for(unsigned nu = 0; nu < 2; nu++)
	      {
		// Add first term contribution to H
		hess_of_transform[alpha](beta,gamma) +=
		  left_hess_inv[alpha](mu,nu)
		  * right_jac(mu,beta) * right_jac(nu,gamma);
	      }
	    }
	  }
	}
      }
      
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Finally, we use jac and hess of transform to add the residual
      // contribution from the constraint

      // First three (Lagrange) dofs are equal
      for(unsigned i_dof = 0; i_dof < 3; i_dof++)
      {
	residuals[i_dof] = 
      }

      // If flag, then add the jacobian contribution
    }
    
    /// Pointer to the left node (before the vertex when anticlockwise)
    Node* Left_node_pt;
    
    /// Pointer to the right node (after the vertex when anticlockwise)
    Node* Right_node_pt;

    /// Pointer to the left node's boundary parametrisation
    CurvilineGeomObject* Left_boundary_pt;

    /// Pointer to the right node's boundary parametrisation
    CurvilineGeomObject* Right_boundary_pt;

    /// Coordinate of the left node on the left boundary
    Vector<double> Left_node_coord;

    /// Coordinate of the left node on the left boundary
    Vector<double> Right_node_coord;
  };

  //========= start_of_point_force_and_torque_wrapper======================
  /// Class to impose point force and torque to (wrapped) Fvk element
  //=======================================================================
  template<class ELEMENT> 
  class FvKPointForceAndSourceElement : public virtual ELEMENT
  {

  public:

    /// Constructor
    FvKPointForceAndSourceElement()
    {
      // Add internal Data to store forces (for now)
      oomph_info << "# of internal Data objects before: " <<
	this->ninternal_data() << std::endl;
    }
 
    /// Destructor (empty)
    ~FvKPointForceAndSourceElement(){}
 
    /// Set local coordinate and point force and_torque
    void setup(const Vector<double>& s_point_force_and_torque)
    {
      S_point_force_and_torque=s_point_force_and_torque;
    }
 
 
    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double> &residuals)
    {
      // Call the generic residuals function with flag set to 0 hierher why no arg?
      // using a dummy matrix argument
      ELEMENT::fill_in_generic_residual_contribution_foeppl_von_karman(
	residuals,
	GeneralisedElement::Dummy_matrix,
	0);
   
      // fill_in_contribution_to_residuals(residuals);

      // Add point force_and_torque contribution
      fill_in_point_force_and_torque_contribution_to_residuals(residuals);
    }

 

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double> &residuals,
					  DenseMatrix<double> &jacobian)
    {
      //Call the generic routine with the flag set to 1 hierher why no arg?
      // ELEMENT::fill_in_contribution_to_jacobian(residuals,
      //                                           jacobian);
      ELEMENT::fill_in_generic_residual_contribution_foeppl_von_karman(
	residuals,
	jacobian,
	1);
   
      // Add point force_and_torque contribution
      fill_in_point_force_and_torque_contribution_to_residuals(residuals);
    }
 

  private:


 
    /// Add the point force_and_torque contribution to the residual vector
    void fill_in_point_force_and_torque_contribution_to_residuals(Vector<double> &residuals)
    {
      // No further action
      if (S_point_force_and_torque.size()==0)
      {
	oomph_info << "bailing" << std::endl;
	return;
      }
   

      //Find out how many nodes there are
      const unsigned n_node = this->nnode();
   
      //Set up memory for the shape/test functions
      Shape psi(n_node);
   
      //Integers to store the local equation and unknown numbers
      // int local_eqn_real=0;
      // int local_eqn_imag=0; 
   
      // Get shape/test fcts
      this->shape(S_point_force_and_torque,psi);
   
      //  // Assemble residuals
      //  //--------------------------------
   
      //  // Loop over the test functions
      //  for(unsigned l=0;l<n_node;l++)
      //   {
      //    // first, compute the real part contribution 
      //    //-------------------------------------------
     
      //    //Get the local equation
      //    local_eqn_real = this->nodal_local_eqn(l,this->u_index_helmholtz().real());
     
      //    /*IF it's not a boundary condition*/
      //    if(local_eqn_real >= 0)
      //     {
      //      residuals[local_eqn_real] += Point_force_and_torque_magnitude.real()*psi(l);
      //     }     
     
      //    // Second, compute the imaginary part contribution 
      //    //------------------------------------------------
     
      //    //Get the local equation
      //    local_eqn_imag = this->nodal_local_eqn(l,this->u_index_helmholtz().imag());
     
      //    /*IF it's not a boundary condition*/
      //    if(local_eqn_imag >= 0)
      //     {
      //      // Add body force/force_and_torque term and Helmholtz bit
      //      residuals[local_eqn_imag] += Point_force_and_torque_magnitude.imag()*psi(l);
      //     }     
      //   }

    }

 
    /// Local coordinates of point at which point force_and_torque is applied
    Vector<double> S_point_force_and_torque;
 
  };
 


  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FvKPointForceAndSourceElement<ELEMENT> > 
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };


  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as 
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<FvKPointForceAndSourceElement<ELEMENT> > >
    : public virtual FaceGeometry<FaceGeometry<ELEMENT> >
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
  };


}

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//========================================================================
/// Namespace for problem parameters
//========================================================================
namespace Parameters
{

  // // Enumeration of cases
  // enum{
  //   Clamped_validation,
  //   Axisymmetric_shear_buckling,
  //   Nonaxisymmetric_shear_buckling
  // };

  /// Which case are we doing
  //  unsigned Problem_case = ; // Nonaxisymmetric_shear_buckling; // Axisymmetric_shear_buckling;
 
  /// Upper ellipse x span
  double A1 = 0.5;

  /// Upper ellipse y span
  double B1 = 1.5;
  
  /// Lower ellipse x span
  double A2 = 0.6;

  /// Lower ellipse y span
  double B2 = 0.7;

  /// Char array containing the condition on each boundary. Character array
  /// index corresponds to boundary enumeration and the entry to the contition
  /// type: // [zdec] use this at some point
  ///   - 'c' for clamped
  ///   - 'p' for pinned
  ///   - 's' for sliding
  ///   - 'f' for free
  char Boundary_conditions[3] = "cc";  
  
  /// Poisson ratio
  double Nu = 0.5;
  
  /// FvK parameter
  double Eta = 120.0e3;
 
  /// Pressure magnitude
  double P_mag = 0.0;
 
  /// In-plane traction magnitude
  double T_mag = 0.0;

  /// Order of the polynomial interpolation of the boundary
  unsigned Boundary_order = 3;
  
  //----------------------------------------------------------------------------
  // [zdec] DEPENDENT VARIABLES -- need some sort of update hook

  /// x-component of intersection (positive value)
  double X_intersect = sqrt(A1*A1*A2*A2*(B1*B1-B2*B2)
			    / (A2*A2*B1*B1 - A1*A1*B2*B2));
  
  /// y-component of boundary intersection
  double Y_intersect = sqrt( B1*B1
			     * ( 1.0 - X_intersect*X_intersect/(A1*A1) ) );

  /// Theta of intersect 1
  // (shift by -pi/2 as the ellipse uses the angle about positive y)
  double Theta1 = atan2(Y_intersect/B1,X_intersect/A1) - Pi/2.0;

  /// Theta of intersect 2
  // (shift by -pi/2 as the ellipse uses the angle about positive y)
  double Theta2 = atan2(Y_intersect/B2,X_intersect/A2) - Pi/2.0;
  

  // Boundary info
  //                       __
  //                     -    -
  //                   -        -   *Upper ellipse arc*
  //                 /            \
  //               /                \
  //             /                    \
  //           /                        \
  //          /         ________         \
  //         /       --         --        \
  //       /       --              --       \
  //      /      -                    -      \
  //     /    /  *Lower ellipse arc*     \    \
  //    /  /                                \  \
  //   / /                                    \ \
  //   X(Theta2)                                X(Theta1)
  
  /// Parametric curve for the upper elliptical boundary arc
  // (true->anticlockwise parametrisation)
  CurvilineEllipseTop Upper_parametric_elliptical_curve(A1,B1,false);

  /// Parametric curve for the lower elliptical boundary arc
  // (true->clockwise parametrisation)
  CurvilineEllipseTop Lower_parametric_elliptical_curve(A2,B2,true);

  /// Vector of parametric boundaries
  Vector<CurvilineGeomObject*> Parametric_curve_pt =
  {
     &Upper_parametric_elliptical_curve,
    &Lower_parametric_elliptical_curve
  };
  
  
  
  /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    pressure = P_mag;
  }


 
  /// In plane forcing (shear stress) depending on the position (x,y)
  void get_in_plane_force(const Vector<double>& x, Vector<double>& tau)
  {
    tau[0] = 0.0;
    tau[1] = 0.0;
  }


  // hierher: kill but check with Aidan first
 
  // // This metric will flag up any non--axisymmetric parts
  // void axiasymmetry_metric(const Vector<double>& x,
  //        		  const Vector<double>& u,
  //        		  const Vector<double>& u_exact,
  //        		  double& error,
  //        		  double& norm)
  // {
  //  // We use the theta derivative of the out of plane deflection
  //  error = pow((-x[1]*u[1] + x[0]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
  //  norm  = pow(( x[0]*u[1] + x[1]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
  // }

  // Get the null function for applying homogenous BCs
  void get_null_fct(const Vector<double>& X, double& exact_w)
  {
    exact_w = 0.0;
  }

  // Get the unit function for applying homogenous BCs
  void get_unit_fct(const Vector<double>& X, double& exact_w)
  {
    exact_w = 1.0;
  }

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Problem definition
//====================================================================
using Parameters::Parametric_curve_pt;

template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

  /// Constructor
  UnstructuredFvKProblem(double const& element_area = 0.09);

  /// Destructor
  ~UnstructuredFvKProblem()
  {
    // Close trace file
    Trace_file.close();
  
    // Clean up memory
    delete Bulk_mesh_pt;
    delete Outer_boundary_pt;
    delete Outer_curvilinear_boundary_pt[0];
    delete Outer_curvilinear_boundary_pt[1];
  };

  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve: empty
  void actions_before_newton_solve(){}

  /// Doc the solution
  void doc_solution(const std::string& comment="");

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return Bulk_mesh_pt; //dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt());
  }


private:

  // [zdec] this is not necessarily inside the domain for this geometry, needs
  // fixing
  //
  // /// Pin all displacements and rotation (dofs 0-4) at the centre
  // void pin_all_displacements_and_rotation_at_centre_node();

  /// Setup and build the mesh
  void build_mesh();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of  all elements
  void complete_problem_setup();

  /// Loop over all curved edges, then loop over elements and upgrade
  /// them to be curved elements
  void upgrade_edge_elements_to_curve(const unsigned &b);

  /// Duplicate nodes at corners in order to properly apply boundary
  /// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
  /// problem in order to constrain continuous interpolation here across its (8)
  /// vertex dofs. (Note "corner" here refers to the meeting point of any two
  /// sub-boundaries in the closed external boundary)
  void duplicate_corner_nodes();

  /// Lookup scheme for duplicated nodes -- left meaning before the vertex when
  /// traversing the boundary in an anticlockwise manner
  Vector<Node*> Left_duplicated_node_pt;

  /// Lookup scheme for duplicated nodes -- right meaning after the vertex when
  /// traversing the boundary in an anticlockwise manner
  Vector<Node*> Right_duplicated_node_pt;
  
  /// Vector of pointers to the lagrange data associated with each pair of
  /// duplicated nodes one Lagrange multiplier per dof type resulting in six.
  /// If the first node uses coordinates (s_1,s_2) and the second uses
  /// coordinates (t_1, t_2) then we have one equation for displacement:
  ///     0 = L_1 * (w_1-w_2)
  /// two equations constraining gradient:
  ///     0 = L_{1+\alpha} * (dw_1/ds_\alpha - J_{\alpha\beta} dw_2/dt_\beta)
  /// and three equations constraining curvature (beta>=alpha):
  ///     0 = L_{2+\alpha+\beta} * (d^2w_1/ds_\alpha ds_\beta
  ///            - J_{\alpha\gamma} * J_{\beta\delta} * d^2w_2/dt_\gamma dt_\delta
  ///            - H_{\gamma\alpha\beta} * dw_2/dt_gamma)
  Vector<Data*> Lagrange_multiplier_pt;

  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in Parameters
  void rotate_edge_degrees_of_freedom();

  /// Delete traction elements and wipe the surface mesh
  void delete_traction_elements(Mesh* const &surface_mesh_pt);
  
  /// Trace file to document norm of solution
  ofstream Trace_file;
  
  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1
  };

  /// Target element area
  double Element_area;

  /// Pointer to "surface" mesh
  Mesh* Surface_mesh_pt;

  /// Doc info object for labeling output
  DocInfo Doc_info;
  
  // /// Outer boundary Geom Object
  // Ellipse* Outer_boundary_ellipse_pt;
 
  /// The outer curves
  Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;
 
  /// The close outer boundary
  TriangleMeshClosedCurve* Outer_boundary_pt;
 
}; // end_of_problem_class



//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(const double& element_area)
  :
  Element_area(element_area)
{
  // Build the mesh
  build_mesh();

  // Curved Edge upgrade
  upgrade_edge_elements_to_curve(Outer_boundary0);
  upgrade_edge_elements_to_curve(Outer_boundary1);

  // Add extra nodes at boundaries and constrain the dofs there.
  duplicate_corner_nodes();
  
  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom();

  // Store number of bulk elements
  complete_problem_setup();

  // Set directory
  Doc_info.set_directory("RESLT");
 
  // Open trace file
  char filename[100];
  sprintf(filename, "RESLT/trace.dat");
  Trace_file.open(filename);

  // Assign equation numbers
  oomph_info << "Number of equations: "
  << assign_eqn_numbers() << '\n';
 
} // end Constructor




//======================================================================
/// Set up and build the mesh
//======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
  Vector<double> zeta(1);
  Vector<double> posn(2);
  
  //Outer boundary
  //--------------
  double theta1=Parameters::Theta1;
  double theta2=Parameters::Theta2;
  
  // [zdec] debug
  oomph_info << theta1 << " " << theta2 << std::endl;
  
  //First bit
  double zeta_start = theta1;
  double zeta_end = -theta1;
  unsigned nsegment = (int)(2.0*(theta2-theta1)/sqrt(Element_area))+2;
  
  Outer_curvilinear_boundary_pt.resize(2);
  Outer_curvilinear_boundary_pt[0] =
    new TriangleMeshCurviLine(Parameters::Parametric_curve_pt[0], zeta_start,
			      zeta_end, nsegment, Outer_boundary0);

  //Second bit (Even though we go from 2->1, the parametrisation is flipped)
  zeta_start = theta2;
  zeta_end = -theta2;
  Outer_curvilinear_boundary_pt[1] =
    new TriangleMeshCurviLine(Parameters::Parametric_curve_pt[1], zeta_start,
			      zeta_end, nsegment, Outer_boundary1);
  
  // Combine
  Outer_boundary_pt =
    new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);

  // [zdec] debug: output mesh boundary
  ofstream mesh_debug;
  mesh_debug.open("boundary_file.dat");
  Outer_boundary_pt->output(mesh_debug, 200);
  mesh_debug.close();
  
  //Create mesh parameters object
  TriangleMeshParameters mesh_parameters(Outer_boundary_pt);

  // Element area
  mesh_parameters.element_area() = Element_area;
  
  // Build an assign bulk mesh
  Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);

  // Reset the non-vertex node positions
  for(unsigned e=0; e<Bulk_mesh_pt->nelement(); e++)
  {
    // [zdec] debug
    oomph_info << "Repairing for element " << e << std::endl;
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    el_pt->repair_lagrange_node_positions();
  }
  
  // Split elements that have two boundary edges
  TimeStepper* time_stepper_pt = Bulk_mesh_pt->Time_stepper_pt;
  Bulk_mesh_pt->
    template split_elements_with_multiple_boundary_edges<ELEMENT>(time_stepper_pt);
    
  //Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);
  
  // Combine submeshes into a single Mesh (over the top; could just have
  // assigned bulk mesh directly.
  build_global_mesh();
 
}// end build_mesh




// //==start_of_pin_all_displacements_and_rotation_at_centre_node======================
// /// pin all displacements and rotations in the centre
// //==============================================================================
// template<class ELEMENT>
// void UnstructuredFvKProblem<ELEMENT>::pin_all_displacements_and_rotation_at_centre_node()
// {


//   // Choose non-centre node on which we'll supress
//   // the rigid body rotation around the z axis.
//   double max_x_potentially_pinned_node=-DBL_MAX;
//   Node* pinned_rotation_node_pt=0;
  
//   // Pin the node that is at the centre in the domain
//   unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
//   for (unsigned inod=0;inod<num_int_nod;inod++)
//   {
//     // Get node point
//     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);

//     // Check which coordinate increases along this boundary
//     // oomph_info << "node: "
//     //            << nod_pt->x(0) << " "
//     //            << nod_pt->x(1) << " "
//     //            << std::endl;

//     // Find the node with the largest x coordinate
//     if (fabs(nod_pt->x(0))>max_x_potentially_pinned_node)
//     {
//       max_x_potentially_pinned_node=fabs(nod_pt->x(0));
//       pinned_rotation_node_pt=nod_pt;
//     }
   
        
//     // If the node is on the other internal boundary too
//     if( nod_pt->is_on_boundary(3))
//     {
//       // Pin it! It's the centre of the domain!
//       // In-plane dofs are always 0 and 1
//       // Out of plane displacement is 2, x and y derivatives are 3 and 4.
//       nod_pt->pin(0);
//       nod_pt->set_value(0,0.0);
//       nod_pt->pin(1);
//       nod_pt->set_value(1,0.0);
//       nod_pt->pin(2);
//       nod_pt->set_value(2,0.0);
//       nod_pt->pin(3);
//       nod_pt->set_value(3,0.0);
//       nod_pt->pin(4);
//       nod_pt->set_value(4,0.0);


//     }
//   }


//   oomph_info << "rotation pinning node: "
//   << pinned_rotation_node_pt->x(0) << " "
//   << pinned_rotation_node_pt->x(1) << " "
//   << std::endl;
//   // Pin y displacement 
//   pinned_rotation_node_pt->pin(1);
 
// }



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{  
  // Set the boundary conditions
  apply_boundary_conditions();
 
  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
  {   
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    
    // Set the pressure function pointers and the physical constants
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->in_plane_forcing_fct_pt() = &Parameters::get_in_plane_force;

    // hierher Aidan: kill this error metric thing in element
    // There is no error metric in this case
    // el_pt->error_metric_fct_pt() = &Parameters::axiasymmetry_metric;

    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_pt() = &Parameters::Eta;
  }

}



//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
  // Clamp it
  //if (Parameters::Problem_case==Parameters::Clamped_validation)
  {
    // Set the boundary conditions
    unsigned n_bound = 2;
    for(unsigned b=0;b<n_bound;b++)
    {
      const unsigned n_b_element = Bulk_mesh_pt->nboundary_element(b);

      // Get all nodes on the boundary by looping over boundary elements and
      // using helper functions
      for(unsigned e=0;e<n_b_element;e++)
      {
	// Get pointer to bulk element adjacent to b
	ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

	// Resting pin so set [0, 2, 5]
	// z-displacement
	el_pt->fix_out_of_plane_displacement_dof(0,b,Parameters::get_null_fct);
	// dz/dn
	// el_pt->fix_out_of_plane_displacement_dof(1,b,Parameters::get_null_fct);
	// dz/dt
	el_pt->fix_out_of_plane_displacement_dof(2,b,Parameters::get_null_fct);
	// d2z/dn2
	// el_pt->fix_out_of_plane_displacement_dof(3,b,Parameters::get_null_fct);
	// d2z/dndt
	// el_pt->fix_out_of_plane_displacement_dof(4,b,Parameters::get_null_fct);
	// d2z/dt2
	el_pt->fix_out_of_plane_displacement_dof(5,b,Parameters::get_null_fct);
      } // End loop over elements [e]

      // On boundary b, apply appropriate boundary conditions at the vertex
      // node between b and [b+1]%n_bound
      // (modulo so we find n_bound-1 -- 0 corner)
      unsigned n_b_node = Bulk_mesh_pt->nboundary_node(b);
      for(unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
      {
	Node* node_pt = Bulk_mesh_pt->boundary_node_pt(b,i_b_node);
	// If also on the next boundary then we are at a corner
	if( node_pt->is_on_boundary((b+1)%n_bound) )
	{
	  // Homogenous pin both sides so all dofs pinned and zero
	  for(unsigned j_w_dof = 0; j_w_dof < 6; j_w_dof++)
	  {
	    // Pin the dof (+2 due to in-plane dofs) and set to 0.0
	    unsigned j_dof = j_w_dof+2;
	    node_pt->pin(j_dof);
	    node_pt->set_value(j_dof, 0.0);
	  } // End loop over each dof [j_w_dof]
	} // End if( at corner )
      } // End loop over boundary nodes on b [i_b_node]
      
    } // End loop over boundaries [b]
    
  }
  // // All other cases: simply pin and stop rotation via the centre
  // else
  // {
  //   pin_all_displacements_and_rotation_at_centre_node();
  // }
 
} // end set bc



//==============================================================================
/// A function that upgrades straight sided elements to be curved. This involves
/// Setting up the parametric boundary, F(s) and the first derivative F'(s)
/// We also need to set the edge number of the upgraded element and the positions
/// of the nodes j and k (defined below) and set which edge (k) is to be exterior
/*            @ k                                                             */
/*           /(                                                               */
/*          /. \                                                              */
/*         /._._)                                                             */
/*      i @     @ j                                                           */
/// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
/// representation. That is we need to have a continuous 2nd derivative defined
/// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
/// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
/// as well.
//==start_of_upgrade_edge_elements==============================================
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT >::
upgrade_edge_elements_to_curve(const unsigned &ibound)
{

  // [zdec] debug
  ofstream debug_stream;
  debug_stream.open("mesh_debug.dat");
  Bulk_mesh_pt->output(debug_stream);
  debug_stream.close();
  
  // Parametric curve describing the boundary
  // hierher code share with thing that defines the
  // boundaries of the mesh!
  CurvilineGeomObject* parametric_curve_pt=0;
  switch (ibound)
  {
  case 0:
    parametric_curve_pt = &Parameters::Upper_parametric_elliptical_curve;
    break;
   
  case 1:
    parametric_curve_pt = &Parameters::Lower_parametric_elliptical_curve;
    break;
   
  default:
    throw OomphLibError("Unexpected boundary number.",
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
    break;
  } // end parametric curve switch
  
  
  // Loop over the bulk elements adjacent to boundary ibound
  const unsigned n_els=Bulk_mesh_pt->nboundary_element(ibound);
  for(unsigned e=0; e<n_els; e++)
  {
    // Get pointer to bulk element adjacent to b
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(ibound,e));
 
    // hierher what is that? why "My"?
    // Initialise enum for the curved edge
    MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

    // Loop over all (three) nodes of the element and record boundary nodes
    unsigned index_of_interior_node = 3;
    unsigned nnode_on_neither_boundary = 0;
    const unsigned nnode = 3;


    // hierher what does this comment mean?
    // Fill in vertices' positions (this step should be moved inside the curveable
    // Bell element)
     Vector<Vector<double> > xn(nnode,Vector<double>(2,0.0));
    for(unsigned n=0;n<nnode;++n)
    {
      Node* nod_pt = bulk_el_pt->node_pt(n);
      xn[n][0]=nod_pt->x(0);
      xn[n][1]=nod_pt->x(1);
     
      // Check if it is on the outer boundaries
      if(!(nod_pt->is_on_boundary(ibound)))
      {
	index_of_interior_node = n;
	++nnode_on_neither_boundary;
      }
    }// end record boundary nodes


    // hierher: ouch! This seems to map (x,y) to zeta! This is at best possible to within
    // a tolerance. Needs a redesign!
   
    // s at the next (cyclic) node after interior
    const double s_ubar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
   
    // s at the previous (cyclic) node before interior
    const double s_obar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);
   
    // Assign edge case
    edge = static_cast<MyC1CurvedElements::Edge>(index_of_interior_node);

    // Check nnode_on_neither_boundary
    if(nnode_on_neither_boundary == 0)
    {
      throw OomphLibError(
	"No interior nodes. One node per CurvedElement must be interior.",
	OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else if (nnode_on_neither_boundary > 1)
    {
      throw OomphLibError(
 	"Multiple interior nodes. Only one node per CurvedElement can be interior.",
	OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Check for inverted elements
    if (s_ubar>s_obar)
    {
      throw OomphLibError(
	"Decreasing parametric coordinate. Parametric coordinate must increase as the edge is traversed anti-clockwise.",
	OOMPH_CURRENT_FUNCTION,
	OOMPH_EXCEPTION_LOCATION);
    } // end checks

    // Upgrade it
    bulk_el_pt->upgrade_element_to_curved(edge,s_ubar,s_obar,parametric_curve_pt,5);
  }
}// end_upgrade_elements



//==============================================================================
/// Duplicate nodes at corners in order to properly apply boundary
/// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
/// problem in order to constrain continuous interpolation here across its (8)
/// vertex dofs. (Note "corner" here refers to the meeting point of any two
/// sub-boundaries in the closed external boundary)
//============================================================================== 
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT >::duplicate_corner_nodes()
{
  // [zdec] WRITEME

  // // Loop over the sections of the external boundary
  // for()
  // {
  //   // Get node and el pts at the boundary vertex
  //   Node* old_node = ;
  //   Node* new_node = NEW Node*;
  //   ELEMENT* left_element = ;
  //   ELEMENT* right_element = ;
  //
  //   // Add new node to right element
  //   
  //   // Update boudary info
  //
  //   // Create Lagrange multiplier dof
  //
  //   // Add nodes and Lagrange multipliers lookup scheme
  // }
}



//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
//======================================================================
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::rotate_edge_degrees_of_freedom()
{
  
  std::ofstream debugstream;
  debugstream.open("test_file");
  debugstream << "x y nx ny tx ty nx_x nx_y ny_x ny_y tx_x tx_y ty_x ty_y"
              << std::endl;
  debugstream.close();
		
  // Get the number of boundaries
  unsigned n_bound = 2;
  // Loop over the bulk elements
  unsigned n_element = Bulk_mesh_pt-> nelement();
  for(unsigned e=0; e<n_element; e++)
  {
    // Get pointer to bulk element adjacent
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
    // Loop over each boundary and add the boundary parametrisation to the
    // relevant nodes' boundary data
    for(unsigned b=0; b<n_bound; b++)
    {
      // Calculate nodes on the relevant boundaries
      const unsigned nnode=3;
      // Count the number of boundary nodes on external boundaries
      Vector<unsigned> boundary_node;
      // Store the boundary coordinates of nodes on the boundaries
      Vector<double> boundary_coordinate_of_node;
      for (unsigned n=0; n<nnode;++n)
      {
	// If on external boundary b
	if (el_pt->node_pt(n)->is_on_boundary(b))
	{
	  boundary_node.push_back(n);
	  double coord = Parametric_curve_pt[b]
	    ->get_zeta(el_pt->node_pt(n)->position());
	  boundary_coordinate_of_node.push_back(coord);
	}
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if(!boundary_node.empty())
      {
	// Rotate the nodes by passing the index of the nodes and the
	// normal / tangent vectors to the element
	el_pt->
	  rotated_boundary_helper_pt()->
	  add_nodal_boundary_parametrisation(boundary_node,
					     boundary_coordinate_of_node,
					     Parametric_curve_pt[b]);
      }
    }
  }
}// end rotate_edge_degrees_of_freedom

 

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const
						   std::string& comment)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 30;

  sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
  << comment << "\"\n";
  some_file.close();
 

  // // Find the solution at r=0
  // // ----------------------

  // // hierher precompute
  // MeshAsGeomObject Mesh_as_geom_obj(Bulk_mesh_pt);
  // Vector<double> s(2);
  // GeomObject* geom_obj_pt=0;
  // Vector<double> r(2,0.0);
  // Mesh_as_geom_obj.locate_zeta(r,geom_obj_pt,s);

  // // Compute the interpolated displacement vector
  // Vector<double> u_0(12,0.0);
  // u_0=dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_foeppl_von_karman(s);
 
  // oomph_info << "w in the middle: " <<std::setprecision(15) << u_0[0] << std::endl;
 
  // Trace_file << Parameters::P_mag << " " << u_0[0] << '\n';

  // Increment the doc_info number
  Doc_info.number()++;

} // end of doc



//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  // Define possible command line arguments and parse the ones that
  // were actually specified


  // Clamped boundary conditions?
  CommandLineArgs::specify_command_line_flag("--use_clamped_bc");

  // Poisson Ratio
  CommandLineArgs::specify_command_line_flag("--nu",
					     &Parameters::Nu);

  // Applied Pressure
  CommandLineArgs::specify_command_line_flag("--p",
					     &Parameters::P_mag);

  // FvK prameter
  CommandLineArgs::specify_command_line_flag("--eta",
					     &Parameters::Eta);

  // Element Area 
  double element_area=0.05;
  CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

  
  // Order of the polynomial interpolation of the boundary
  CommandLineArgs::specify_command_line_flag("--boundary_order",
					     &Parameters::Boundary_order);

  
  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();


  // Build problem
  // UnstructuredFvKProblem<NON_WRAPPED_ELEMENT
  //UnstructuredFvKProblem<FvKPointForceAndSourceElement<NON_WRAPPED_ELEMENT>>
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>> 
    problem(element_area);


  double dp_mag=0.1;
  double dt_mag=0.00;
  unsigned n_step=10;

  for(unsigned i_step = 0; i_step<n_step; i_step++)
  {
    // Increment control parameters
    Parameters::P_mag+=dp_mag;
    Parameters::T_mag+=dt_mag;

    // Solve
    problem.newton_solve();

    // Document
    problem.doc_solution();
  }
   
} //End of main
