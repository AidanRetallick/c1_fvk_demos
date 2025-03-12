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

// Strings
#include <string.h>
// Floating point environment
#include <fenv.h>
// Directory parsing
#include <dirent.h>

// Generic oomph-lib routines
#include "generic.h"
// The equations
#include "c1_foeppl_von_karman.h"
// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;


//===========================================================================
/// Namespace for problem parameters
//===========================================================================
namespace Parameters
{
  /// Square edge length
  double L = 1.0;

  /// Plate vertices.
  ///     L/2*{(1,1), (-1,1), (-1,-1), (1,-1)}
  Vector<Vector<double>> Vertices =
  {
    { L/2.0, 0.0   },
    { L/2.0, L/2.0 },
    {   0.0, L/2.0 },
    {   0.0,   0.0 }
  };

  /// The plate thickness
  double Thickness = 0.01;

  /// Poisson ratio
  double Nu = 0.5;

  /// FvK parameter (slightly dangerous; should really
  /// be computed as a dependent parameter to accommodate
  /// changes in Nu and Thickness.
  double Eta = 0.0;

  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Magnitude of shear stress
  double T_mag = 0.0;

  /// Element size
  double Element_area=0.01;


  //------- Forcing functions --------------------------------------------------
  /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    pressure = P_mag;
  }

  /// Shear stress depending on the position (x,y)
  void get_in_plane_force(const Vector<double>& x, Vector<double>& tau)
  {
    tau[0]= T_mag;
    tau[1]= T_mag;
  }


  //---------------------------------------------------------------------
  // Eigenmode data and functions

  /// Characters to store boundary conditions on straight edges:
  /// f: free
  /// p: resting pin
  /// s: sliding clamp
  /// c: true clamp
  /// n: none set yet
  char Bc_char[3] = "nn";

  /// Eigenvalue
  complex<double> Lambda = 0.0+0.0i;

  /// Eigenvector
  Vector<complex<double>> B(4, 0.0+0.0i);

  /// Theta dependent part of the mode
  complex<double> F(const double& theta)
  {
    return
      B[0] * sin((Lambda + 1.0) * theta)
      + B[1] * cos((Lambda + 1.0) * theta)
      + B[2] * sin((Lambda - 1.0) * theta)
      + B[3] * cos((Lambda - 1.0) * theta);
  }

  /// First derivative of the theta dependent part of the mode
  complex<double> dF(const double& theta)
  {
    return
      (Lambda + 1.0) * (B[0] * cos((Lambda + 1.0) * theta)
			- B[1] * sin((Lambda + 1.0) * theta))
      + (Lambda - 1.0) * (B[2] * cos((Lambda - 1.0) * theta)
			  - B[3] * sin((Lambda - 1.0) * theta));
  }

  /// Second derivative of the theta dependent part of the mode
  complex<double> d2F(const double& theta)
  {
    return
      - pow(Lambda + 1.0, 2) * (B[0] * sin((Lambda + 1.0) * theta)
				+ B[1] * cos((Lambda + 1.0) * theta))
      - pow(Lambda - 1.0, 2) * (B[2] * sin((Lambda - 1.0) * theta)
				+ B[3] * cos((Lambda - 1.0) * theta));
  }

  /// Get the analytic biharmonic eigenmode at x
  void eigenmode(const Vector<double>& x,
		 double& w)
  {
    // Polar coordinates
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double theta = atan2(x[1], x[0]);

    // Angle dependence
    complex<double> f = F(theta);

    // Set mode at x: w=r^{Lambda+1}*F(theta)
    w = real( pow(r,Lambda+1.0) * f );
  }

  /// Get the theta derivative of the analytic biharmonic eigenmode at x
  void deigenmode_dtheta(const Vector<double>& x,
			 double& w)
  {
    // Polar coordinates
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double theta = atan2(x[1], x[0]);

    // Angle dependence
    complex<double> df = dF(theta);
    // Set mode at x: w=r^{Lambda+1}*F(theta)
    w = real( pow(r,Lambda+1.0) * df );
  }

  /// Get the second theta derivative of the analytic biharmonic eigenmode at x
  void d2eigenmode_dtheta2(const Vector<double>& x,
			   double& w)
  {
    // Polar coordinates
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double theta = atan2(x[1], x[0]);

    // Angle dependence
    complex<double> d2f = d2F(theta);

    // Set mode at x: w=r^{Lambda+1}*F(theta)
    w = real( pow(r,Lambda+1.0) * d2f );
  }

  /// Get the r derivative of the analytic biharmonic eigenmode at x
  void deigenmode_dr(const Vector<double>& x,
		     double& w)
  {
    // Polar coordinates
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double theta = atan2(x[1], x[0]);

    // Angle dependence
    complex<double> f = F(theta);

    // Set mode at x: w=r^{Lambda+1}*F(theta)
    w = real( (Lambda + 1.0) * pow(r,Lambda) * f );
  }

  /// Get the second r derivative of the analytic biharmonic eigenmode at x
  void d2eigenmode_dr2(const Vector<double>& x,
		       double& w)
  {
    // Polar coordinates
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double theta = atan2(x[1], x[0]);

    // Angle dependence
    complex<double> f = F(theta);

    // Set mode at x: w=r^{Lambda+1}*F(theta)
    w = real( (Lambda + 1.0) * Lambda * pow(r,Lambda-1.0) * f );
  }

  /// Get the mixed derivative of the analytic biharmonic eigenmode at x
  void d2eigenmode_drdtheta(const Vector<double>& x,
			    double& w)
  {
    // Polar coordinates
    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
    double theta = atan2(x[1], x[0]);

    // Angle dependence
    complex<double> df = dF(theta);

    // Set mode at x: w=r^{Lambda+1}*F(theta)
    w = real( (Lambda + 1.0) * pow(r,Lambda) * df );
  }

  /// Vector wrapper to the eigenmode function to use as error function
  void exact_soln(const Vector<double>& x, Vector<double>& w)
  {
    // Fill out the displacement for the second
    eigenmode(x,w[0]);
  }

  /// Null function for any zero (homogenous) BCs
  void get_null_fct(const Vector<double>& x, double& value)
  {
    value = 0.0;
  }

  //----------------------------------------------------------------------
  // Coordinate transform Jacobian (d_polar/d_cart)
  void jac_polar_cart(const Vector<double>& x,
		      DenseMatrix<double>& jac)
  {
    double r = sqrt(x[0]*x[0]+x[1]*x[1]);

    jac(0,0) = x[0] / r;
    jac(0,1) = x[1] / r;
    jac(1,0) =-x[1] / (r*r);
    jac(1,1) = x[0] / (r*r);
  }

  /// Coordinate transform Hessian (d2_polar/d_cart2)
  void hess_polar_cart(const Vector<double>& x,
		       Vector<DenseMatrix<double>>& hess)
  {
    double r = sqrt(x[0]*x[0]+x[1]*x[1]);

    hess[0](0,0) = x[1]*x[1] / (r*r*r);
    hess[0](0,1) =-x[0]*x[1] / (r*r*r);
    hess[0](1,0) =-x[0]*x[1] / (r*r*r);
    hess[0](1,1) = x[0]*x[0] / (r*r*r);

    hess[1](0,0) = 2.0*x[0]*x[1] / (r*r*r*r);
    hess[1](0,1) = x[1]*x[1]-x[0]*x[0] / (r*r*r*r);
    hess[1](1,0) = x[1]*x[1]-x[0]*x[0] / (r*r*r*r);
    hess[1](1,1) =-2.0*x[0]*x[1] / (r*r*r*r);
  }

  /// The alpha-th first derivative of the analytic biharmonic eigenmode at x
  void deigenmode(const Vector<double>& x,
		  const unsigned& alpha,
		  double& dw)
  {
    // r and theta derivatives
    double de_dr = 0.0;
    double de_dt = 0.0;

    // Fill out the polar derivatives of the eigenmode
    deigenmode_dr(x,de_dr);
    deigenmode_dtheta(x,de_dt);

    // Get the Jacobian to convert between polar and cartesian derivatives
    DenseMatrix<double> jac(2,2,0.0);
    jac_polar_cart(x, jac);

    // Set mode at x: w=r^{Lambda+1}*F(theta)
    dw = jac(0,alpha)*de_dr + jac(1,alpha)*de_dt;
  }

  /// The alpha,beta-th second derivative of the analytic biharmonic eigenmode
  /// at x
  void d2eigenmode(const Vector<double>& x,
		   const unsigned& alpha,
		   const unsigned& beta,
		   double& d2w)
  {
    // r and theta derivatives
    double de_dr    = 0.0;
    double de_dt    = 0.0;
    double d2e_dr2  = 0.0;
    double d2e_drdt = 0.0;
    double d2e_dt2  = 0.0;

    // Fill out the polar derivatives of the eigenmode
    deigenmode_dr(x,de_dr);
    deigenmode_dtheta(x,de_dt);
    d2eigenmode_dr2(x, d2e_dr2);
    d2eigenmode_drdtheta(x, d2e_drdt);
    d2eigenmode_dtheta2(x, d2e_dt2);

    // And stick them in vectors to make looping easier
    Vector<double> d_polar = {de_dr, de_dt};
    Vector<double> d2_polar = {d2e_dr2, d2e_drdt, d2e_dt2};

    // Get the Jacobian and Hessian to convert between polar and cartesian
    // derivatives
    DenseMatrix<double>
      jac(2, 2, 0.0);
    Vector<DenseMatrix<double>> hess(2,DenseMatrix<double>(2,2,0.0));
    jac_polar_cart(x, jac);
    hess_polar_cart(x, hess);

    // Add contributions to second derivative
    d2w = 0.0;
    for(unsigned gamma = 0; gamma < 2; gamma++)
    {
      // Add contributions from gradient
      d2w += d_polar[gamma] * hess[gamma](alpha,beta);
      for(unsigned delta = 0; delta < 2; delta++)
      {
	// Add contributions from second derivatives
	d2w += d2_polar[gamma+delta] * jac(gamma,alpha) * jac(delta,beta);
      }
    }
  }

  /// Get the x derivative of the analytic biharmonic eigenmode at x
  void deigenmode_dx(const Vector<double>& x,
		     double& dw)
  {
    deigenmode(x,0,dw);
  }

  /// Get the y derivative of the analytic biharmonic eigenmode at x
  void deigenmode_dy(const Vector<double>& x,
		     double& dw)
  {
    deigenmode(x,1,dw);
  }

  /// Get the second x derivative of the analytic biharmonic eigenmode at x
  void d2eigenmode_dx2(const Vector<double>& x,
		       double& d2w)
  {
    d2eigenmode(x,0,0,d2w);
  }

  /// Get the xy derivative of the analytic biharmonic eigenmode at x
  void d2eigenmode_dxdy(const Vector<double>& x,
			double& d2w)
  {
    d2eigenmode(x,0,1,d2w);
  }

  /// Get the second y derivative of the analytic biharmonic eigenmode at x
  void d2eigenmode_dy2(const Vector<double>& x,
		       double& d2w)
  {
    d2eigenmode(x,1,1,d2w);
  }

} // end of Parameters


///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

  /// Constructor
  UnstructuredFvKProblem();

  /// Destructor
  ~UnstructuredFvKProblem()
  {
    // Close the trace file as we are done with the problem
    Trace_file.close();

    // Clean up memory
    delete Bulk_mesh_pt;
    delete Boundary_pt;
    delete Boundary0_pt;
    delete Boundary1_pt;
    delete Boundary2_pt;
    delete Boundary3_pt;
  };


  /// Print information about the parameters we are trying to solve for.
  void actions_before_newton_solve()
  {
    oomph_info << "-------------------------------------------------------"
    << std::endl;
    oomph_info << "Solving for P = " << Parameters::P_mag << std::endl;
    oomph_info << "Solving for T = " << Parameters::T_mag << std::endl;
    oomph_info << "         step = " << Doc_info.number() << std::endl;
    oomph_info << "-------------------------------------------------------"
    << std::endl;
  }

  /// Update after solve (empty)
  void actions_after_newton_solve() {}


  /// Update the boundary conditions along with the constraints and eqn
  /// numbering
  void update_boundary_conditions()
  {
    // Reset by unpinning all the nodal dofs in the mesh
    unsigned n_node = Bulk_mesh_pt->nnode();
    for(unsigned i_node = 0; i_node < n_node; i_node++)
    {
      Bulk_mesh_pt->node_pt(i_node)->unpin_all();
    }

    // Apply the new boundary conditions
    apply_boundary_conditions();

    // Pin in-plane displacements throughout the bulk
    if(Solve_linear_bending)
    {
      pin_all_in_plane_displacements();
    }

    // Update the equation numbering
    assign_eqn_numbers();

  //   // Output the new dofs
  //   describe_dofs();
  }


  /// Doc the solution
  void doc_solution(const std::string& comment="");

  // [zdec] This doesn't work, dynamic cast always fails -- returns 0
  // /// Overloaded version of the problem's access function to
  // /// the mesh. Recasts the pointer to the base Mesh object to
  // /// the actual mesh type.
  // TriangleMesh<ELEMENT>* mesh_pt()
  // {
  //   oomph_info << Problem::mesh_pt() << std::endl;
  //   oomph_info << dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt()) << std::endl;
  //   return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt());
  // }

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return Bulk_mesh_pt;
  }

private:

  /// Setup and build the mesh
  void build_mesh();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

  /// Pin all in-plane displacements in the domain
  /// (for solving the linear problem)
  void pin_all_in_plane_displacements();

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// Triangle mesh parameters
  TriangleMeshParameters* Triangle_mesh_parameters_pt;

  /// Closed outer boundary
  TriangleMeshClosedCurve* Boundary_pt;

  /// Polyline defining boundary 0
  TriangleMeshPolyLine* Boundary0_pt;

  /// Polyline defining boundary 1
  TriangleMeshPolyLine* Boundary1_pt;

  /// Polyline defining boundary 2
  TriangleMeshPolyLine* Boundary2_pt;

  /// Polyline defining boundary 3
  TriangleMeshPolyLine* Boundary3_pt;

  /// Doc info object for labeling output
  DocInfo Doc_info;

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Flag to turn the elements into linear bending (biharmonic)
  bool Solve_linear_bending;

  /// Keep track of boundary ids, (b)ottom, (r)ight, (t)op, (l)eft
  // (slightly redundant in this example)
  // ((after rotation this naming convention is unhelpful))
  enum
  {
    Boundary_b_bnum = 0,
    Boundary_r_bnum = 1,
    Boundary_t_bnum = 2,
    Boundary_l_bnum = 3
  };

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to "surface" mesh
  Mesh* Constraint_mesh_pt;

  /// Pointer to element that contains the central point
  GeomObject* Central_element_geom_obj_pt;

  /// Local coordinate in element pointed to by Central_element_pt
  /// that contains central point
  Vector<double> Central_point_local_coord;

}; // end_of_problem_class


//==start_of_problem_constructor===========================================
/// Constructor
//=========================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem()
{

  // Set output directory
  Doc_info.set_directory("RESLT");

  // Step number
  Doc_info.number()=0;

  // Build the mesh
  build_mesh();

  // Complete problem setup
  complete_problem_setup();


  // Output parameters
  oomph_info << "Problem parameters:\n"
  << "L            " << Parameters::L         << std::endl
  << "thickness    " << Parameters::Thickness << std::endl
  << "nu           " << Parameters::Nu        << std::endl
  << "eta          " << Parameters::Eta       << std::endl
  << "Element area " << Parameters::Element_area << std::endl;


  // Open trace file
  char filename[100];
  strcpy(filename, (Doc_info.directory()+"/trace.dat").c_str());
  Trace_file.open(filename);


  // Assign equation numbers
  oomph_info << "Number of equations: "
  << assign_eqn_numbers() << '\n';

} // end Constructor


//==start_of_build_mesh====================================================
/// Build the rectangular mesh
//=========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
  // *********************************************************************
  // ********** THIS DIAGRAM IS WRONG NOW THAT WE USE L & Alpha **********
  // *********************************************************************
  //=================================//
  //  Rotated square mesh boundary   //
  //                                 //
  //                X (0,A)          //
  //              /   \              //
  //         e1 /       \ e0         //
  //          /           \          //
  // -(A,0) X               X (A,0)  //
  //          \           /          //
  //         e2 \       / e3         //
  //              \   /              //
  //                X-(0,A)          //
  //                                 //
  //                                 //
  //=================================//


  // Get the vertices from the parameters
  Vector<double>
    vertex0 = Parameters::Vertices[0],
    vertex1 = Parameters::Vertices[1],
    vertex2 = Parameters::Vertices[2],
    vertex3 = Parameters::Vertices[3];

  // Declare the edges...
  Vector<Vector<double>>
    edge0(2,Vector<double>(2,0.0)),
    edge1(2,Vector<double>(2,0.0)),
    edge2(2,Vector<double>(2,0.0)),
    edge3(2,Vector<double>(2,0.0));

  // ...and assign their endpoints
  edge0[0] = vertex0;
  edge0[1] = vertex1;
  edge1[0] = vertex1;
  edge1[1] = vertex2;
  edge2[0] = vertex2;
  edge2[1] = vertex3;
  edge3[0] = vertex3;
  edge3[1] = vertex0;

  // Define boundaries from edges
  Boundary0_pt = new TriangleMeshPolyLine(edge0, 0);
  Boundary1_pt = new TriangleMeshPolyLine(edge1, 1);
  Boundary2_pt = new TriangleMeshPolyLine(edge2, 2);
  Boundary3_pt = new TriangleMeshPolyLine(edge3, 3);

  // Create closed outer boundary
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
  boundary_polyline_pt[0] = Boundary0_pt;
  boundary_polyline_pt[1] = Boundary1_pt;
  boundary_polyline_pt[2] = Boundary2_pt;
  boundary_polyline_pt[3] = Boundary3_pt;
  Boundary_pt = new TriangleMeshClosedCurve(boundary_polyline_pt);


  // Define mesh parameters
  TriangleMeshParameters Triangle_mesh_parameters(Boundary_pt);

  // Set the maximum element area
  Triangle_mesh_parameters.element_area()=Parameters::Element_area ;

  // Build  bulk mesh
  Bulk_mesh_pt=new TriangleMesh<ELEMENT>(Triangle_mesh_parameters);

  //Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);

  // Combine submeshes into a single Mesh (bit over the top here; could
  // have assigned bulk mesh to mesh_pt() directly).
  build_global_mesh();

}// end build_mesh



//==start_of_complete======================================================
/// Set boundary conditions and complete the build of
/// all elements
//=========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{

  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the pressure & temperature function pointers and the physical constants
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->in_plane_forcing_fct_pt() = &Parameters::get_in_plane_force;

    // Assign the parameter function pointers for the element
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_pt() = &Parameters::Eta;
  }
  // Set the boundary conditions
  apply_boundary_conditions();

} // end of complete




//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
  //----------------------------------------------------------------------------
  // Sets of possible boundary conditions

  // The free boundary condition is completely unpinned
  static const Vector<Vector<unsigned>> free{{},{}};

  // Out-of-plane dofs:
  // |  0  |  1  |  2  |  3  |  4  |  5  |
  // |  w  | w_n | w_t | w_nn| w_nt| w_tt|
  // Possible boundary conditions for the out-of-plane displacement
  static const Vector<Vector<unsigned>> resting_pin_dofs{{0, 1, 3},
							 {0, 2, 5}};
  static const Vector<Vector<unsigned>> sliding_clamp_dofs{{2, 4},
							   {1, 4}};
  static const Vector<Vector<unsigned>> true_clamp_dofs{{0, 1, 2, 3, 4},
                                                        {0, 1, 2, 4, 5}};

  //----------------------------------------------------------------------------
  // Storage for boundary conditions to each straight edge
  Vector<Vector<unsigned>> straight_edge_pinned_w_dofs(2, Vector<unsigned>{});
  for(unsigned i_straight_b = 0; i_straight_b < 2; i_straight_b++)
  {
    switch(Parameters::Bc_char[i_straight_b])
    {
    case 'f':
      cout << "Edge " << i_straight_b << " is f" << std::endl;
      straight_edge_pinned_w_dofs[i_straight_b] = free[i_straight_b];
      break;
    case 'p':
      cout << "Edge " << i_straight_b << " is p" << std::endl;
      straight_edge_pinned_w_dofs[i_straight_b] = resting_pin_dofs[i_straight_b];
      break;
    case 's':
      cout << "Edge " << i_straight_b << " is s" << std::endl;
      straight_edge_pinned_w_dofs[i_straight_b] = sliding_clamp_dofs[i_straight_b];
      break;
    case 'c':
      cout << "Edge " << i_straight_b << " is c" << std::endl;
      straight_edge_pinned_w_dofs[i_straight_b] = true_clamp_dofs[i_straight_b];
      break;
    }
  }

  //----------------------------------------------------------------------
  // Loop over the real straight corner boundaries
  for (unsigned i_bound = 0; i_bound < 2; i_bound++)
  {
    unsigned bound;
    switch(i_bound)
    {
    case 0:
      bound = 3;
      break;
    case 1:
      bound = 2;
      break;
    }
    // Loop over straight side elements and apply homogenous BCs
    unsigned n_b_element = Bulk_mesh_pt->nboundary_element(bound);
    unsigned n_pinned_w_dofs = straight_edge_pinned_w_dofs[i_bound].size();
    for(unsigned e=0;e<n_b_element;e++)
    {
      // Get pointer to bulk element adjacent to b
      ELEMENT* el_pt =
	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(bound,e));

      // Pin out-of-plane dofs
      for(unsigned j_dof = 0; j_dof < n_pinned_w_dofs; j_dof++)
      {
	unsigned dof = straight_edge_pinned_w_dofs[i_bound][j_dof];
	cout << "On " << bound << " pinning " << dof << endl;
	el_pt->fix_out_of_plane_displacement_dof(dof,
						 bound,
						 Parameters::get_null_fct);
      } // End loop over out-of-plane dofs [j_dof]
    } // End loop over boundary elements [e]
  } // End loop over boundaries [i_bound]

  //----------------------------------------------------------------------
  // Loop over the two "internal" edges
  for (unsigned i_bound = 0; i_bound < 2; i_bound++)
  {
    unsigned n_b_element = Bulk_mesh_pt->nboundary_element(i_bound);
    for (unsigned e = 0; e < n_b_element; e++)
    {
      // Get pointer to bulk element adjacent to curved arc
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
	Bulk_mesh_pt->boundary_element_pt(i_bound, e));

      // // Pin in-plane dofs
      // for(unsigned i=0; i<n_pinned_u_dofs; i++)
      // {
      //   unsigned idof=circular_arc_pinned_u_dofs[i];
      //   el_pt->fix_in_plane_displacement_dof(idof,
      // 					   Circular_arc_bnum,
      // 					   Parameters::get_null_fct);
      // }
      // Pin out-of-plane dofs (resting pin -- only set tangent derivatives)
      for (unsigned idof = 0; idof < 6; ++idof)
      {
        switch (idof)
        {
            // [hierher] Make function of arclength rather than global x
          case 0:
            el_pt->fix_out_of_plane_displacement_dof(
              idof, i_bound, Parameters::eigenmode);
            break;
          case 1:
            el_pt->fix_out_of_plane_displacement_dof(
              idof, i_bound, Parameters::deigenmode_dx);
            break;
          case 2:
            el_pt->fix_out_of_plane_displacement_dof(
              idof, i_bound, Parameters::deigenmode_dy);
            break;
          case 3:
	    el_pt->fix_out_of_plane_displacement_dof(
	      idof, i_bound, Parameters::d2eigenmode_dx2);
           	break;
          case 4:
            el_pt->fix_out_of_plane_displacement_dof(
              idof, i_bound, Parameters::d2eigenmode_dxdy);
            break;
          case 5:
            el_pt->fix_out_of_plane_displacement_dof(
              idof, i_bound, Parameters::d2eigenmode_dy2);
            break;
          default:
            // Leave free
            break;
        } // End of switch-case [idof]
      } // End loop over dofs [idof]
    } // End loop over boundary elements [e]
  } // End loop over internal boundaries

} // end set bc



//==start_of_pin_all_in_plane_displacements=====================================
/// Pin the in-plane displacements
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_all_in_plane_displacements()
{
  unsigned nnode = Bulk_mesh_pt->nnode();
  for(unsigned inode=0; inode<nnode; inode++)
  {
    Bulk_mesh_pt->node_pt(inode)->pin(0);
    Bulk_mesh_pt->node_pt(inode)->set_value(0,0.0);
    Bulk_mesh_pt->node_pt(inode)->pin(1);
    Bulk_mesh_pt->node_pt(inode)->set_value(1,0.0);
  }
}



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const
						   std::string& comment)
{
  ofstream some_file;
  char filename[100];


  // Number of plot points for fine (full) output. Lots of plot points so
  // we see the goodness of the high-order Hermite/Bell
  // interpolation.
  unsigned npts = 10;
  sprintf(filename, "%s/soln_%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
  << comment << "\"\n";
  some_file.close();


  // interpolated_foeppl_von_karman_values(s)
  // and tidy it up; it's a total mess!!!
  //
  // However, it appears to return:
  //
  // interpolated_f[0] = interpolated_w[0]; // w
  // interpolated_f[1] = interpolated_dwdxi(0, 0); // dwdx
  // interpolated_f[2] = interpolated_dwdxi(0, 1); // dwdy
  // interpolated_f[3] = interpolated_d2wdxi2(0, 0); // d2wdx2
  // interpolated_f[4] = interpolated_d2wdxi2(0, 1); // d2wdxdy
  // interpolated_f[5] = interpolated_d2wdxi2(0, 2); // d2wdy2
  // interpolated_f[6] = interpolated_u[0]; // ux
  // interpolated_f[7] = interpolated_u[1]; // uy
  // interpolated_f[8] = interpolated_dudxi(0, 0); // duxdx
  // interpolated_f[9] = interpolated_dudxi(0, 1); // duxdy
  // interpolated_f[10] = interpolated_dudxi(1, 0); // duydx
  // interpolated_f[11] = interpolated_dudxi(1, 1); // duydy
  //
  // Make this clear in the doxygen-documented comments on top
  // of the function.

  // Doc error and return of the square of the L2 error
  //---------------------------------------------------
  //double error,norm,dummy_error,zero_norm;
  double error, zero_norm;
  sprintf(filename,"RESLT/error%i.csv",Doc_info.number());
  some_file.open(filename);

  Bulk_mesh_pt->compute_error(some_file,
			      Parameters::exact_soln,
			      error,
			      zero_norm);
  some_file.close();

  // Doc L2 error and norm of solution
  oomph_info << "Absolute norm of computed solution: " << sqrt(error)
	     << std::endl;

  oomph_info << "Norm of computed solution: " << sqrt(zero_norm)
	     << std::endl;

  // Find the solution at r=0
  //   // ----------------------
  MeshAsGeomObject* Mesh_as_geom_obj_pt=
    new MeshAsGeomObject(Bulk_mesh_pt);
  Vector<double> s(2);
  GeomObject* geom_obj_pt=0;
  Vector<double> r(2,0.0);
  Mesh_as_geom_obj_pt->locate_zeta(r,geom_obj_pt,s);
  // Compute the interpolated displacement vector
  Vector<double> u_0(12,0.0);
  u_0=dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_foeppl_von_karman(s);

  oomph_info << "w in the middle: " <<std::setprecision(15) << u_0[0] << std::endl;

  Trace_file << u_0[0] << '\n';

  // Doc error and return of the square of the L2 error
  //---------------------------------------------------
  sprintf(filename,"RESLT/L2-norm%i.dat",
	  Doc_info.number());
  some_file.open(filename);

  some_file<<"### L2 Norm\n";
  some_file<<"##  Format: err^2 norm^2 \n";
  // Print error in prescribed format
  some_file<< error <<" "<< zero_norm <<"\n";
  some_file.close();

  // Bump
  Doc_info.number()++;

  // Clean up
  delete Mesh_as_geom_obj_pt;

} // end of doc



//=======start_of_main========================================
/// Driver code for demo of unstructured C1 Foeppl-von Karman
/// elements
//============================================================
int main(int argc, char **argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

 // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  // Define possible command line arguments and parse the ones that
  // were actually specified
  // Directory for solution
  string output_dir="RESLT";
  CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

  // Poisson Ratio
  CommandLineArgs::specify_command_line_flag("--nu", &Parameters::Nu);

  // Applied Pressure
  CommandLineArgs::specify_command_line_flag("--eta", &Parameters::Eta);

  // Element Area
  CommandLineArgs::specify_command_line_flag("--element_area",
					     &Parameters::Element_area);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>>
    problem;

  problem.newton_solver_tolerance() = 1.0e-10;

  // Pre mess doc
  problem.doc_solution();

  // The name of the directory we expect to find our modes
  char dirname[100];
  sprintf(dirname, "Modes_square");
  DIR* dir = opendir(dirname);

  oomph_info << "About to start the loop" << std::endl;

  // Loop over the eigenmodes in the Modes directory and validate that we are
  // close to them.
  for(dirent* dent; (dent = readdir(dir)) != NULL; )
  {
    // Get this entries name and take the first and third letters to be the bcs
    // on the first and second straight edges respectively
    std::string filename = dent->d_name;
    // Skip the directory links
    if(filename=="." || filename=="..")
    {
      continue;
    }
    Parameters::Bc_char[0]=filename[0];
    Parameters::Bc_char[1]=filename[2];
    oomph_info << Parameters::Bc_char << std::endl;

    // Open the file and copy the first four lines into the eigenvector
    std::ifstream fin((string)(dirname)+'/'+filename);
    std::string line;
    for(unsigned i=0; i<4; i++)
    {
      // Get a line (complex number)
      std::getline(fin, line);
      // Copy the line to the complex evec entry via a stringstream
      stringstream linestream('('+line+')');
      oomph_info << linestream.str() << std::endl;//Parameters::Lambda;
      linestream >> Parameters::B[i];
    }
    // The last line is the eigenvalue
    std::getline(fin, line);
    // Copy the line to the complex eval entry via a stringstream
    stringstream linestream('('+line+')');
    oomph_info << linestream.str() << std::endl;
    linestream >> Parameters::Lambda;

    // Apply the appropriate boundary conditions for this test
    problem.update_boundary_conditions();
    // Now solve and doc the solution with these bcs
    problem.steady_newton_solve();
    problem.doc_solution();
  }

  // Print success
  oomph_info<<"Exiting Normally\n";
} //End of main
