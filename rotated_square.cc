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

  /// Angle of rotation (radians obviously)
  double Alpha = Pi/5.0;
  
  /// Plate vertices. For the unrotated square, these coincide with:
  ///     L/2*{(1,1), (-1,1), (-1,-1), (1,-1)}
  Vector<Vector<double>> Vertices =
  {
    { L/2.0*cos(Alpha)-L/2.0*sin(Alpha),  L/2.0*sin(Alpha)+L/2.0*cos(Alpha)},
    {-L/2.0*cos(Alpha)-L/2.0*sin(Alpha), -L/2.0*sin(Alpha)+L/2.0*cos(Alpha)},
    {-L/2.0*cos(Alpha)+L/2.0*sin(Alpha), -L/2.0*sin(Alpha)-L/2.0*cos(Alpha)},
    { L/2.0*cos(Alpha)+L/2.0*sin(Alpha),  L/2.0*sin(Alpha)-L/2.0*cos(Alpha)}
  };
 
  /// The plate thickness
  double Thickness = 0.01;
 
  /// Poisson ratio
  double Nu = 0.5;

  /// FvK parameter (slightly dangerous; should really
  /// be computed as a dependent parameter to accommodate
  /// changes in Nu and Thickness.
  double Eta = 12.0*(1.0-Nu*Nu)*Thickness*Thickness;
 
  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Magnitude of shear stress
  double T_mag = 0.0;
 
  /// Element size
  double Element_area=0.01;

  
  // ---- Parametric boundaries ------------------------------------------------

  /// The type of function pointer that set_up_rotated_dofs() expects
  typedef void Norm_and_tan_func(const Vector<double>&,
				 Vector<double>&,
				 Vector<double>&,
				 DenseMatrix<double>&,
				 DenseMatrix<double>&);
  
  /// The normal and tangential directions for a straight line. We need the
  /// derivatives so we can form the Hessian and the Jacobian of the rotation
  void get_normal_and_tangent_straight_boundary_for_straight_line(const Vector<double>& x0,
								  const Vector<double>& x1,
								  const Vector<double>& x,
								  Vector<double>& n,
								  Vector<double>& t,
								  DenseMatrix<double>& dn,
								  DenseMatrix<double>& dt)
  {
    double dx = x1[0] - x0[0];
    double dy = x1[1] - x0[1];
    double mag = sqrt(dx*dx + dy*dy);
    
    // // Flip the normal and tangent if we are below the centreline 
    // if(x[0]+x[1]<7.0*A/8.0)
    //   {mag*=-1.0;}
    
    // Fill in the normal and derivatives of the normal
    n[0] =-dy/mag;
    n[1] = dx/mag;
    
    // Fill in the tangent and derivatives of the tangent
    t[0] = dx/mag;
    t[1] = dy/mag;
    
    // Zero derivatives for straight lines
    dn(0,0) = 0.0;
    dn(1,0) = 0.0;
    dn(0,1) = 0.0;
    dn(1,1) = 0.0;
    
    dt = dn;
  }
  
  //--------------------------------------------------------------------------
  // Wrappers to the generic straight line boundary rotation

  /// The normal and tangential directions for boundary 0. We need the
  /// derivatives so we can form the Hessian and the Jacobian of the rotation
  void get_normal_and_tangent_straight_boundary_0(const Vector<double>& x,
						  Vector<double>& n,
						  Vector<double>& t,
						  DenseMatrix<double>& dn,
						  DenseMatrix<double>& dt)
  {
    // Get the endpoints of boundary 0
    const Vector<double> x0 = Vertices[0];
    const Vector<double> x1 = Vertices[1];

    // Fill in the normal and tangent vectors using the generic function for
    // straight lines
    get_normal_and_tangent_straight_boundary_for_straight_line(x0, x1, x,
							       n, t, dn, dt);
  }
  
  /// The normal and tangential directions for boundary 1. We need the
  /// derivatives so we can form the Hessian and the Jacobian of the rotation
  void get_normal_and_tangent_straight_boundary_1(const Vector<double>& x,
						  Vector<double>& n,
						  Vector<double>& t,
						  DenseMatrix<double>& dn,
						  DenseMatrix<double>& dt)
  {
    // Get the endpoints of boundary 1
    const Vector<double> x0 = Vertices[1];
    const Vector<double> x1 = Vertices[2];

    // Fill in the normal and tangent vectors using the generic function for
    // straight lines
    get_normal_and_tangent_straight_boundary_for_straight_line(x0, x1, x,
							       n, t, dn, dt);
  }
  
  /// The normal and tangential directions for boundary 2. We need the
  /// derivatives so we can form the Hessian and the Jacobian of the rotation
  void get_normal_and_tangent_straight_boundary_2(const Vector<double>& x,
						  Vector<double>& n,
						  Vector<double>& t,
						  DenseMatrix<double>& dn,
						  DenseMatrix<double>& dt)
  {
    // Get the endpoints of boundary 2
    const Vector<double> x0 = Vertices[2];
    const Vector<double> x1 = Vertices[3];

    // Fill in the normal and tangent vectors using the generic function for
    // straight lines
    get_normal_and_tangent_straight_boundary_for_straight_line(x0, x1, x,
							       n, t, dn, dt);
  }
  
  /// The normal and tangential directions for boundary 3. We need the
  /// derivatives so we can form the Hessian and the Jacobian of the rotation
  void get_normal_and_tangent_straight_boundary_3(const Vector<double>& x,
						  Vector<double>& n,
						  Vector<double>& t,
						  DenseMatrix<double>& dn,
						  DenseMatrix<double>& dt)
  {
    // Get the endpoints of boundary 3
    const Vector<double> x0 = Vertices[3];
    const Vector<double> x1 = Vertices[0];

    // Fill in the normal and tangent vectors using the generic function for
    // straight lines
    get_normal_and_tangent_straight_boundary_for_straight_line(x0, x1, x,
							       n, t, dn, dt);
  }
  

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


  //-------- Boundary conditions -----------------------------------------------
  /// Helper function to specify boundary conditions (here all homogeneous)
  /// as a function of both coordinates (This is convenient for this problem;
  /// other interfaces that specify boundary conditions in terms of
  /// boundary coordinate exist).
  void get_null_fct(const Vector<double>& x, double& value)
  {
    value = 0.0;
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
  
  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in Parameters
  void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);
  
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

  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom(Bulk_mesh_pt);
  
  // Complete problem setup
  complete_problem_setup();


  // Output parameters
  oomph_info << "Problem parameters:\n"
  << "L            " << Parameters::L         << std::endl
  << "Alpha        " << Parameters::Alpha     << std::endl
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




//==============================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
/// begin rotate_edge_degrees_of_freedom
//==============================================================================
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
  // Number of boundaries in our domain
  unsigned n_boundaries = 4;

  // Store the edge parametrisations in a vector
  Vector<Parameters::Norm_and_tan_func*> boundary_parametrisation_pt(n_boundaries);
  boundary_parametrisation_pt[Boundary_b_bnum] =
    &Parameters::get_normal_and_tangent_straight_boundary_0;
  boundary_parametrisation_pt[Boundary_r_bnum] =
    &Parameters::get_normal_and_tangent_straight_boundary_1;
  boundary_parametrisation_pt[Boundary_t_bnum] =
    &Parameters::get_normal_and_tangent_straight_boundary_2;
  boundary_parametrisation_pt[Boundary_l_bnum] =
    &Parameters::get_normal_and_tangent_straight_boundary_3;

  // Loop over the boundaries
  for(unsigned b=0; b<n_boundaries; b++)
  {
    // Loop over the bulk elements
    unsigned n_element = bulk_mesh_pt-> nelement();
    for(unsigned e=0; e<n_element; e++)
    {
      // Get pointer to bulk element
      ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

      // Calculate which nodes are on boundary b
      const unsigned nnode=3;

      // Count the number of boundary nodes on boundary b
      Vector<unsigned> boundary_nodes;
      for (unsigned n=0; n<nnode;++n)
      {
	// Rotate nodes if on boundary b
	if (el_pt->node_pt(n)->is_on_boundary(b))
	{ boundary_nodes.push_back(n); }
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if(!boundary_nodes.empty())
      {
	// Rotate the nodes by passing the index of the nodes and the
	// normal / tangent vectors to the element
	el_pt->set_up_rotated_dofs(boundary_nodes.size(),
				   boundary_nodes,
				   boundary_parametrisation_pt[b]);
      }
    } // End loop over elements [e]
  } // End loop over boundaries [b]
}// end rotate_edge_degrees_of_freedom





//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{


  //------------------------------------------------------------------
  //------------------------------------------------------------------
  // Boundary conditions for FvK elements are complicated and we
  // provide an illustration of how to apply all physically relevant
  // boundary conditions for problems with axis-aligned boundaries
  // here. Other tutorials/driver codes explain what to do about
  // (i) straight boundaries that are not aligned with the coordinate axes
  // and (ii) curvilinear boundaries.
  //
  // FvK elements have two different types of degrees of freedom:
  // (1) The ones associated with the two in-plane displacements, u_x and u_y.
  //     These are interpolated with standard C0 continuous Lagrange
  //     interpolation between all the nodes in the underlying
  //     TElement<2,NNODE_1D>. Each node stores the values
  //     of the two in-plane displacements. We enumerate these dofs
  //     as 0 for the x displacement, and 1 for the y displacement.
  // (2) The dofs associated with the out-of-plane displacement, w.
  //     This is interpolated with Bell (Hermite) interpolants
  //     which involve six different types of degree of freedom,
  //     enumerated from 0 to 5 in the order w, w_n, w_t, w_nn, w_nt,
  //     w_tt. These values are only stored at three vertices of the element.
  //
  // Given that the book-keeping for which node stores which type of
  // degree of freedom is complicated, we let the element do the
  // relevant assigments for us. For this purpose the FvK elements provide
  // two member functions:
  //
  //     fix_in_plane_displacement_dof(idof, b, fct_pt)
  //
  // assigns boundary conditions for the in-plane displacements via
  // the specification of
  //
  // idof   : the enumeration of the dof in the scheme listed above,
  //          so idof can take values 0 or 1.
  // b      : the mesh boundary along which the boundary condition is
  //          to be applied
  // fct_pt : a function pointer to a global function with arguments
  //          (const Vector<double> x, double& value) which computes
  //          the value for the relevant in-plane displacement as a
  //          function of the coordinate, x, a 2D vector.
  //
  // So, if the function fix_in_plane_displacement_dof(idof, b, fct_pt) is called with idof=1
  // and b=3, say, the y-in-plane displacement is pinned for all the element's
  // nodes (if any) that are located on mesh boundary 3. The value of the y-in-plane 
  // displacement is set to whatever the function pointed to by fct_pt computes when
  // evaluated at the nodal coordinate. 
  //
  // Similarly,
  //
  //      fix_out_of_plane_displacement_dof(idof, b, fct_pt);
  //
  // hierher complete once Aidan has signed off the explanation above.
  // [zdec] "This is all good." -- Aidan
  //
  //
  // Using the conventions introduced above, the following vectors identify
  // the in-plane and out-of-plane degrees of freedom to be pinned for
  // various physically meaningful boundary conditions:

 
  // [zdec] NOTE THAT THE IN-PLANE DISPLACEMENTS ARENT ROTATED
  // In-plane dofs:
  //---------------
  // |  0  |  1  |
  // |  ux |  uy |


  // Case: Pin x in-plane displacement only:
  const Vector<unsigned> pin_ux_pinned_dof{0};

  // Case: Pin y in-plane displacement only:
  const Vector<unsigned> pin_uy_pinned_dof{1};

  // Case: Pin both in-plane displacements
  const Vector<unsigned> pin_ux_and_uy_pinned_dof{0,1};


 
  // Out-of-plane dofs:
  //-------------------
  // |  0  |  1  |  2  |  3  |  4  |  5  |
  // |  w  | w_n | w_t | w_nn| w_nt| w_tt|

 
  // Possible boundary conditions for out-of-plane displacements:
  // Given that the out-of-plane displacements feature in the fourth-order
  // biharmonic operator, we can apply boundary conditions on w and
  // dw/dn, where n is the coordinate direction normal to the (assumed to be
  // axis aligned!) boundary. However if w is given along the entire
  // boundary (parametrised by the tangential coordinate, t) we also know what
  // dw/dt and d^2w/dt^2 are. Similarly, if we are given dw/dn along the 
  // boundary, (parametrised by the tangential coordinate, t) we also know 
  // what d^2w/dndt is. In the various cases below we identify physical
  // scenarios of a pinned edge (w given, dw/dn left free); a vertically
  // sliding edge (w left free; dw/dn given) and fully clamped (w and dw/dn
  // given). The four possible combinations of boundary condition are:
  //   fully free    -- nothing given,  
  //   pinned edge   -- only w(t) given,
  //   sliding clamp -- only dwdt(t) given,
  //   fully clamped -- both w(t) and dwdt(t) given.

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dt and d^2w/dt^2
  const Vector<unsigned> pinned_edge_pinned_dof{0,2,5};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) dw/dn and d^2w/dndt
  const Vector<unsigned> sliding_clamp_pinned_dof{1,4};
 
  // Case: The plate is clamped (w given, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dn,
  // dw/dt, d^2w/dndt and d^2w/dt^2
  const Vector<unsigned> fully_clamped_pinned_dof{0,1,2,4,5};


  //------------------------------------------------------------------
  //------------------------------------------------------------------

 
  // Vectors to store which boundary conditions we are applying to each edge.
  Vector<Vector<unsigned>> pinned_u_dofs(4);
  Vector<Vector<unsigned>> pinned_w_dofs(4);
 
  // Pin both in-plane displacements everywhere
  pinned_u_dofs[0] = pin_ux_and_uy_pinned_dof;
  pinned_u_dofs[1] = pin_ux_and_uy_pinned_dof;
  pinned_u_dofs[2] = pin_ux_and_uy_pinned_dof;
  pinned_u_dofs[3] = pin_ux_and_uy_pinned_dof;

  // Use pinned edge boundary conditions for the out-of-plane
  // displacements. Boundaries 0 and 2 are are constant y,
  // boundaries 1 and 3 are constant x:
  pinned_w_dofs[0] = pinned_edge_yn_pinned_dof;
  pinned_w_dofs[1] = pinned_edge_xn_pinned_dof;
  pinned_w_dofs[2] = pinned_edge_yn_pinned_dof;
  pinned_w_dofs[3] = pinned_edge_xn_pinned_dof;

 
  // Loop over all the boundaries in our bulk mesh
  unsigned n_bound = Bulk_mesh_pt->nboundary();
  for(unsigned b=0;b<n_bound;b++)
  {
    // Number of elements on b
    const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
   
    // Number of dofs we are pinning on boundary b
    const unsigned n_pinned_u_dofs = pinned_u_dofs[b].size();
    const unsigned n_pinned_w_dofs = pinned_w_dofs[b].size();

    // Loop over the elements on boundary b
    for(unsigned e=0; e<nb_element; e++)
    {
      // Get pointer to bulk element adjacent to b
      ELEMENT* el_pt =
	dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

      // Pin in-plane dofs (enumerated as explained above) for
      // all nodes on boundary b. Here we're applying homogeneous
      // BCs so all pinned values are simply set to zero.
      for(unsigned i=0; i<n_pinned_u_dofs; i++)
      {
	unsigned idof_to_be_pinned=pinned_u_dofs[b][i];
	el_pt->fix_in_plane_displacement_dof(idof_to_be_pinned, b,
					     Parameters::get_null_fct);
      }
      // Pin out-of-plane dofs (enumerated as explained above) for all
      // nodes on boundary b. Here we're applying homogeneous BCs so
      // all pinned values are simply set to zero.
      for(unsigned i=0; i<n_pinned_w_dofs; i++)
      {
	unsigned idof_to_be_pinned=pinned_w_dofs[b][i];
	el_pt->fix_out_of_plane_displacement_dof(idof_to_be_pinned, b,
						 Parameters::get_null_fct);
      }
    } // end for loop over elements on b
  } // end for loop over boundaries

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const
						   std::string& comment)
{
  ofstream some_file;
  char filename[100];
 

  // Number of plot points for coarse output (just showing the
  // element outline)
  unsigned npts = 2;
  sprintf(filename, "%s/coarse_soln_%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	    << comment << "\"\n";
  some_file.close();
 
  // Number of plot points for fine (full) output. Lots of plot points so
  // we see the goodness of the high-order Hermite/Bell
  // interpolation. 
  npts = 10;
  sprintf(filename, "%s/soln_%i.dat",
	  Doc_info.directory().c_str(),
	  Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
  << comment << "\"\n";
  some_file.close();


 
  // Write the load magnitude and centrepoint displacements to trace file
  Vector<double> w_centre(12,0.0);

  // hierher: Aidan: rename this function to
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
 
  w_centre = dynamic_cast<ELEMENT*>(Central_element_geom_obj_pt)
    ->interpolated_u_foeppl_von_karman(Central_point_local_coord);
   

  Trace_file << Parameters::P_mag  << " "
  << Parameters::T_mag  << " "
	     << w_centre[0]        << " " 
	     << w_centre[6]        << " " 
	     << w_centre[7]        << " "
  << Doc_info.number()  << endl;

  // Bump
  Doc_info.number()++;

} // end of doc



//=======start_of_main========================================
/// Driver code for demo of unstructured C1 Foeppl-von Karman
/// elements
//============================================================
int main(int argc, char **argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
 
  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>>
    problem;  

  // Set pressure
  Parameters::P_mag=10.0;
 
  // Solve the system
  problem.newton_solve();
 
  // Document the current solution
  problem.doc_solution();
 
} //End of main
