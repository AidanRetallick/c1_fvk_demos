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
#include "C1_foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"
#include "C1_basis/my_geom_object.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

//                      OUTLINE OF PROBLEM CONSTRUCTION
// The basic constuction is much the same as the usual order of things in a
// problem. Underneath is the order of actions (with stars next to non actions
// that are unique to these types of problems).
// 1.  Setup mesh parameters
// 2.  Build the mesh
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
// 3.* (For 5 order boundary representation) The second tangential derivative
//     of the parametric function defining the curve section.
// 4.  A unit normal and tangent to each curve section and corresponding
//     derivatives, to allow the rotation of boundary coordinates.
// It also convenient to define:
// 1.  An inverse function (x,y) -> s (the arc coordinate) to help in setting
//     up the nodal positions in terms of this parametric coordinate.

namespace Parameters
{
 /// Shape of the domain
 double A = 1.0;
 double B = 1.0;
 /// The plate thickness
 double Thickness = 0.1;
 /// Poisson ratio
 double Nu = 0.5;
 /// Membrane coupling coefficient
 double Eta = 12*(1-Nu*Nu)/(Thickness*Thickness);

 /// Pressure magnitude
 double P_mag = 0.1;
 /// In-plane traction magnitude
 double T_mag = 0.01;


 //                     PARAMETRIC BOUNDARY DEFINITIONS
 /// Here we create the geom objects for the Parametric Boundary Definition
 CurvilineCircleTop parametric_curve_top;
 CurvilineCircleBottom parametric_curve_bottom;

 /// The normal and tangential directions. We need the derivatives so we can form
 /// The Hessian and the Jacobian of the rotation
 void get_normal_and_tangent(const Vector<double>& x,
			     Vector<double>& n,
			     Vector<double>& t,
			     DenseMatrix<double>& Dn,
			     DenseMatrix<double>& Dt)
 {
  double mag = sqrt(x[0]*x[0] + x[1]*x[1]);
  
  // Fill in the normal and derivatives of the normal
  n[0] = x[0]/mag;
  n[1] = x[1]/mag;

  // The (x,y) derivatives of the (x,y) components
  Dn(0,0) = x[1]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
  Dn(1,0) =-x[1]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
  Dn(0,1) =-x[0]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
  Dn(1,1) = x[0]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);

  // Fill in the tangent and derivatives of the tangent
  t[0] =-x[1]/mag;
  t[1] = x[0]/mag;

  Dt(0,0) =-Dn(1,0);
  Dt(1,0) = Dn(0,0);
  Dt(0,1) =-Dn(1,1);
  Dt(1,1) = Dn(0,1);
 }

 //                           PROBLEM DEFINITIONS
 /// Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& x, double& pressure)
 {
  pressure = P_mag;
 }

 /// Pressure wrapper so we can output the pressure function
 void get_pressure(const Vector<double>& X, Vector<double>& pressure)
 {
  pressure.resize(1);
  get_pressure(X,pressure[0]);
 }

 /// Assigns the value of in plane forcing depending on the position (x,y)
 void get_in_plane_force(const Vector<double>& x, Vector<double>& grad)
 {
  grad[0] = T_mag * (1.0 - x[0]*x[0] - x[1]*x[1]);
  grad[1] = T_mag * (1.0 - x[0]*x[0] - x[1]*x[1]);
 }

 // This metric will flag up any non--axisymmetric parts
 void axiasymmetry_metric(const Vector<double>& x,
			  const Vector<double>& u,
			  const Vector<double>& u_exact,
			  double& error,
			  double& norm)
 {
  // We use the theta derivative of the out of plane deflection
  error = pow((-x[1]*u[1] + x[0]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
  norm  = pow(( x[0]*u[1] + x[1]*u[2])/sqrt(x[0]*x[0]+x[1]*x[1]),2);
 }

 // Get the null function for applying homogenous BCs
 void get_null_fct(const Vector<double>& X, double& exact_w)
 {
  exact_w = 0.0;
 }

 // Get the exact solution(s)
 void dummy_exact_w(const Vector<double>& x, Vector<double>& exact_w)
 {
  // Do nothing -> no exact solution in this case
 }

}

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
 UnstructuredFvKProblem(double element_area = 0.09);

 /// Destructor
 ~UnstructuredFvKProblem()
 {
  Trace_file.close();
  delete (Surface_mesh_pt);
  delete (Bulk_mesh_pt);
  // Clean up memory
  delete Outer_boundary_pt;
  delete Outer_boundary_ellipse_pt;
  delete Outer_curvilinear_boundary_pt[0];
  delete Outer_curvilinear_boundary_pt[1];
  delete Inner_open_boundaries_pt[0];
  delete Inner_open_boundaries_pt[1];
  delete Boundary2_pt;
  delete Boundary3_pt;
 };

 /// Setup and build the mesh
 void build_mesh();

 /// Update after solve (empty)
 void actions_after_newton_solve()
 {
  // No actions before newton solve
 }

 /// Pin the in-plane displacements (dofs 0,1) and set to zero at centre
 void pin_in_plane_displacements_at_centre_node();

 /// Pin all displacements and rotation (dofs 0-4) at the centre
 void pin_all_displacements_and_rotation_at_centre_node();
 
 /// Update the problem specs before solve: Re-apply boundary conditions
 /// Empty as the boundary conditions stay fixed
 void actions_before_newton_solve()
 {
  // Reapply boundary conditions
  apply_boundary_conditions();
 }

 /// Doc the solution
 void doc_solution(const std::string& comment="");

 /// Overloaded version of the problem's access function to
 /// the mesh. Recasts the pointer to the base Mesh object to
 /// the actual mesh type.
 TriangleMesh<ELEMENT>* mesh_pt()
 {
  return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt());
 }

 /// Doc info object for labeling output
 DocInfo Doc_info;

private:

 /// Triangle Mesh Parameter Data
 // This is the data used to set-up the mesh, we need to store the pointers
 // HERE otherwise we will not be able to clean up the memory once we have
 // finished the problem.
 Ellipse* Outer_boundary_ellipse_pt;
 /// The outer curves
 Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;
 /// The Internal curves
 Vector<TriangleMeshOpenCurve *> Inner_open_boundaries_pt;
 /// The close outer boundary
 TriangleMeshClosedCurve* Outer_boundary_pt;
 /// The first of the internal boundaries
 TriangleMeshPolyLine* Boundary2_pt;
 /// The second of the internal boundaries
 TriangleMeshPolyLine* Boundary3_pt;

 /// Actions to be performed after read in of meshes
 void actions_after_read_unstructured_meshes()
 {
  // Curved Edges need to be upgraded after the rebuild
  upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
  upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom(Bulk_mesh_pt);
  // Make the problem fully functional
  complete_problem_setup();
  // Apply any boundary conditions
  apply_boundary_conditions();
 }

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup();

 /// Trace file to document norm of solution
 ofstream Trace_file;

 // Keep track of boundary ids
 enum
  {
   Outer_boundary0 = 0,
   Outer_boundary1 = 1,
   Inner_boundary0 = 2,
   Inner_boundary1 = 3
  };

 double Element_area;

 /// Loop over all curved edges, then loop over elements and upgrade
 /// them to be curved elements
 void upgrade_edge_elements_to_curve(const unsigned &b,
				     Mesh* const &bulk_mesh_pt);

 /// Loop over all edge elements and rotate the Hermite degrees of freedom
 /// to be in the directions of the two in-plane vectors specified in Parameters
 void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

 /// Delete traction elements and wipe the surface mesh
 void delete_traction_elements(Mesh* const &surface_mesh_pt);

 /// Pointer to "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to "surface" mesh
 Mesh* Surface_mesh_pt;

}; // end_of_problem_class



//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
 :
 Element_area(element_area)
{
 // Build the mesh
 build_mesh();

 // Curved Edge upgrade
 upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
 upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);

 // Rotate degrees of freedom
 rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

 // Store number of bulk elements
 complete_problem_setup();

 char filename[100];
 sprintf(filename, "RESLT/trace.dat");
 Trace_file.open(filename);

 oomph_info << "Number of equations: "
	    << assign_eqn_numbers() << '\n';
} // end Constructor

/// Set up and build the mesh
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
 Vector<double> zeta(1);
 Vector<double> posn(2);
 
 //Outer boundary
 //--------------
 
 double A = 1.0;
 double B = 1.0;
 Outer_boundary_ellipse_pt = new Ellipse(A, B);
 
 //First bit
 double zeta_start = 0.0;
 double zeta_end = MathematicalConstants::Pi;
 unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(Element_area));
 
 Outer_curvilinear_boundary_pt.resize(2);
 Outer_curvilinear_boundary_pt[0] =
  new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
			    zeta_end, nsegment, Outer_boundary0);
 
 //Second bit
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment = (int)(MathematicalConstants::Pi/sqrt(Element_area));
 Outer_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(Outer_boundary_ellipse_pt, zeta_start,
			    zeta_end, nsegment, Outer_boundary1);
 
 Outer_boundary_pt =
  new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);
 
 // Internal open boundaries
 // Total number of open curves in the domain
 unsigned n_open_curves = 2;
 // We want internal open curves
 Inner_open_boundaries_pt.resize(n_open_curves);
 
 // Internal bit - this means we can have a boundary which is just the centre
 // We start by creating the internal boundaries
 // The boundary 2 is defined by its two vertices
 // Open curve 1
 Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
 vertices[0][0] =-0.5;
 vertices[0][1] = 0.0;

 vertices[1][0] = 0.5;
 vertices[1][1] = 0.0;
 unsigned boundary_id = Inner_boundary0;

 Boundary2_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);
 // Open Curve 2
 vertices[0][0] = 0.0;
 vertices[0][1] =-0.5;

 vertices[1][0] = 0.0;
 vertices[1][1] = 0.5;
 boundary_id = Inner_boundary1;

 Boundary3_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);

 // Each internal open curve is defined by a vector of
 // TriangleMeshCurveSection,
 // on this example we only need one curve section for each internal boundary
 Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
 internal_curve_section1_pt[0] = Boundary2_pt;

 Vector<TriangleMeshCurveSection *> internal_curve_section2_pt(1);
 internal_curve_section2_pt[0] = Boundary3_pt;

 // The open curve that define this boundary is composed of just one
 // curve section
 Inner_open_boundaries_pt[0] =
  new TriangleMeshOpenCurve(internal_curve_section1_pt);

 Inner_open_boundaries_pt[1] =
  new TriangleMeshOpenCurve(internal_curve_section2_pt);

 //Create the mesh
 //---------------
 //Create mesh parameters object
 TriangleMeshParameters mesh_parameters(Outer_boundary_pt);
 
 mesh_parameters.element_area() = Element_area;
 
 // Specify the internal open boundaries
 mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;
 
 // Build an assign bulk mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);
 
 // Create "surface mesh" that will contain only the prescribed-traction
 // elements. The constructor creates the mesh without adding any nodes
 // elements etc.
 Surface_mesh_pt =  new Mesh;
 
 //Add two submeshes to problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 
 // Combine submeshes into a single Mesh
 build_global_mesh();
}// end build_mesh



//==start_of_pin_in_plane_displacements_at_centre_node==========================
/// Set boundary condition exactly, and complete the build of
/// all elements
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_in_plane_displacements_at_centre_node()
{
 // Pin the node that is at the centre in the domain!
 // Get the num of nods on internal_boundary 2
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
   // If the node is on the other internal boundary too
   if( nod_pt->is_on_boundary(3))
    {
     // Pin it! It's the centre of the domain!
     // In-plane dofs are always 0 and 1
     nod_pt->pin(0);
     nod_pt->set_value(0,0.0);
     nod_pt->pin(1);
     nod_pt->set_value(1,0.0);
    }
  }
}




//==start_of_pin_all_displacements_and_rotation_at_centre_node======================
/// Set boundary condition exactly, and complete the build of
/// all elements
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_all_displacements_and_rotation_at_centre_node()
{
 // Pin the node that is at the centre in the domain!
 // Get the num of nods on internal_boundary 2
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
 for (unsigned inod=0;inod<num_int_nod;inod++)
  {
   // Get node point
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
   // If the node is on the other internal boundary too
   if( nod_pt->is_on_boundary(3))
    {
     // Pin it! It's the centre of the domain!
     // In-plane dofs are always 0 and 1
     // Out of plane displacement is 2, x and y derivatives are 3 and 4.
     nod_pt->pin(0);
     nod_pt->set_value(0,0.0);
     nod_pt->pin(1);
     nod_pt->set_value(1,0.0);
     nod_pt->pin(2);
     nod_pt->set_value(2,0.0);
     nod_pt->pin(3);
     nod_pt->set_value(3,0.0);
     nod_pt->pin(4);
     nod_pt->set_value(4,0.0);
    }
  }
}



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
 // The node at the very centre should be pinned when we have "do-nothing"
 // conditions on the in-plane displacements (i.e stress free b/c)
 pin_in_plane_displacements_at_centre_node();

 // Complete the build of all elements so they are fully functional
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the pressure function pointers and the physical constants
   el_pt->pressure_fct_pt() = &Parameters::get_pressure;
   el_pt->in_plane_forcing_fct_pt() = &Parameters::get_in_plane_force;
   // There is no error metric in this case
   el_pt->error_metric_fct_pt() = &Parameters::axiasymmetry_metric;
   el_pt->nu_pt() = &Parameters::Nu;
   el_pt->eta_pt() = &Parameters::Eta;
  }

 // Set the boundary conditions
 apply_boundary_conditions();
}



//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
 // Set the boundary conditions
 // Clamp half the disk
 unsigned nbound = Outer_boundary0 + 1;
 for(unsigned b=0;b<nbound;b++)
  {
   const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<nb_element;e++)
    {
     // Get pointer to bulk element adjacent to b
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));
     // A true clamp, so we set everything except the second normal to zero
     for(unsigned idof=0; idof<6; ++idof)
      {
       // Cannot set second normal derivative
       if(idof!=3)
	{
	 el_pt->fix_out_of_plane_displacement_dof(idof,b,Parameters::get_null_fct);
	}
      }
    }
  }
} // end set bc



//==============================================================================
/// A function that upgrades straight sided elements to be curved. This involves
// Setting up the parametric boundary, F(s) and the first derivative F'(s)
// We also need to set the edge number of the upgraded element and the positions
// of the nodes j and k (defined below) and set which edge (k) is to be exterior
/*            @ k                                                             */
/*           /(                                                               */
/*          /. \                                                              */
/*         /._._)                                                             */
/*      i @     @ j                                                           */
// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
// representation. That is we need to have a continuous 2nd derivative defined
// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
// as well.
//==start_of_upgrade_edge_elements==============================================
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT >::
upgrade_edge_elements_to_curve(const unsigned &ibound, Mesh* const &bulk_mesh_pt)
{
 // These depend on the boundary we are on
 CurvilineGeomObject* parametric_curve_pt;

 // Define the functions for each part of the boundary
 switch (ibound)
  {
  case 0:
   parametric_curve_pt = &Parameters::parametric_curve_top;
   break;
  case 1:
   parametric_curve_pt = &Parameters::parametric_curve_bottom;
   break;
  default:
   throw OomphLibError("Unexpected boundary number. Please add additional \
curved boundaries as required.", OOMPH_CURRENT_FUNCTION,
		       OOMPH_EXCEPTION_LOCATION);
   break;
  } // end parametric curve switch

 // Loop over the bulk elements adjacent to boundary ibound
 const unsigned n_els=bulk_mesh_pt->nboundary_element(ibound);
 for(unsigned e=0; e<n_els; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
						bulk_mesh_pt->boundary_element_pt(ibound,e));

   // Initialise enum for the curved edge
   MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

   // Loop over all (three) nodes of the element and record boundary nodes
   unsigned index_of_interior_node=3,nnode_on_neither_boundary = 0;
   const unsigned nnode = 3;
   // Fill in vertices' positions (this step should be moved inside the curveable
   // Bell element)
   Vector<Vector<double> > xn(nnode,Vector<double>(2,0.0));
   for(unsigned n=0;n<nnode;++n)
    {
     Node* nod_pt = bulk_el_pt->node_pt(n);
     xn[n][0]=nod_pt->x(0);
     xn[n][1]=nod_pt->x(1);

     // Check if it is on the outer boundaries
     if(!(nod_pt->is_on_boundary(Outer_boundary0) ||
          nod_pt->is_on_boundary(Outer_boundary1)))
      {
       index_of_interior_node = n;
       ++nnode_on_neither_boundary;
      }
    }// end record boundary nodes

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
			 "Decreasing parametric coordinate. Parametric coordinate must increase \
as the edge is traversed anti-clockwise.",
			 OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    } // end checks

   // Upgrade it
   bulk_el_pt->upgrade_element_to_curved(edge,s_ubar,s_obar,parametric_curve_pt,3);
  }
}// end_upgrade_elements



//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
/// begin rotate_edge_degrees_of_freedom
//======================================================================
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // Loop over the bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Calculate nodes on the relevant boundaries
   const unsigned nnode=3;
   // Count the number of boundary nodes on external boundaries
   Vector<unsigned> boundary_nodes;
   for (unsigned n=0; n<nnode;++n)
    {
     // If on_external_boundary
     if (el_pt->node_pt(n)->is_on_boundary(0))
      { boundary_nodes.push_back(n); }
     else if (el_pt->node_pt(n)->is_on_boundary(1))
      { boundary_nodes.push_back(n); }
    }

   // If the element has nodes on the boundary, rotate the Hermite dofs
   if(!boundary_nodes.empty())
    {
     // Rotate the nodes by passing the index of the nodes and the
     // normal / tangent vectors to the element
     el_pt->set_up_rotated_dofs(boundary_nodes.size(),
				boundary_nodes,
				&Parameters::get_normal_and_tangent);
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
 unsigned npts = 5;

 sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();
 
 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 //double error,norm,dummy_error,zero_norm;
 double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);

 Bulk_mesh_pt->compute_error(some_file,Parameters::dummy_exact_w,
			     dummy_error,zero_norm);
 some_file.close();

 // Doc L2 error and norm of solution
 oomph_info << "Absolute norm of computed solution: " << sqrt(dummy_error)
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

 Trace_file << Parameters::P_mag
            << " " << u_0[0] << '\n';

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 sprintf(filename,"RESLT/L2-norm%i-%f.dat",
	 Doc_info.number(),
	 Element_area);
 some_file.open(filename);

 some_file<<"### L2 Norm\n";
 some_file<<"##  Format: err^2 norm^2 \n";
 // Print error in prescribed format
 some_file<< dummy_error <<" "<< zero_norm <<"\n";
 some_file.close();

 // Increment the doc_info number
 Doc_info.number()++;

 // Clean up
 delete Mesh_as_geom_obj_pt;
} // end of doc



//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>
::delete_traction_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements



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
 // Directory for solution
 string output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &Parameters::Nu);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &Parameters::P_mag);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--eta", &Parameters::Eta);

 // Element Area (no larger element than 0.09)
 double element_area=0.09;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4> >
  problem(element_area);

 // Set up some problem paramters
 problem.max_residuals()=1e3;
 problem.max_newton_iterations()=20;

 // Do the newton solves
 oomph_info<<"Solving for p=" << Parameters::P_mag << "\n";
 problem.newton_solve();

 // Document
 problem.doc_solution();
 oomph_info << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")"
            << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << std::endl;
 // Print success
 oomph_info<<"Exiting Normally\n";
} //End of main
