// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#include <fenv.h>

// Generic routines
#include "generic.h"

// The equations
#include "c1_foeppl_von_karman.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;


//===========================================================================
/// Namespace for problem parameters
//===========================================================================
namespace Parameters
{
  /// Square edge length
  double L = 1.0;

  /// The plate thickness
  double Thickness = 0.001;

  /// Poisson ratio
  double Nu = 0.5;

  /// FvK parameter
  double Eta = 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);

  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Magnitude of shear stress
  double T_mag = 0.0;

  /// Element size
  double Element_area = 0.02;

  /// Function call to update dependent parameters (Eta)
  void update_dependent_parameters()
  {
    Eta = 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);
  }

  /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    pressure = P_mag;
  }

  /// Shear stress depending on the position (x,y)
  void get_in_plane_force(const Vector<double>& x, Vector<double>& tau)
  {
    tau[0] = T_mag;
    tau[1] = T_mag;
  }


  //-------- Boundary conditions -----------------------------------------------
  /// Function to specify boundary conditions (here all homogeneous) as a
  /// function of both coordinates (This is convenient for this problem; other
  /// interfaces that specify boundary conditions in terms of boundary
  /// coordinate exist).
  void get_null_fct(const Vector<double>& x, double& value)
  {
    value = 0.0;
  }

} // namespace Parameters


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


  /// Actions to complete before each Newton solve
  void actions_before_newton_solve()
  {
    // Always update the dependent parameters before a newton solve
    Parameters::update_dependent_parameters();

    // Print a solve header to output
    oomph_info << "-------------------------------------------------------"
               << std::endl;
    oomph_info << "Solving for P = " << Parameters::P_mag << std::endl;
    oomph_info << "Solving for T = " << Parameters::T_mag << std::endl;
    oomph_info << "         step = " << Doc_info.number() << std::endl;
    oomph_info << "-------------------------------------------------------"
               << std::endl;
  }


  /// Make the problem linear (biharmonic) by pinning all in-plane dofs and
  /// setting eta=0
  void make_linear()
  {
    // Remove stretching coupling
    Parameters::Eta = 0.0;

    // Pin all in-plane displacements
    unsigned n_node = Bulk_mesh_pt->nnode();
    for (unsigned i_node = 0; i_node < n_node; i_node++)
    {
      Bulk_mesh_pt->node_pt(i_node)->pin(0);
      Bulk_mesh_pt->node_pt(i_node)->set_value(0, 0.0);
      Bulk_mesh_pt->node_pt(i_node)->pin(1);
      Bulk_mesh_pt->node_pt(i_node)->set_value(1, 0.0);
    }
  } // End make_linear()

  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

private:
  /// Setup and build the mesh
  void build_mesh();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

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
  Doc_info.number() = 0;

  // Build the mesh
  build_mesh();

  // Complete problem setup
  complete_problem_setup();


  // Output parameters
  oomph_info << "Problem parameters:\n"
             << "L            " << Parameters::L << std::endl
             << "thickness    " << Parameters::Thickness << std::endl
             << "nu           " << Parameters::Nu << std::endl
             << "eta          " << Parameters::Eta << std::endl
             << "Element area " << Parameters::Element_area << std::endl;


  // Open trace file
  char filename[100];
  strcpy(filename, (Doc_info.directory() + "/trace.dat").c_str());
  Trace_file.open(filename);


  // Assign equation numbers
  oomph_info << "Number of equations: " << assign_eqn_numbers() << '\n';

} // end Constructor


//==start_of_build_mesh====================================================
/// Build the rectangular mesh
//=========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
  /*================================
    Rectangular mesh boundary

    V3       E2       V2
    O-----------------O     ^
    |                 |     |
    |        (0,0)    |     |
 E3 |        x        | E1  1
    |                 |     |
    |                 |     |
    O-----------------O     v
    V0       E0       V1
    <--------L-------->
    ================================*/


  double length = Parameters::L;


  // Declare the vertices...
  Vector<double> vertex0(2, 0.0), vertex1(2, 0.0), vertex2(2, 0.0),
    vertex3(2, 0.0);

  // ...and place them at the corners of the rectangle
  vertex0[0] = -0.5 * length;
  vertex0[1] = -0.5;
  vertex1[0] = 0.5 * length;
  vertex1[1] = -0.5;
  vertex2[0] = 0.5 * length;
  vertex2[1] = 0.5;
  vertex3[0] = -0.5 * length;
  vertex3[1] = 0.5;

  // Declare the edges...
  Vector<Vector<double>> edge0(2, Vector<double>(2, 0.0)),
    edge1(2, Vector<double>(2, 0.0)), edge2(2, Vector<double>(2, 0.0)),
    edge3(2, Vector<double>(2, 0.0));

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
  Triangle_mesh_parameters.element_area() = Parameters::Element_area;

  // Build  bulk mesh
  Bulk_mesh_pt = new TriangleMesh<ELEMENT>(Triangle_mesh_parameters);

  // Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);

  // Combine submeshes into a single Mesh (bit over the top here; could
  // have assigned bulk mesh to mesh_pt() directly).
  build_global_mesh();

} // end build_mesh



//==start_of_complete======================================================
/// Set boundary conditions and complete the build of
/// all elements
//=========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the pressure & temperature function pointers and the physical
    // constants
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->in_plane_forcing_fct_pt() = &Parameters::get_in_plane_force;

    // Assign the parameter pointers for the element
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_pt() = &Parameters::Eta;
  }

  // Set the boundary conditions
  apply_boundary_conditions();

  // Find element (and local coordinate within it), that contains the
  // central point
  Vector<double> origin(2, 0.0);

  // Create a MeshAsGeomObject from the Mesh:
  MeshAsGeomObject mesh_as_geom_object(Bulk_mesh_pt);

  // Make space for local coordinate in element pointed to by Central_element_pt
  // that contains central point
  Central_point_local_coord.resize(2);

  // Find it!
  mesh_as_geom_object.locate_zeta(
    origin, Central_element_geom_obj_pt, Central_point_local_coord);

  // Did that work out?
  if (Central_element_geom_obj_pt == 0)
  {
    throw OomphLibError("couldn't find central element",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

} // end of complete


//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
  //------------------------------------------------------------------
  // Boundary conditions for FvK elements are complicated and we provide an
  // illustration of how to apply all physically relevant boundary conditions
  // for problems with axis-aligned boundaries here. Other tutorials/driver
  // codes explain what to do about (i) straight boundaries that are not aligned
  // with the coordinate axes and (ii) curvilinear boundaries.
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
  //     enumerated from 0 to 5 in the order w, w_x, w_y, w_xx, w_xy,
  //     w_yy. These values are only stored at three vertices of the element.
  //
  // Given that the book-keeping for which node stores which type of degree of
  // freedom is complicated, we let the element do the relevant assigments for
  // us. For this purpose the FvK elements provide two member functions:
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
  //          function of the coordinate, x, a 2D vector. Note that,
  //          since the boundary is assumed to be aligned with the
  //          coordinate axes, one of the two coordinates will be
  //          irrelevant, but it is still passed to the function.
  //
  // So, if the function fix_in_plane_displacement_dof(idof, b, fct_pt) is
  // called with idof=1 and b=3, say, the y-in-plane displacement is pinned for
  // all the element's nodes (if any) that are located on mesh boundary b. The
  // value of the y-in-plane displacement is set to whatever the function
  // pointed to by fct_pt computes when evaluated at the nodal coordinate.
  //
  // Similarly,
  //
  //      fix_out_of_plane_displacement_dof(idof, b, fct_pt);
  //
  // assigns boundary conditions for the out-of-plane displacements via
  // the specification of
  //
  // idof   : the enumeration of the dof in the scheme listed above,
  //          so idof can take values 0, 1. 2, 3, 4 or 5
  // b      : the mesh boundary along which the boundary condition is
  //          to be applied
  // fct_pt : a function pointer to a global function with arguments
  //          (const Vector<double> x, double& value) which computes
  //          the value for the relevant out-of-plane displacement (or 
  //          its derivative; depending on what the dof represents) as a
  //          function of the coordinate, x, a 2D vector. 
  //
  // So, if the function fix_out_of_plane_displacement_dof(idof, b, fct_pt) is called 
  //  with idof=2 and b=3, say, the y-derivative of the out-of-plane displacement 
  // is pinned for all the element's nodes (if any) that are located on mesh boundary 3. 
  // The value of this derivative is set to whatever the function pointed to by fct_pt 
  // computes when evaluated at the nodal coordinate. 
  //
  // Using the conventions introduced above, the following vectors identify the
  // in-plane and out-of-plane degrees of freedom to be pinned for various
  // physically meaningful boundary conditions:

  // In-plane dofs:
  //---------------
  // |  0  |  1  |
  // |  ux |  uy |

  // Case: Pin x in-plane displacement only:
  const Vector<unsigned> pin_ux_pinned_dof{0};

  // Case: Pin y in-plane displacement only:
  const Vector<unsigned> pin_uy_pinned_dof{1};

  // Case: Pin both in-plane displacements
  const Vector<unsigned> pin_ux_and_uy_pinned_dof{0, 1};


  // Possible boundary conditions for out-of-plane displacements: Given that the
  // out-of-plane displacements feature in the fourth-order biharmonic operator,
  // we can apply boundary conditions on w and dw/dn, where n is the coordinate
  // direction normal to the (assumed to be axis aligned!) boundary. However if
  // w is given along the entire boundary (parametrised by the tangential
  // coordinate, t) we also know what dw/dt and d^2w/dt^2 are. Likewise if dw/dn
  // is known along the whole boundary we also know d^2w/dndt. In the various
  // cases below we identify physical scenarios of a pinned edge (w given, dw/dn
  // left free); a vertically sliding edge (w left free; dw/dn given) and fully
  // clamped (w and dw/dn given). Together with the two possible orientations of
  // the axis aligned boundaries (x aligned or y aligned) we get six different
  // cases:

  // Out-of-plane dofs:
  //-------------------
  // |  0  |  1  |  2  |  3  |  4  |  5  |
  // |  w  | w_x | w_y | w_xx| w_xy| w_yy|

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary where
  // the outer unit normal points in the postive or negative x direction, so x
  // is constant and y varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dy and d^2w/dy^2
  const Vector<unsigned> pinned_edge_xn_dof{0, 2, 5};

  // Case: The plate is pinned (w given, dw/dn left free) along a boundary where
  // the outer unit normal points in the postive or negative y direction, so y
  // is constant and x varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx and d^2w/dx^2
  const Vector<unsigned> pinned_edge_yn_dof{0, 1, 3};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary
  // where the outer unit normal points in the postive or negative x direction,
  // so x is constant and y varies along the boundary.  We therefore have to pin
  // (and assign values for) dw/dx and d^2w/dxdy
  const Vector<unsigned> sliding_clamp_xn_dof{1, 4};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary
  // where the outer unit normal points in the postive or negative y direction,
  // so y is constant and x varies along the boundary.  We therefore have to pin
  // (and assign values for) dw/dy and d^2w/dxdy
  const Vector<unsigned> sliding_clamp_yn_dof{2, 4};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary where
  // the outer unit normal points in the postive or negative x direction, so x
  // is constant and y varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx, dw/dy, d^2w/dxdy and d^2w/dy^2
  const Vector<unsigned> fully_clamped_xn_dof{0, 1, 2, 4, 5};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary where
  // the outer unit normal points in the postive or negative y direction, so y
  // is constant and x varies along the boundary.  We therefore have to pin (and
  // assign values for) w, dw/dx, dw/dy, d^2w/dx^2 and d^2w/dxdy
  const Vector<unsigned> fully_clamped_yn_dof{0, 1, 2, 3, 4};


  //------------------------------------------------------------------
  //------------------------------------------------------------------


  // Vector containers to store which boundary conditions we are applying to
  // each edge. (outlined above)
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
  pinned_w_dofs[0] = pinned_edge_yn_dof;
  pinned_w_dofs[1] = pinned_edge_xn_dof;
  pinned_w_dofs[2] = pinned_edge_yn_dof;
  pinned_w_dofs[3] = pinned_edge_xn_dof;


  // Loop over all the boundaries in our bulk mesh
  unsigned n_bound = Bulk_mesh_pt->nboundary();
  for (unsigned b = 0; b < n_bound; b++)
  {
    // Number of elements on b
    const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);

    // Number of dofs we are pinning on boundary b
    const unsigned n_pinned_u_dofs = pinned_u_dofs[b].size();
    const unsigned n_pinned_w_dofs = pinned_w_dofs[b].size();

    // Loop over the elements on boundary b
    for (unsigned e = 0; e < nb_element; e++)
    {
      // Get pointer to bulk element adjacent to b
      ELEMENT* el_pt =
        dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

      // Pin in-plane dofs (enumerated as explained above) for all nodes on
      // boundary b. Here we're applying homogeneous BCs so all pinned values
      // are simply set to zero using Parameters::get_null_fct.
      for (unsigned i = 0; i < n_pinned_u_dofs; i++)
      {
        unsigned idof_to_be_pinned = pinned_u_dofs[b][i];
        el_pt->fix_in_plane_displacement_dof(
          idof_to_be_pinned, b, Parameters::get_null_fct);
      }
      // Pin out-of-plane dofs (enumerated as explained above) for all nodes on
      // boundary b. Here we're applying homogeneous BCs so all pinned values
      // are simply set to zero using Parameters::get_null_fct.
      for (unsigned i = 0; i < n_pinned_w_dofs; i++)
      {
        unsigned idof_to_be_pinned = pinned_w_dofs[b][i];
        el_pt->fix_out_of_plane_displacement_dof(
          idof_to_be_pinned, b, Parameters::get_null_fct);
      }
    } // end for loop over elements on b
  } // end for loop over boundaries

} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  ofstream some_file;
  char filename[100];


  // Number of plot points for coarse output (just showing the
  // element outline)
  unsigned npts = 2;
  sprintf(filename,
          "%s/coarse_soln_%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();

  // Number of plot points for fine (full) output. Lots of plot points so
  // we see the goodness of the high-order Hermite/Bell
  // interpolation.
  npts = 10;
  sprintf(filename,
          "%s/soln_%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();


  // Write the load magnitude and centrepoint displacements to trace file
  Vector<double> w_centre(12, 0.0);

  w_centre = dynamic_cast<ELEMENT*>(Central_element_geom_obj_pt)
               ->interpolated_fvk_disp_and_deriv(Central_point_local_coord);

  Trace_file << Parameters::P_mag << " " << Parameters::T_mag << " "
             << w_centre[0] << " " << w_centre[1] << " " << w_centre[2] << " "
             << Doc_info.number() << endl;

  // Bump
  Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
/// Driver code for demo of unstructured C1 Foeppl-von Karman
/// elements
//============================================================
int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

  // Create the problem, using FvK elements derived from TElement<2,4>
  // elements (with 4 nodes per element edge and 10 nodes overall).
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem;

  problem.max_residuals() = 1.0e3;
  problem.max_newton_iterations() = 30;

  // Set pressure
  Parameters::P_mag = 1.0e-5;
  // Set the Poisson ratio
  Parameters::Nu = 0.5;

  // Document the initial state
  problem.doc_solution();
  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

} // End of main
