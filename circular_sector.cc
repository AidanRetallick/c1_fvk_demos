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

// Floating point environment
#include <fenv.h>

// Generic oomph-lib routines
#include "generic.h"
// The equations
#include "c1_foeppl_von_karman.h"
// The mesh
#include "meshes/triangle_mesh.h"

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
  /// Opening angle of the domain corner
  double Alpha = Pi / 4.0;

  /// The plate thickness
  double Thickness = 0.01;

  /// Poisson ratio
  double Nu = 0.5;

  /// Membrane coupling coefficient
  double Eta = 12.0 * (1.0 - Nu * Nu) / (Thickness * Thickness);

  /// Membrane coupling coefficient for linear bending
  double Eta_linear = 0.0;

  /// Boundary wave amplitude
  double Boundary_amp = 0.1;

  /// Magnitude of pressure
  double P_mag = 10.0;

  //                     PARAMETRIC BOUNDARY DEFINITIONS
  // Here we create the geom objects for the Parametric Boundary Definition
  // (needed to update boundary elements to curved)

  /// Curved arc of sector
  CurvilineCircleTop Circular_arc;

  /// The lower (x-axis) radius
  CurvilineLine Lower_line(0.0, false);

  /// The upper radial line
  CurvilineLine Upper_line(Alpha, true);

  Vector<CurvilineGeomObject*> Parametric_curve_pt = {
    &Lower_line, &Circular_arc, &Upper_line};

  //                           PROBLEM DEFINITIONS
  /// Assigns the value of pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    // No pressure
    pressure = P_mag;
  }

  /// Pressure wrapper so we can output the pressure function
  void get_pressure(const Vector<double>& X, Vector<double>& pressure)
  {
    pressure.resize(1);
    get_pressure(X, pressure[0]);
  }

  /// Assigns the value of in plane forcing depending on the position (x,y)
  void get_in_plane_force(const Vector<double>& x, Vector<double>& grad)
  {
    // No in plane force
    grad[0] = 0.0;
    grad[1] = 0.0;
  }

  // This metric will flag up any non--axisymmetric parts
  void axiasymmetry_metric(const Vector<double>& x,
                           const Vector<double>& u,
                           const Vector<double>& u_exact,
                           double& error,
                           double& norm)
  {
    // We use the theta derivative of the out of plane deflection
    error =
      pow((-x[1] * u[1] + x[0] * u[2]) / sqrt(x[0] * x[0] + x[1] * x[1]), 2);
    norm =
      pow((x[0] * u[1] + x[1] * u[2]) / sqrt(x[0] * x[0] + x[1] * x[1]), 2);
  }

  /// Get the exact solution(s)
  void dummy_exact_w(const Vector<double>& x, Vector<double>& exact_w)
  {
    // Do nothing -> no exact solution in this case
  }

  //---------------------------------------------------------------------
  // Functions to assign boundary conditions
  // (e.g. sin along arc 0<theta<Alpha)

  /// Sin along circular arc for w BC
  void get_w_along_arc(const Vector<double>& x, double& value)
  {
    value = Boundary_amp * sin(2.0 * Pi * atan2(x[1], x[0]) / Alpha);
  }

  /// Cos along circular arc for dw/dt BC
  void get_dwdt_along_arc(const Vector<double>& x, double& value)
  {
    value = 2.0 * Pi / Alpha * Boundary_amp *
            cos(2.0 * Pi * atan2(x[1], x[0]) / Alpha);
  }

  /// -Sin along circular arc for d2w/dt2 BC
  void get_d2wdt2_along_arc(const Vector<double>& x, double& value)
  {
    value = -pow(2.0 * Pi / Alpha, 2.0) * Boundary_amp *
            sin(2.0 * Pi * atan2(x[1], x[0]) / Alpha);
  }

  /// Null function for any zero (homogenous) BCs
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
  UnstructuredFvKProblem(double element_area = 0.09);

  /// Destructor
  ~UnstructuredFvKProblem()
  {
    Trace_file.close();
    delete (Constraint_mesh_pt);
    delete (Bulk_mesh_pt);
    // Clean up memory
    delete Boundary_pt;
    delete Outer_boundary_ellipse_pt;
    delete Outer_boundary_curvilines_pt[0];
    delete Outer_boundary_curvilines_pt[1];
    delete Outer_boundary_curvilines_pt[2];
  };

  /// Setup and build the mesh
  void build_mesh();

  /// Duplicate corner nodes and create constraint elements at those corners
  void duplicate_corner_nodes();

  /// Update after solve (empty)
  void actions_after_newton_solve()
  {
    // No actions before newton solve
  }

  /// Pin the in-plane displacements and set to zero at centre
  void pin_in_plane_displacements_at_centre_node();

  /// Update the problem specs before solve: Re-apply boundary conditions
  /// Empty as the boundary conditions stay fixed
  void actions_before_newton_solve()
  {
    // Reapply boundary conditions
    apply_boundary_conditions();
  }

  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

  /// Doc info object for labeling output
  DocInfo Doc_info;

private:
  // This is the data used to set-up the mesh, we need to store the pointers
  // HERE otherwise we will not be able to clean up the memory once we have
  // finished the problem.
  /// Parametrised boundary geometric object
  Ellipse* Outer_boundary_ellipse_pt;
  /// The outer boundary component curves
  Vector<TriangleMeshCurveSection*> Outer_boundary_curvilines_pt;
  /// The closed outer boundary
  TriangleMeshClosedCurve* Boundary_pt;

  /// Actions to be performed after read in of meshes
  void actions_after_read_unstructured_meshes()
  {
    // Curved Edges need to be upgraded after the rebuild
    upgrade_edge_elements_to_curve(Circular_arc_bnum, Bulk_mesh_pt);
    // Rotate degrees of freedom
    rotate_edge_degrees_of_freedom();
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

  /// Pin all in-plane displacements in the domain
  /// (for solving the linear problem)
  void pin_all_in_plane_displacements();

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Keep track of boundary ids
  enum
  {
    Straight_edge_0_bnum = 0,
    Circular_arc_bnum = 1,
    Straight_edge_1_bnum = 2
  };

  /// Maximum element area size
  double Element_area;

  /// Flag to turn the elements into linear beinding
  bool Solve_linear_bending;

  /// Loop over all curved edges, then loop over elements and upgrade
  /// them to be curved elements
  void upgrade_edge_elements_to_curve(const unsigned& b,
                                      Mesh* const& bulk_mesh_pt);

  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in
  /// Parameters
  void rotate_edge_degrees_of_freedom();

  /// Delete traction elements and wipe the surface mesh
  void delete_traction_elements(Mesh* const& surface_mesh_pt);

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to "surface" mesh
  Mesh* Constraint_mesh_pt;

}; // end_of_problem_class


//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
  : Element_area(element_area), Solve_linear_bending(true)
{
  // Build the mesh
  build_mesh();

  // Curved Edge upgrade
  upgrade_edge_elements_to_curve(Circular_arc_bnum, Bulk_mesh_pt);

  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom();

  // Store number of bulk elements
  complete_problem_setup();

  char filename[100];
  sprintf(filename, "RESLT/trace.dat");
  Trace_file.open(filename);

  oomph_info << "Number of equations: " << assign_eqn_numbers() << '\n';

  // Doc the equatiuons
  describe_dofs();

} // end Constructor


/// Set up and build the mesh
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
  Vector<double> zeta(1);
  Vector<double> posn(2);

  // Opening angle
  double alpha = Parameters::Alpha;

  // Vertices
  Vector<Vector<double>> endpoints(2);
  Vector<double> corner(2, 0.0);
  Vector<double> arc_start(2, 0.0);
  Vector<double> arc_end(2, 0.0);
  arc_start[0] = 1.0;
  arc_end[0] = cos(alpha);
  arc_end[1] = sin(alpha);

  // Sector of a circle with radius 1
  double A = 1.0;
  double B = 1.0;
  Outer_boundary_ellipse_pt = new Ellipse(A, B);

  // Three boundary lines
  Outer_boundary_curvilines_pt.resize(3);

  // Straight boundary 0 (lower)
  endpoints[0] = corner;
  endpoints[1] = arc_start;
  Outer_boundary_curvilines_pt[Straight_edge_0_bnum] =
    new TriangleMeshPolyLine(endpoints, Straight_edge_0_bnum);

  // Straight boundary 1 (upper)
  endpoints[0] = arc_end;
  endpoints[1] = corner;
  Outer_boundary_curvilines_pt[Straight_edge_1_bnum] =
    new TriangleMeshPolyLine(endpoints, Straight_edge_1_bnum);

  // Curved arc boundary
  double zeta_start = 0.0;
  double zeta_end = alpha;
  unsigned nsegment = 1.0; //(int)(Pi/sqrt(Element_area));
  Outer_boundary_curvilines_pt[Circular_arc_bnum] =
    new TriangleMeshCurviLine(Outer_boundary_ellipse_pt,
                              zeta_start,
                              zeta_end,
                              nsegment,
                              Circular_arc_bnum);

  // Form closed curve from components
  Boundary_pt = new TriangleMeshClosedCurve(Outer_boundary_curvilines_pt);


  // Create the mesh
  //---------------
  // Create mesh parameters object
  TriangleMeshParameters mesh_parameters(Boundary_pt);

  mesh_parameters.element_area() = Element_area;

  // Build an assign bulk mesh
  Bulk_mesh_pt = new TriangleMesh<ELEMENT>(mesh_parameters);
  Bulk_mesh_pt->setup_boundary_element_info();

  // Build mesh to contain constraint elements
  Constraint_mesh_pt = new Mesh();

  // Split elements that have two boundary edges
  TimeStepper* time_stepper_pt = Bulk_mesh_pt->Time_stepper_pt;
  Bulk_mesh_pt->template split_elements_with_multiple_boundary_edges<ELEMENT>(
    time_stepper_pt);

  // Add extra nodes at boundaries and constrain the dofs there.
  duplicate_corner_nodes();

  // Add two submeshes to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

  // Combine submeshes into a single Mesh
  build_global_mesh();

} // end build_mesh


//==============================================================================
/// Duplicate nodes at corners in order to properly apply boundary
/// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
/// problem in order to constrain continuous interpolation here across its (8)
/// vertex dofs. (Note "corner" here refers to the meeting point of any two
/// sub-boundaries in the closed external boundary)
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::duplicate_corner_nodes()
{
  // Loop over the sections of the external boundary
  unsigned n_bound = Bulk_mesh_pt->nboundary();
  for (unsigned i_bound = 0; i_bound < n_bound; i_bound++)
  {
    // Store the index of the next boundary
    unsigned ip1_bound = (i_bound + 1) % n_bound;
    // Storage for node and el pts at the boundary vertex
    Node* old_node_pt = 0;
    Node* new_node_pt = 0;
    // FiniteElement* left_element_pt = 0;
    FiniteElement* right_element_pt = 0;

    // To find the node between boundaries i and i+1, we loop over all nodes on
    // boundary i until we find the one that is also on i+1, we then loop over
    // all boundary elements on i and i+1 until we find the elements that sit
    // either side of the corner. (there might be a smarter way but this was the
    // first idea I had -- Aidan)

    //----------------------------------------------------------------------
    // First, find corner the node
    unsigned n_b_node = Bulk_mesh_pt->nboundary_node(i_bound);
    for (unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
    {
      // Store the node we are checking
      Node* node_pt = Bulk_mesh_pt->boundary_node_pt(i_bound, i_b_node);

      // If it is on the next boundary we have found the corner node
      if (node_pt->is_on_boundary(ip1_bound))
      {
        // [zdec] debug
        oomph_info << "Found a corner node at " << std::endl
                   << "  (" << node_pt->position(0) << ","
                   << node_pt->position(1) << ")" << std::endl;
        old_node_pt = node_pt;
        break;
      }
    }

    //----------------------------------------------------------------------
    // Find the right (i+1th boundary) side element
    unsigned n_b_el = Bulk_mesh_pt->nboundary_element(ip1_bound);
    for (unsigned i_b_el = 0; i_b_el < n_b_el; i_b_el++)
    {
      // Get the element pointer
      FiniteElement* el_pt =
        Bulk_mesh_pt->boundary_element_pt(ip1_bound, i_b_el);
      // If the corner node pt is in the element we have found the right
      // element
      if (el_pt->get_node_number(old_node_pt) != -1)
      {
        right_element_pt = el_pt;
        break;
      }
    }

    //----------------------------------------------------------------------
    // Now we need to create a new node and substitute the right elements
    // old corner node for this new one
    new_node_pt = right_element_pt->construct_boundary_node(
      right_element_pt->get_node_number(old_node_pt));
    // Copy the position and other info from the old node into the new node
    // [debug]
    oomph_info << "About to copy node data" << std::endl;
    new_node_pt->x(0) = old_node_pt->x(0);
    new_node_pt->x(1) = old_node_pt->x(1);
    oomph_info << "Copied node data" << std::endl;
    // Then we add this node to the mesh
    Bulk_mesh_pt->add_node_pt(new_node_pt);
    // Then replace the old node for the new one on the right boundary
    Bulk_mesh_pt->remove_boundary_node(ip1_bound, old_node_pt);
    Bulk_mesh_pt->add_boundary_node(ip1_bound, new_node_pt);

    //----------------------------------------------------------------------
    // The final job is to constrain this duplication using the specialised
    // Lagrange multiplier elements which enforce equality of displacement and
    // its derivatives either side of this corner.
    CurvilineGeomObject* left_parametrisation_pt =
      Parameters::Parametric_curve_pt[i_bound];
    CurvilineGeomObject* right_parametrisation_pt =
      Parameters::Parametric_curve_pt[ip1_bound];

    // Get the coordinates on each node on their respective boundaries
    Vector<double> left_boundary_coordinate = {
      left_parametrisation_pt->get_zeta(old_node_pt->position())};
    Vector<double> right_boundary_coordinate = {
      right_parametrisation_pt->get_zeta(new_node_pt->position())};

    // Create the constraining element
    DuplicateNodeConstraintElement* constraint_element_pt =
      new DuplicateNodeConstraintElement(old_node_pt,
                                         new_node_pt,
                                         left_parametrisation_pt,
                                         right_parametrisation_pt,
                                         left_boundary_coordinate,
                                         right_boundary_coordinate);

    // Add the constraining element to the mesh
    Constraint_mesh_pt->add_element_pt(constraint_element_pt);
  }
}


//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
  // Complete the build of all elements so they are fully functional
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // Set the pressure function pointers and the physical constants
    el_pt->pressure_fct_pt() = &Parameters::get_pressure;
    el_pt->in_plane_forcing_fct_pt() = &Parameters::get_in_plane_force;
    // There is no error metric in this case
    el_pt->error_metric_fct_pt() = &Parameters::axiasymmetry_metric;
    el_pt->nu_pt() = &Parameters::Nu;
    if (Solve_linear_bending)
    {
      el_pt->eta_pt() = &Parameters::Eta_linear;
    }
    else
    {
      el_pt->eta_pt() = &Parameters::Eta;
    }
  }

  // Set the boundary conditions
  apply_boundary_conditions();

  // Pin in-plane displacements throughout the bulk
  if (Solve_linear_bending)
  {
    pin_all_in_plane_displacements();
  }

  // Update the corner constraintes based on boundary conditions
  unsigned n_el = Constraint_mesh_pt->nelement();
  for (unsigned i_el = 0; i_el < n_el; i_el++)
  {
    dynamic_cast<DuplicateNodeConstraintElement*>(
      Constraint_mesh_pt->element_pt(i_el))
      ->validate_and_pin_redundant_constraints();
  }
}


//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
  //----------------------------------------------------------------------------
  // Sets of possible boundary conditions

  // The free boundary condition is completely unpinned
  static const Vector<unsigned> free{};

  // In-plane dofs:
  // |  0  |  1  |
  // |  un |  ut |
  // Possible boundary conditions on the in-plane displacement
  static const Vector<unsigned> pin_un_dofs{0};
  static const Vector<unsigned> pin_ut_dofs{1};
  static const Vector<unsigned> pin_inplane_dofs{0, 1};

  // Out-of-plane dofs:
  // |  0  |  1  |  2  |  3  |  4  |  5  |
  // |  w  | w_n | w_t | w_nn| w_nt| w_tt|
  // Possible boundary conditions for the out-of-plane displacement
  static const Vector<unsigned> resting_pin_dofs{0, 2, 5};
  static const Vector<unsigned> sliding_clamp_dofs{1, 4};
  static const Vector<unsigned> true_clamp_dofs{0, 1, 2, 4, 5};

  //----------------------------------------------------------------------------
  // Assign boundary conditions to each edge (do circular arc manually)
  Vector<unsigned> straight_edge_0_pinned_u_dofs = pin_inplane_dofs;
  Vector<unsigned> straight_edge_1_pinned_u_dofs = pin_inplane_dofs;
  Vector<unsigned> circular_arc_pinned_u_dofs = pin_inplane_dofs;
  Vector<unsigned> straight_edge_0_pinned_w_dofs = resting_pin_dofs;
  Vector<unsigned> straight_edge_1_pinned_w_dofs = free; // resting_pin_dofs;
  // Vector<unsigned> circular_arc_pinned_w_dofs = true_clamp_dofs;

  // Allocate storage for the number of elements and dofs
  unsigned n_b_element = 0;
  unsigned n_pinned_u_dofs = 0;
  unsigned n_pinned_w_dofs = 0;

  // Loop over the circular arc elements
  n_b_element = Bulk_mesh_pt->nboundary_element(Circular_arc_bnum);
  n_pinned_u_dofs = circular_arc_pinned_u_dofs.size();
  // n_pinned_w_dofs = circular_arc_pinned_w_dofs.size();
  for (unsigned e = 0; e < n_b_element; e++)
  {
    // Get pointer to bulk element adjacent to curved arc
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(Circular_arc_bnum, e));

    // Pin in-plane dofs
    for (unsigned i = 0; i < n_pinned_u_dofs; i++)
    {
      unsigned idof = circular_arc_pinned_u_dofs[i];
      el_pt->fix_in_plane_displacement_dof(
        idof, Circular_arc_bnum, Parameters::get_null_fct);
    }
    // Pin out-of-plane dofs (resting pin -- only set tangent derivatives)
    for (unsigned idof = 0; idof < 6; ++idof)
    {
      switch (idof)
      {
          // [hierher] Make function of arclength rather than global x
        case 0:
          // w
          el_pt->fix_out_of_plane_displacement_dof(
            idof, Circular_arc_bnum, Parameters::get_w_along_arc);
          // Parameters::get_null_fct);
          break;
        case 2:
          // dwdt
          el_pt->fix_out_of_plane_displacement_dof(
            idof, Circular_arc_bnum, Parameters::get_dwdt_along_arc);
          // Parameters::get_null_fct);
          break;
        case 5:
          // d2wdt2
          el_pt->fix_out_of_plane_displacement_dof(
            idof, Circular_arc_bnum, Parameters::get_d2wdt2_along_arc);
          // Parameters::get_null_fct);
          break;
        default:
          // Leave free
          break;
      } // End of switch-case [idof]
    } // End loop over dofs [idof]
  } // End loop over boundary elements [e]

  // Loop over straight side 0 elements and apply homogenous BCs
  n_b_element = Bulk_mesh_pt->nboundary_element(Straight_edge_0_bnum);
  n_pinned_u_dofs = straight_edge_0_pinned_u_dofs.size();
  n_pinned_w_dofs = straight_edge_0_pinned_w_dofs.size();
  for (unsigned e = 0; e < n_b_element; e++)
  {
    // Get pointer to bulk element adjacent to b
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(Straight_edge_0_bnum, e));

    // Pin in-plane dofs
    for (unsigned i = 0; i < n_pinned_u_dofs; i++)
    {
      unsigned idof = straight_edge_0_pinned_u_dofs[i];
      el_pt->fix_in_plane_displacement_dof(
        idof, Straight_edge_0_bnum, Parameters::get_null_fct);
    } // End loop over in-plane dofs [i]
      // Pin out-of-plane dofs
    for (unsigned i = 0; i < n_pinned_w_dofs; i++)
    {
      unsigned idof = straight_edge_0_pinned_w_dofs[i];
      el_pt->fix_out_of_plane_displacement_dof(
        idof, Straight_edge_0_bnum, Parameters::get_null_fct);
    } // End loop over out-of-plane dofs [i]
  } // End loop over boundary elements [e]

  // Loop over straight side 1 elements and apply homogenous BCs
  n_b_element = Bulk_mesh_pt->nboundary_element(Straight_edge_1_bnum);
  n_pinned_u_dofs = straight_edge_1_pinned_u_dofs.size();
  n_pinned_w_dofs = straight_edge_1_pinned_w_dofs.size();
  for (unsigned e = 0; e < n_b_element; e++)
  {
    // Get pointer to bulk element adjacent to b
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(Straight_edge_1_bnum, e));

    // Pin in-plane dofs
    for (unsigned i = 0; i < n_pinned_u_dofs; i++)
    {
      unsigned idof = straight_edge_1_pinned_u_dofs[i];
      el_pt->fix_in_plane_displacement_dof(
        idof, Straight_edge_1_bnum, Parameters::get_null_fct);
    } // End loop over in-plane dofs [i]
      // Pin out-of-plane dofs
    for (unsigned i = 0; i < n_pinned_w_dofs; i++)
    {
      unsigned idof = straight_edge_1_pinned_w_dofs[i];
      el_pt->fix_out_of_plane_displacement_dof(
        idof, Straight_edge_1_bnum, Parameters::get_null_fct);
    } // End loop over out-of-plane dofs [i]
  } // End loop over boundary elements [e]


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
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::upgrade_edge_elements_to_curve(
  const unsigned& ibound, Mesh* const& bulk_mesh_pt)
{
  // These depend on the boundary we are on
  CurvilineGeomObject* parametric_curve_pt = 0;

  // Define the functions for each part of the boundary
  switch (ibound)
  {
    case Circular_arc_bnum:
      parametric_curve_pt = &Parameters::Circular_arc;
      break;
    default:
      throw OomphLibError("Unexpected boundary number. Please add additional \
curved boundaries as required.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      break;
  } // end parametric curve switch

  // Loop over the bulk elements adjacent to boundary ibound
  const unsigned n_els = bulk_mesh_pt->nboundary_element(ibound);
  for (unsigned e = 0; e < n_els; e++)
  {
    // Get pointer to bulk element adjacent to b
    ELEMENT* bulk_el_pt =
      dynamic_cast<ELEMENT*>(bulk_mesh_pt->boundary_element_pt(ibound, e));

    // Initialise enum for the curved edge
    MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

    // Loop over all (three) nodes of the element and record boundary nodes
    unsigned index_of_interior_node = 3;
    unsigned nnode_not_on_curved_boundary = 0;
    const unsigned nnode = 3;
    // Fill in vertices' positions (this step should be moved inside the
    // curveable Bell element)
    Vector<Vector<double>> xn(nnode, Vector<double>(2, 0.0));
    for (unsigned n = 0; n < nnode; ++n)
    {
      Node* nod_pt = bulk_el_pt->node_pt(n);
      xn[n][0] = nod_pt->x(0);
      xn[n][1] = nod_pt->x(1);

      // Check if it is on the curved boundary
      if (!(nod_pt->is_on_boundary(Circular_arc_bnum)))
      {
        index_of_interior_node = n;
        nnode_not_on_curved_boundary++;
      }
    } // end record boundary nodes

    // s at the next (cyclic) node after interior
    const double s_ubar =
      parametric_curve_pt->get_zeta(xn[(index_of_interior_node + 1) % 3]);
    // s at the previous (cyclic) node before interior
    const double s_obar =
      parametric_curve_pt->get_zeta(xn[(index_of_interior_node + 2) % 3]);
    // Assign edge case
    edge = static_cast<MyC1CurvedElements::Edge>(index_of_interior_node);

    // Check nnode_not_on_curved_boundary
    if (nnode_not_on_curved_boundary == 0)
    {
      throw OomphLibError(
        "No interior nodes. One node per CurvedElement must be interior.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    else if (nnode_not_on_curved_boundary > 1)
    {
      throw OomphLibError("Multiple interior nodes. Only one node per "
                          "CurvedElement can be interior.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Check for inverted elements
    if (s_ubar > s_obar)
    {
      throw OomphLibError(
        "Decreasing parametric coordinate. Parametric coordinate must increase \
as the edge is traversed anti-clockwise.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    } // end checks

    // Upgrade it
    bulk_el_pt->upgrade_element_to_curved(
      edge, s_ubar, s_obar, parametric_curve_pt, 5);
  }
} // end_upgrade_elements


//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to
/// set up physical boundary conditions on a curved boundary with Hermite type
/// dofs. For example if we know w(n,t) = f(t) (where n and t are the normal and
/// tangent to a boundary) we ALSO know dw/dt and d2w/dt2. NB no rotation is
/// needed if the edges are completely free!
//======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::rotate_edge_degrees_of_freedom()
{
  // Get the number of boundaries
  unsigned n_bound = 3;

  // Loop over the bulk elements
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Get pointer to bulk element adjacent
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    // [zdec] debug
    oomph_info << "In element " << e << " (" << el_pt << ")" << std::endl;

    // Loop over each boundary and add the boundary parametrisation to the
    // relevant nodes' boundary data
    for (unsigned b = 0; b < n_bound; b++)
    {
      // Calculate nodes on the relevant boundaries
      const unsigned nnode = 3;
      // Count the number of boundary nodes on external boundaries
      Vector<unsigned> boundary_node;
      // Store the boundary coordinates of nodes on the boundaries
      Vector<double> boundary_coordinate_of_node;
      for (unsigned n = 0; n < nnode; ++n)
      {
        // If on external boundary b
        if (el_pt->node_pt(n)->is_on_boundary(b))
        {
          boundary_node.push_back(n);
          double coord = Parameters::Parametric_curve_pt[b]->get_zeta(
            el_pt->node_pt(n)->position());
          boundary_coordinate_of_node.push_back(coord);
        }
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if (!boundary_node.empty())
      {
        // [zdec] debug
        oomph_info << " Nodes ";
        for (unsigned n : boundary_node)
        {
          oomph_info << n << " ";
        }
        oomph_info << " are on boundary " << b << std::endl;
        // Rotate the nodes by passing the index of the nodes and the
        // normal / tangent vectors to the element
        el_pt->rotated_boundary_helper_pt()->set_nodal_boundary_parametrisation(
          boundary_node,
          boundary_coordinate_of_node,
          Parameters::Parametric_curve_pt[b]);
      }
    }
  }
} // end rotate_edge_degrees_of_freedom


//==start_of_pin_all_in_plane_displacements=====================================
/// Pin the in-plane displacements
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_all_in_plane_displacements()
{
  unsigned nnode = Bulk_mesh_pt->nnode();
  for (unsigned inode = 0; inode < nnode; inode++)
  {
    Bulk_mesh_pt->node_pt(inode)->pin(0);
    Bulk_mesh_pt->node_pt(inode)->set_value(0, 0.0);
    Bulk_mesh_pt->node_pt(inode)->pin(1);
    Bulk_mesh_pt->node_pt(inode)->set_value(1, 0.0);
  }
}


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 5;

  sprintf(filename, "RESLT/soln%i-%f.dat", Doc_info.number(), Element_area);
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();

  // Doc error and return of the square of the L2 error
  //---------------------------------------------------
  // double error,norm,dummy_error,zero_norm;
  double dummy_error, zero_norm;
  sprintf(filename, "RESLT/error%i-%f.dat", Doc_info.number(), Element_area);
  some_file.open(filename);

  Bulk_mesh_pt->compute_error(
    some_file, Parameters::dummy_exact_w, dummy_error, zero_norm);
  some_file.close();

  // Doc L2 error and norm of solution
  oomph_info << "Absolute norm of computed solution: " << sqrt(dummy_error)
             << std::endl;

  oomph_info << "Norm of computed solution: " << sqrt(zero_norm) << std::endl;

  // Find the solution at r=0
  //   // ----------------------
  MeshAsGeomObject* Mesh_as_geom_obj_pt = new MeshAsGeomObject(Bulk_mesh_pt);
  Vector<double> s(2);
  GeomObject* geom_obj_pt = 0;
  Vector<double> r(2, 0.0);
  Mesh_as_geom_obj_pt->locate_zeta(r, geom_obj_pt, s);
  // Compute the interpolated displacement vector
  Vector<double> u_0(12, 0.0);
  u_0 =
    dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_foeppl_von_karman(s);

  oomph_info << "w in the middle: " << std::setprecision(15) << u_0[0]
             << std::endl;

  Trace_file << u_0[0] << '\n';

  // Doc error and return of the square of the L2 error
  //---------------------------------------------------
  sprintf(filename, "RESLT/L2-norm%i-%f.dat", Doc_info.number(), Element_area);
  some_file.open(filename);

  some_file << "### L2 Norm\n";
  some_file << "##  Format: err^2 norm^2 \n";
  // Print error in prescribed format
  some_file << dummy_error << " " << zero_norm << "\n";
  some_file.close();

  // Increment the doc_info number
  Doc_info.number()++;

  // Clean up
  delete Mesh_as_geom_obj_pt;
} // end of doc


//=======start_of_main========================================
/// Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char** argv)
{
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Define possible command line arguments and parse the ones that
  // were actually specified
  // Directory for solution
  string output_dir = "RESLT";
  CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

  // Opening angle
  CommandLineArgs::specify_command_line_flag("--alpha", &Parameters::Alpha);

  // Poisson Ratio
  CommandLineArgs::specify_command_line_flag("--nu", &Parameters::Nu);

  // Applied Pressure
  CommandLineArgs::specify_command_line_flag("--eta", &Parameters::Eta);

  // Element Area (no larger element than 0.09)
  double element_area = 1.0;
  CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
    element_area);

  // Set up some problem paramters
  problem.max_residuals() = 1e3;
  problem.max_newton_iterations() = 20;

  // Pre mess doc
  problem.doc_solution();

  // Do the newton solve
  problem.steady_newton_solve();

  // Document
  problem.doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << "Solution number (" << problem.Doc_info.number() - 1 << ")"
             << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;
  // Print success
  oomph_info << "Exiting Normally\n";
} // End of main
