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
using MathematicalConstants::Pi;


//===========================================================================
/// Namespace for problem parameters
//===========================================================================
namespace Parameters
{
  /// Square edge length
  double L = 1.0;

  /// Angle of rotation (radians obviously)
  double Alpha = Pi / 5.0;

  /// Plate vertices. For the unrotated square, these coincide with:
  ///     L/2*{(1,1), (-1,1), (-1,-1), (1,-1)}
  Vector<Vector<double>> Vertices = {
    {L / 2.0 * cos(Alpha) - L / 2.0 * sin(Alpha),
     L / 2.0 * sin(Alpha) + L / 2.0 * cos(Alpha)},
    {-L / 2.0 * cos(Alpha) - L / 2.0 * sin(Alpha),
     -L / 2.0 * sin(Alpha) + L / 2.0 * cos(Alpha)},
    {-L / 2.0 * cos(Alpha) + L / 2.0 * sin(Alpha),
     -L / 2.0 * sin(Alpha) - L / 2.0 * cos(Alpha)},
    {L / 2.0 * cos(Alpha) + L / 2.0 * sin(Alpha),
     L / 2.0 * sin(Alpha) - L / 2.0 * cos(Alpha)}};

  /// The plate thickness
  double Thickness = 0.01;

  /// Poisson ratio
  double Nu = 0.5;

  /// FvK parameter (slightly dangerous; should really
  /// be computed as a dependent parameter to accommodate
  /// changes in Nu and Thickness.
  double Eta = 10000.0; // 12.0*(1.0-Nu*Nu)/(Thickness*Thickness);

  /// Magnitude of pressure
  double P_mag = 0.0;

  /// Magnitude of shear stress
  double T_mag = 0.0;

  /// Element size
  double Element_area = 0.01;

  /// Order of the polynomial boundary interpolation
  unsigned Boundary_order = 5;

  // ---- Parametric boundaries ------------------------------------------------

  /// Straight edge 0
  CurvilineLine Edge_0(Vertices[0], Vertices[1]);

  /// Straight edge 1
  CurvilineLine Edge_1(Vertices[1], Vertices[2]);

  /// Straight edge 2
  CurvilineLine Edge_2(Vertices[2], Vertices[3]);

  /// Straight edge 3
  CurvilineLine Edge_3(Vertices[3], Vertices[0]);

  /// Vector container of addresses for iterating over the edges
  Vector<CurvilineGeomObject*> Curviline_edge_pt = {
    &Edge_0, &Edge_1, &Edge_2, &Edge_3};

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
  /// Helper function to specify boundary conditions (here all homogeneous)
  /// as a function of both coordinates (This is convenient for this problem;
  /// other interfaces that specify boundary conditions in terms of
  /// boundary coordinate exist).
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

    // Update the corner constraintes based on boundary conditions
    unsigned n_el = Constraint_mesh_pt->nelement();
    for (unsigned i_el = 0; i_el < n_el; i_el++)
    {
      dynamic_cast<DuplicateNodeConstraintElement*>(
        Constraint_mesh_pt->element_pt(i_el))
        ->validate_and_pin_redundant_constraints();
    }

    // Reassign the equation numbers
    assign_eqn_numbers();

  } // End make_linear()


  /// Doc the solution
  void doc_solution(const std::string& comment = "");

  // [zdec] This doesn't work, dynamic cast always fails -- returns 0
  // /// Overloaded version of the problem's access function to
  // /// the mesh. Recasts the pointer to the base Mesh object to
  // /// the actual mesh type.
  // TriangleMesh<ELEMENT>* mesh_pt()
  // {
  //   oomph_info << Problem::mesh_pt() << std::endl;
  //   oomph_info << dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt())
  //   << std::endl; return dynamic_cast<TriangleMesh<ELEMENT>*>
  //   (Problem::mesh_pt());
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

  /// Duplicate corner nodes and create constraint elements at those corners
  void duplicate_corner_nodes();

  /// Helper function to (re-)set boundary condition
  /// and complete the build of all elements
  void complete_problem_setup();

  /// Loop over all curved edges, then loop over elements and upgrade
  /// them to be curved elements
  void upgrade_edge_elements_to_curve(const unsigned& b);

  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in
  /// Parameters
  void rotate_edge_degrees_of_freedom(Mesh* const& bulk_mesh_pt);

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

  /// Pointer to constraint mesh
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
  Doc_info.number() = 0;

  // Build the mesh
  build_mesh();

  // Curved Edge upgrade
  upgrade_edge_elements_to_curve(Boundary_b_bnum);
  upgrade_edge_elements_to_curve(Boundary_r_bnum);
  upgrade_edge_elements_to_curve(Boundary_t_bnum);
  upgrade_edge_elements_to_curve(Boundary_l_bnum);


  // Rotate degrees of freedom
  rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

  // Complete problem setup
  complete_problem_setup();


  // Output parameters
  oomph_info << "Problem parameters:\n"
             << "L            " << Parameters::L << std::endl
             << "Alpha        " << Parameters::Alpha << std::endl
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
  // Get the vertices from the parameters
  Vector<double> vertex0 = Parameters::Vertices[0];
  Vector<double> vertex1 = Parameters::Vertices[1];
  Vector<double> vertex2 = Parameters::Vertices[2];
  Vector<double> vertex3 = Parameters::Vertices[3];

  // Declare the edges...
  Vector<Vector<double>> edge0(2, Vector<double>(2, 0.0));
  Vector<Vector<double>> edge1(2, Vector<double>(2, 0.0));
  Vector<Vector<double>> edge2(2, Vector<double>(2, 0.0));
  Vector<Vector<double>> edge3(2, Vector<double>(2, 0.0));

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

  // Split elements that have two boundary edges
  TimeStepper* time_stepper_pt = Bulk_mesh_pt->Time_stepper_pt;
  Bulk_mesh_pt->template split_elements_with_multiple_boundary_edges<ELEMENT>(
    time_stepper_pt);

  // Create the empty constraint element mesh
  Constraint_mesh_pt = new Mesh();

  // Add extra nodes at boundaries and constrain the dofs there.
  duplicate_corner_nodes();

  // Add submesh to problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Constraint_mesh_pt);

  // Combine submeshes into a single Mesh (over the top; could just have
  // assigned bulk mesh directly.
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
  unsigned n_bound = 4;
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
      Parameters::Curviline_edge_pt[i_bound];
    CurvilineGeomObject* right_parametrisation_pt =
      Parameters::Curviline_edge_pt[ip1_bound];

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

    // Assign the parameter function pointers for the element
    el_pt->nu_pt() = &Parameters::Nu;
    el_pt->eta_pt() = &Parameters::Eta;
  }
  // Set the boundary conditions
  apply_boundary_conditions();

} // end of complete


//==============================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to
/// set up physical boundary conditions on a curved boundary with Hermite type
/// dofs. For example if we know w(n,t) = f(t) (where n and t are the normal and
/// tangent to a boundary) we ALSO know dw/dt and d2w/dt2. NB no rotation is
/// needed if the edges are completely free! begin
/// rotate_edge_degrees_of_freedom
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::rotate_edge_degrees_of_freedom(
  Mesh* const& bulk_mesh_pt)
{
  // Get the number of boundaries
  unsigned n_bound = 4;

  // Loop over the bulk elements
  unsigned n_element = Bulk_mesh_pt->nelement();
  for (unsigned e = 0; e < n_element; e++)
  {
    // Get pointer to bulk element adjacent
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

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
          double coord = Parameters::Curviline_edge_pt[b]->get_zeta(
            el_pt->node_pt(n)->position());
          boundary_coordinate_of_node.push_back(coord);
        }
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if (!boundary_node.empty())
      {
        // Rotate the nodes by passing the index of the nodes and the
        // normal / tangent vectors to the element
        el_pt->rotated_boundary_helper_pt()->set_nodal_boundary_parametrisation(
          boundary_node,
          boundary_coordinate_of_node,
          Parameters::Curviline_edge_pt[b]);
      }
    }
  }
} // end rotate_edge_degrees_of_freedom


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
  // So, if the function fix_in_plane_displacement_dof(idof, b, fct_pt) is
  // called with idof=1 and b=3, say, the y-in-plane displacement is pinned for
  // all the element's nodes (if any) that are located on mesh boundary 3. The
  // value of the y-in-plane displacement is set to whatever the function
  // pointed to by fct_pt computes when evaluated at the nodal coordinate.
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
  const Vector<unsigned> pin_ux_and_uy_pinned_dof{0, 1};


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
  const Vector<unsigned> pinned_edge_pinned_dof{0, 2, 5};

  // Case: The plate is sliding (w left free, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) dw/dn and d^2w/dndt
  const Vector<unsigned> sliding_clamp_pinned_dof{1, 4};

  // Case: The plate is clamped (w given, dw/dn given) along a boundary.
  // We therefore have to pin (and assign values for) w, dw/dn,
  // dw/dt, d^2w/dndt and d^2w/dt^2
  const Vector<unsigned> fully_clamped_pinned_dof{0, 1, 2, 4, 5};


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
  pinned_w_dofs[0] = fully_clamped_pinned_dof;
  pinned_w_dofs[1] = fully_clamped_pinned_dof;
  pinned_w_dofs[2] = fully_clamped_pinned_dof;
  pinned_w_dofs[3] = fully_clamped_pinned_dof;


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

      // Pin in-plane dofs (enumerated as explained above) for
      // all nodes on boundary b. Here we're applying homogeneous
      // BCs so all pinned values are simply set to zero.
      for (unsigned i = 0; i < n_pinned_u_dofs; i++)
      {
        unsigned idof_to_be_pinned = pinned_u_dofs[b][i];
        el_pt->fix_in_plane_displacement_dof(
          idof_to_be_pinned, b, Parameters::get_null_fct);
      }
      // Pin out-of-plane dofs (enumerated as explained above) for all
      // nodes on boundary b. Here we're applying homogeneous BCs so
      // all pinned values are simply set to zero.
      for (unsigned i = 0; i < n_pinned_w_dofs; i++)
      {
        unsigned idof_to_be_pinned = pinned_w_dofs[b][i];
        el_pt->fix_out_of_plane_displacement_dof(
          idof_to_be_pinned, b, Parameters::get_null_fct);
      }
    } // end for loop over elements on b
  } // end for loop over boundaries

  // Update the corner constraintes based on boundary conditions
  unsigned n_el = Constraint_mesh_pt->nelement();
  for (unsigned i_el = 0; i_el < n_el; i_el++)
  {
    dynamic_cast<DuplicateNodeConstraintElement*>(
      Constraint_mesh_pt->element_pt(i_el))
      ->validate_and_pin_redundant_constraints();
  }

  // Assign the equation numbers
  assign_eqn_numbers();

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
  const unsigned& ibound)
{
  // Loop over the bulk elements adjacent to boundary ibound
  const unsigned n_els = Bulk_mesh_pt->nboundary_element(ibound);
  for (unsigned e = 0; e < n_els; e++)
  {
    // Get pointer to bulk element adjacent to b
    ELEMENT* bulk_el_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound, e));

    // hierher what is that? why "My"?
    // Initialise enum for the curved edge
    MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

    // Loop over all (three) nodes of the element and record boundary nodes
    unsigned index_of_interior_node = 3;
    unsigned nnode_on_neither_boundary = 0;
    const unsigned nnode = 3;


    // hierher what does this comment mean?
    // Fill in vertices' positions (this step could be moved inside the
    // curveable Bell element)
    Vector<Vector<double>> xn(nnode, Vector<double>(2, 0.0));
    for (unsigned n = 0; n < nnode; ++n)
    {
      Node* nod_pt = bulk_el_pt->node_pt(n);
      xn[n][0] = nod_pt->x(0);
      xn[n][1] = nod_pt->x(1);

      // Check if it is on the outer boundaries
      if (!(nod_pt->is_on_boundary(ibound)))
      {
        index_of_interior_node = n;
        ++nnode_on_neither_boundary;
      }
    } // end record boundary nodes


    // s at the next (cyclic) node after interior
    const double s_ubar = Parameters::Curviline_edge_pt[ibound]->get_zeta(
      xn[(index_of_interior_node + 1) % 3]);

    // s at the previous (cyclic) node before interior
    const double s_obar = Parameters::Curviline_edge_pt[ibound]->get_zeta(
      xn[(index_of_interior_node + 2) % 3]);

    // Assign edge case
    edge = static_cast<MyC1CurvedElements::Edge>(index_of_interior_node);

    // Check nnode_on_neither_boundary
    if (nnode_on_neither_boundary == 0)
    {
      throw OomphLibError(
        "No interior nodes. One node per CurvedElement must be interior.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
    else if (nnode_on_neither_boundary > 1)
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
        "Decreasing parametric coordinate. Parametric coordinate must increase "
        "as the edge is traversed anti-clockwise.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    } // end checks

    // Upgrade it
    bulk_el_pt->upgrade_element_to_curved(edge,
                                          s_ubar,
                                          s_obar,
                                          Parameters::Curviline_edge_pt[ibound],
                                          Parameters::Boundary_order);
  }
} // end_upgrade_elements


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

  Trace_file << Parameters::P_mag << " " << Parameters::T_mag << " "
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

  // Make problem linear
  problem.make_linear();

  // Set pressure
  Parameters::P_mag = 10.0;
  // Set Poisson ratio
  Parameters::Nu = 0.5;
  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

  // Set Poisson ratio
  Parameters::Nu = 0.0;
  // Solve the system
  problem.newton_solve();
  // Document the current solution
  problem.doc_solution();

  // Test the rotated dofs


} // End of main