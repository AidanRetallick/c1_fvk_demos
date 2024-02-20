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

//                      OUTLINE OF PROBLEM CONSTRUCTION
// The basic constuction is much the same as the usual order of things in a
// problem. Underneath is the order of actions (with stars next to non actions
// that are unique to these types of problems).
// 1.  Setup mesh parameters
// 2.  Build the mesh

// hierher: combine 3 and 4?
//   aidan: these are different steps (unupgraded elements may need rotating)

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


  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////


  /// \short Specialisation of CurvilineGeomObject for half a circle.
  class MatthiasCurvilineCircleTop : public CurvilineGeomObject
  {
  public:
    /// \short Constructor: Pass dimension of geometric object (# of Eulerian
    /// coords = # of Lagrangian coords; no time history available/needed)
    MatthiasCurvilineCircleTop()
      : CurvilineGeomObject(), Radius(1.0), Clockwise_zeta(false)
    {
    }

    /// \short Constructor: Pass dimension of geometric object (# of Eulerian
    /// coords = # of Lagrangian coords; no time history available/needed)
    MatthiasCurvilineCircleTop(const double& radius, const bool& clockwise_zeta)
      : CurvilineGeomObject(), Radius(radius), Clockwise_zeta(clockwise_zeta)
    {
    }
    /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
    /// and pointer to time-stepper which is used to handle the
    /// position at previous timesteps and allows the evaluation
    /// of veloc/acceleration etc. in cases where the GeomData
    /// varies with time.
    MatthiasCurvilineCircleTop(TimeStepper* time_stepper_pt)
      : CurvilineGeomObject(time_stepper_pt)
    {
    }

    /// Broken copy constructor
    MatthiasCurvilineCircleTop(const MatthiasCurvilineCircleTop& dummy)
    {
      BrokenCopy::broken_copy("MatthiasCurvilineCircleTop");
    }

    /// Broken assignment operator
    void operator=(const MatthiasCurvilineCircleTop&)
    {
      BrokenCopy::broken_assign("MatthiasCurvilineCircleTop");
    }

    /// (Empty) destructor
    virtual ~MatthiasCurvilineCircleTop() {}

    /// \short Position Vector w.r.t. to zeta:
    virtual void position(const Vector<double>& zeta, Vector<double>& r) const
    {
      r[0] = -Radius * std::sin(zeta[0]);
      r[1] = Radius * std::cos(zeta[0]);
      // Zeta -> - Zeta
      if (Clockwise_zeta)
      {
        r[0] *= -1;
      }
    }

    /// \short Derivative of position Vector w.r.t. to zeta:
    virtual void dposition(const Vector<double>& zeta,
                           Vector<double>& drdzeta) const
    {
      drdzeta[0] = -Radius * std::cos(zeta[0]);
      drdzeta[1] = -Radius * std::sin(zeta[0]);
      // Zeta -> - Zeta
      if (Clockwise_zeta)
      {
        drdzeta[0] *= -1;
      }
    }


    /// \short 2nd derivative of position Vector w.r.t. to coordinates:
    /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
    /// ddrdzeta(alpha,beta,i).
    /// Evaluated at current time.
    virtual void d2position(const Vector<double>& zeta,
                            Vector<double>& drdzeta) const
    {
      drdzeta[0] = Radius * std::sin(zeta[0]);
      drdzeta[1] = -Radius * std::cos(zeta[0]);
      // Zeta -> - Zeta
      if (Clockwise_zeta)
      {
        drdzeta[0] *= -1;
      }
    }

    /// Get s from x for part 0 of the boundary (inverse mapping - for
    /// convenience)
    double get_zeta(const Vector<double>& x) const
    {
      // The arc length (parametric parameter) for the upper semi circular arc
      return (Clockwise_zeta ? atan2(x[0], x[1]) : atan2(-x[0], x[1]));
    }

  private:
    double Radius;
    bool Clockwise_zeta;
  };


  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////


  /// \short Specialisation of CurvilineGeomObject for half a circle.
  class MatthiasCurvilineCircleBottom : public CurvilineGeomObject
  {
  public:
    /// \short Constructor: Pass dimension of geometric object (# of Eulerian
    /// coords = # of Lagrangian coords; no time history available/needed)
    MatthiasCurvilineCircleBottom()
      : CurvilineGeomObject(), Radius(1.0), Clockwise_zeta(false)
    {
    }

    /// \short Constructor: Pass dimension of geometric object (# of Eulerian
    /// coords = # of Lagrangian coords; no time history available/needed)
    MatthiasCurvilineCircleBottom(const double& radius,
                                  const bool& clockwise_zeta)
      : CurvilineGeomObject(), Radius(radius), Clockwise_zeta(clockwise_zeta)
    {
    }

    /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
    /// and pointer to time-stepper which is used to handle the
    /// position at previous timesteps and allows the evaluation
    /// of veloc/acceleration etc. in cases where the GeomData
    /// varies with time.
    MatthiasCurvilineCircleBottom(TimeStepper* time_stepper_pt)
      : CurvilineGeomObject(time_stepper_pt)
    {
    }

    /// Broken copy constructor
    MatthiasCurvilineCircleBottom(const MatthiasCurvilineCircleBottom& dummy)
    {
      BrokenCopy::broken_copy("MatthiasCurvilineCircleBottom");
    }

    /// Broken assignment operator
    void operator=(const MatthiasCurvilineCircleBottom&)
    {
      BrokenCopy::broken_assign("MatthiasCurvilineCircleBottom");
    }

    /// (Empty) destructor
    virtual ~MatthiasCurvilineCircleBottom() {}

    /// \short Position Vector w.r.t. to zeta:
    virtual void position(const Vector<double>& zeta, Vector<double>& r) const
    {
      r[0] = Radius * std::sin(zeta[0]);
      r[1] = -Radius * std::cos(zeta[0]);
      // Zeta -> - Zeta
      if (Clockwise_zeta)
      {
        r[1] *= -1;
      }
    }

    /// \short Derivative of position Vector w.r.t. to zeta:
    virtual void dposition(const Vector<double>& zeta,
                           Vector<double>& drdzeta) const
    {
      drdzeta[0] = Radius * std::cos(zeta[0]);
      drdzeta[1] = Radius * std::sin(zeta[0]);
      if (Clockwise_zeta)
      {
        drdzeta[1] *= -1;
      }
    }


    /// \short 2nd derivative of position Vector w.r.t. to coordinates:
    /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
    /// ddrdzeta(alpha,beta,i).
    /// Evaluated at current time.
    virtual void d2position(const Vector<double>& zeta,
                            Vector<double>& drdzeta) const
    {
      drdzeta[0] = -Radius * std::sin(zeta[0]);
      drdzeta[1] = Radius * std::cos(zeta[0]);
      if (Clockwise_zeta)
      {
        drdzeta[1] *= -1;
      }
    }

    /// Get s from x for part 0 of the boundary (inverse mapping - for
    /// convenience)
    double get_zeta(const Vector<double>& x) const
    {
      // The arc length (parametric parameter) for the upper semi circular arc
      return (Clockwise_zeta ? atan2(x[0], x[1]) : atan2(x[0], -x[1]));
    }

  private:
    double Radius;
    bool Clockwise_zeta;
  };

  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // hierher
  typedef FoepplVonKarmanC1CurvableBellElement<4> NON_WRAPPED_ELEMENT;


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
      oomph_info << "# of internal Data objects before: "
                 << this->ninternal_data() << std::endl;
    }

    /// Destructor (empty)
    ~FvKPointForceAndSourceElement() {}

    /// Set local coordinate and point force and_torque
    void setup(const Vector<double>& s_point_force_and_torque)
    {
      S_point_force_and_torque = s_point_force_and_torque;
    }


    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0 hierher why no
      // arg? using a dummy matrix argument
      ELEMENT::fill_in_generic_residual_contribution_foeppl_von_karman(
        residuals, GeneralisedElement::Dummy_matrix, 0);

      // fill_in_contribution_to_residuals(residuals);

      // Add point force_and_torque contribution
      fill_in_point_force_and_torque_contribution_to_residuals(residuals);
    }


    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1 hierher why no arg?
      //  ELEMENT::fill_in_contribution_to_jacobian(residuals,
      //                                            jacobian);
      ELEMENT::fill_in_generic_residual_contribution_foeppl_von_karman(
        residuals, jacobian, 1);

      // Add point force_and_torque contribution
      fill_in_point_force_and_torque_contribution_to_residuals(residuals);
    }


  private:
    /// Add the point force_and_torque contribution to the residual vector
    void fill_in_point_force_and_torque_contribution_to_residuals(
      Vector<double>& residuals)
    {
      // No further action
      if (S_point_force_and_torque.size() == 0)
      {
        oomph_info << "bailing" << std::endl;
        return;
      }


      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memory for the shape/test functions
      Shape psi(n_node);

      // //Integers to store the local equation and unknown numbers
      // int local_eqn_real=0;
      // int local_eqn_imag=0;

      // Get shape/test fcts
      this->shape(S_point_force_and_torque, psi);

      //  // Assemble residuals
      //  //--------------------------------

      //  // Loop over the test functions
      //  for(unsigned l=0;l<n_node;l++)
      //   {
      //    // first, compute the real part contribution
      //    //-------------------------------------------

      //    //Get the local equation
      //    local_eqn_real =
      //    this->nodal_local_eqn(l,this->u_index_helmholtz().real());

      //    /*IF it's not a boundary condition*/
      //    if(local_eqn_real >= 0)
      //     {
      //      residuals[local_eqn_real] +=
      //      Point_force_and_torque_magnitude.real()*psi(l);
      //     }

      //    // Second, compute the imaginary part contribution
      //    //------------------------------------------------

      //    //Get the local equation
      //    local_eqn_imag =
      //    this->nodal_local_eqn(l,this->u_index_helmholtz().imag());

      //    /*IF it's not a boundary condition*/
      //    if(local_eqn_imag >= 0)
      //     {
      //      // Add body force/force_and_torque term and Helmholtz bit
      //      residuals[local_eqn_imag] +=
      //      Point_force_and_torque_magnitude.imag()*psi(l);
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
  class FaceGeometry<FvKPointForceAndSourceElement<ELEMENT>>
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
  class FaceGeometry<FaceGeometry<FvKPointForceAndSourceElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };


} // namespace oomph

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////


//========================================================================
/// Namespace for problem parameters
//========================================================================
namespace Parameters
{


  // Enumeration of cases
  enum
  {
    Clamped_validation,
    Axisymmetric_shear_buckling,
    Nonaxisymmetric_shear_buckling
  };

  /// Which case are we doing
  unsigned Problem_case = Clamped_validation; // Nonaxisymmetric_shear_buckling;
                                              // // Axisymmetric_shear_buckling;

  // Ellipse half axis
  double A = 1.0;

  // Other ellipse half axis
  double B = 1.0;

  /// Poisson ratio
  double Nu = 0.5;

  // /// Membrane coupling coefficient (this should really be computed
  // /// as a dependent parameter...)
  // double Eta = 12.0*(1.0-Nu*Nu)/(Thickness*Thickness);

  // FvK parameter
  double Eta = 2.39e6;

  /// Pressure magnitude
  double P_mag = 0.0;

  /// In-plane traction magnitude
  double T_mag = 0.0;

  /// Order of boundary interpolation
  unsigned Boundary_order = 3;

  // hierher what are these objects? Shouldn't they be
  // used in the mesh generatino too; surely they encode the
  // same information.

  /// Parametric curve for the upper half boundary
  CurvilineEllipseTop Parametric_curve_top(A, B);

  /// Parametric curve for the lower half boundary
  CurvilineEllipseBottom Parametric_curve_bottom(A, B);

  /// Vector of parametric boundaries
  Vector<CurvilineGeomObject*> Parametric_curve_pt = {&Parametric_curve_top,
                                                      &Parametric_curve_bottom};


  /// Pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
    // Constant pressure for validation case
    if (Parameters::Problem_case == Parameters::Clamped_validation)
    {
      pressure = P_mag;
    }
    // Parabolic pressure distribution with zero mean
    else if (Parameters::Problem_case ==
             Parameters::Axisymmetric_shear_buckling)
    {
      pressure = P_mag * (0.25 - x[0] * x[0] - x[1] * x[1]);
    }
    // Parabolic pressure distribution with zero mean
    else if (Parameters::Problem_case ==
             Parameters::Nonaxisymmetric_shear_buckling)
    {
      pressure = P_mag * (0.25 - x[0] * x[0] - x[1] * x[1]);
    }
    else
    {
      throw OomphLibError("Unexpected problem setup",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }


  /// In plane forcing (shear stress) depending on the position (x,y)
  void get_in_plane_force(const Vector<double>& x, Vector<double>& tau)
  {
    // Zero shear stress for validation case
    if (Parameters::Problem_case == Parameters::Clamped_validation)
    {
      tau[0] = 0.0;
      tau[1] = 0.0;
    }
    // Self balancing purely radially outward shear stress
    else if (Parameters::Problem_case ==
             Parameters::Axisymmetric_shear_buckling)
    {
      double phi = atan2(x[1], x[0]);
      double r_squared = x[0] * x[0] + x[1] * x[1];
      tau[0] = T_mag * r_squared * cos(phi);
      tau[1] = T_mag * r_squared * sin(phi);
    }
    // Self-balancing y shear stress over disk:
    else if (Parameters::Problem_case ==
             Parameters::Nonaxisymmetric_shear_buckling)
    {
      //   tau_y := 1/4 - y^2;
      //
      //   resultant:=int(int(subs(y = r*sin(phi), tau_y)*r, phi = 0 .. 2*Pi), r
      //   = 0 .. 1);
      tau[0] = 0.0;
      tau[1] = T_mag * (0.25 - x[1] * x[1]);
    }
    else
    {
      throw OomphLibError("Unexpected problem setup",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }

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

} // namespace Parameters

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
    delete Outer_boundary_ellipse_pt;
    delete Outer_curvilinear_boundary_pt[0];
    delete Outer_curvilinear_boundary_pt[1];
    delete Inner_open_boundaries_pt[0];
    delete Inner_open_boundaries_pt[1];
    delete Boundary2_pt;
    delete Boundary3_pt;
  };

  /// Update after solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs before solve: empty
  void actions_before_newton_solve() {}

  /// [zdec] temp
  void actions_before_newton_step()
  {
    // Filenames
    char res_filename[100];
    char jac_filename[100];
    sprintf(res_filename, "res_%i_%i", Doc_info.number(), Nnewton_iter_taken);
    sprintf(jac_filename, "jac_%i_%i", Doc_info.number(), Nnewton_iter_taken);
    // Get the jacobian
    LinearAlgebraDistribution* dist = this->dof_distribution_pt();
    DoubleVector res(dist, 0.0);
    CRDoubleMatrix jac(dist);
    get_jacobian(res, jac);
    res.output(res_filename);
    jac.sparse_indexed_output(jac_filename);
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

  /// Overloaded version of the problem's access function to
  /// the mesh. Recasts the pointer to the base Mesh object to
  /// the actual mesh type.
  TriangleMesh<ELEMENT>* mesh_pt()
  {
    return Bulk_mesh_pt; // dynamic_cast<TriangleMesh<ELEMENT>*>
                         // (Problem::mesh_pt());
  }


private:
  /// Pin all displacements and rotation (dofs 0-4) at the centre
  void pin_all_displacements_and_rotation_at_centre_node();

  /// Setup and build the mesh
  void build_mesh();

  /// Duplicate corner nodes and create constraint elements at those corners
  void duplicate_corner_nodes();

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

  /// Loop over all curved edges, then loop over elements and upgrade
  /// them to be curved elements
  void upgrade_edge_elements_to_curve(const unsigned& b);

  /// Loop over all edge elements and rotate the Hermite degrees of freedom
  /// to be in the directions of the two in-plane vectors specified in
  /// Parameters
  void rotate_edge_degrees_of_freedom();

  /// Delete traction elements and wipe the surface mesh
  void delete_traction_elements(Mesh* const& surface_mesh_pt);

  /// Pointer to "bulk" mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Enumeration to keep track of boundary ids
  enum
  {
    Outer_boundary0 = 0,
    Outer_boundary1 = 1,
    Inner_boundary0 = 2,
    Inner_boundary1 = 3
  };

  /// Target element area
  double Element_area;

  /// Pointer to constraint mesh
  Mesh* Constraint_mesh_pt;

  /// Doc info object for labeling output
  DocInfo Doc_info;

  /// Outer boundary Geom Object
  Ellipse* Outer_boundary_ellipse_pt;

  /// The outer curves
  Vector<TriangleMeshCurveSection*> Outer_curvilinear_boundary_pt;

  /// The Internal curves
  Vector<TriangleMeshOpenCurve*> Inner_open_boundaries_pt;

  /// The close outer boundary
  TriangleMeshClosedCurve* Outer_boundary_pt;

  /// The first of the internal boundaries
  TriangleMeshPolyLine* Boundary2_pt;

  /// The second of the internal boundaries
  TriangleMeshPolyLine* Boundary3_pt;

}; // end_of_problem_class


//======================================================================
/// Constructor definition
//======================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(
  const double& element_area)
  : Element_area(element_area)
{
  // Build the mesh
  build_mesh();

  // Curved Edge upgrade
  upgrade_edge_elements_to_curve(Outer_boundary0);
  upgrade_edge_elements_to_curve(Outer_boundary1);

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
  oomph_info << "Number of equations: " << assign_eqn_numbers() << '\n';

} // end Constructor


//======================================================================
/// Set up and build the mesh
//======================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::build_mesh()
{
  Vector<double> zeta(1);
  Vector<double> posn(2);

  // Outer boundary
  //--------------

  double A = Parameters::A;
  double B = Parameters::B;
  Outer_boundary_ellipse_pt = new Ellipse(A, B);

  // First bit
  double zeta_start = 0.0;
  double zeta_end = MathematicalConstants::Pi;
  unsigned nsegment = (int)(MathematicalConstants::Pi / sqrt(Element_area));

  Outer_curvilinear_boundary_pt.resize(2);
  Outer_curvilinear_boundary_pt[0] = new TriangleMeshCurviLine(
    Outer_boundary_ellipse_pt, zeta_start, zeta_end, nsegment, Outer_boundary0);

  // Second bit
  zeta_start = MathematicalConstants::Pi;
  zeta_end = 2.0 * MathematicalConstants::Pi;
  nsegment = (int)(MathematicalConstants::Pi / sqrt(Element_area));
  Outer_curvilinear_boundary_pt[1] = new TriangleMeshCurviLine(
    Outer_boundary_ellipse_pt, zeta_start, zeta_end, nsegment, Outer_boundary1);

  // Combine
  Outer_boundary_pt =
    new TriangleMeshClosedCurve(Outer_curvilinear_boundary_pt);

  // Internal open boundaries
  //-------------------------
  // Total number of open curves in the domain
  unsigned n_open_curves = 2;
  // We want internal open curves
  Inner_open_boundaries_pt.resize(n_open_curves);

  // Internal bit - this means we can have a boundary which is just the centre
  // We start by creating the internal boundaries

  // Open curve 1
  Vector<Vector<double>> vertices(2, Vector<double>(2, 0.0));
  vertices[0][0] = -0.5;
  vertices[0][1] = 0.0;

  vertices[1][0] = 0.5;
  vertices[1][1] = 0.0;
  unsigned boundary_id = Inner_boundary0;

  Boundary2_pt = new TriangleMeshPolyLine(vertices, boundary_id);

  // Open Curve 2
  vertices[0][0] = 0.0;
  vertices[0][1] = -0.5;

  vertices[1][0] = 0.0;
  vertices[1][1] = 0.5;
  boundary_id = Inner_boundary1;

  Boundary3_pt = new TriangleMeshPolyLine(vertices, boundary_id);

  // Each internal open curve is defined by a vector of
  // TriangleMeshCurveSections
  Vector<TriangleMeshCurveSection*> internal_curve_section1_pt(1);
  internal_curve_section1_pt[0] = Boundary2_pt;

  Vector<TriangleMeshCurveSection*> internal_curve_section2_pt(1);
  internal_curve_section2_pt[0] = Boundary3_pt;

  // The open curve that defines this boundary
  Inner_open_boundaries_pt[0] =
    new TriangleMeshOpenCurve(internal_curve_section1_pt);

  Inner_open_boundaries_pt[1] =
    new TriangleMeshOpenCurve(internal_curve_section2_pt);


  // Create mesh parameters object
  TriangleMeshParameters mesh_parameters(Outer_boundary_pt);

  // Element area
  mesh_parameters.element_area() = Element_area;

  // Specify the internal open boundaries
  mesh_parameters.internal_open_curves_pt() = Inner_open_boundaries_pt;

  // Build an assign bulk mesh
  Bulk_mesh_pt = new TriangleMesh<ELEMENT>(mesh_parameters);

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
  unsigned n_bound = 2;
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
        // // [zdec] debug
        // oomph_info << "Found a corner node at " << std::endl << "  ("
        //            << node_pt->position(0) << "," << node_pt->position(1) <<
        //            ")"
        //            << std::endl;
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
    new_node_pt->x(0) = old_node_pt->x(0);
    new_node_pt->x(1) = old_node_pt->x(1);
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


//==start_of_pin_all_displacements_and_rotation_at_centre_node======================
/// pin all displacements and rotations in the centre
//==============================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<
  ELEMENT>::pin_all_displacements_and_rotation_at_centre_node()
{
  // Choose non-centre node on which we'll supress
  // the rigid body rotation around the z axis.
  double max_x_potentially_pinned_node = -DBL_MAX;
  Node* pinned_rotation_node_pt = 0;

  // Pin the node that is at the centre in the domain
  unsigned num_int_nod = Bulk_mesh_pt->nboundary_node(2);
  for (unsigned inod = 0; inod < num_int_nod; inod++)
  {
    // Get node point
    Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(2, inod);

    // Check which coordinate increases along this boundary
    // oomph_info << "node: "
    //            << nod_pt->x(0) << " "
    //            << nod_pt->x(1) << " "
    //            << std::endl;

    // Find the node with the largest x coordinate
    if (fabs(nod_pt->x(0)) > max_x_potentially_pinned_node)
    {
      max_x_potentially_pinned_node = fabs(nod_pt->x(0));
      pinned_rotation_node_pt = nod_pt;
    }


    // If the node is on the other internal boundary too
    if (nod_pt->is_on_boundary(3))
    {
      // Pin it! It's the centre of the domain!
      // In-plane dofs are always 0 and 1
      // Out of plane displacement is 2, x and y derivatives are 3 and 4.
      nod_pt->pin(0);
      nod_pt->set_value(0, 0.0);
      nod_pt->pin(1);
      nod_pt->set_value(1, 0.0);
      nod_pt->pin(2);
      nod_pt->set_value(2, 0.0);
      nod_pt->pin(3);
      nod_pt->set_value(3, 0.0);
      nod_pt->pin(4);
      nod_pt->set_value(4, 0.0);
    }
  }

  oomph_info << "rotation pinning node: " << pinned_rotation_node_pt->x(0)
             << " " << pinned_rotation_node_pt->x(1) << " " << std::endl;
  // Pin y displacement
  pinned_rotation_node_pt->pin(1);
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

    // Assign the parameter pointers for the element
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
  // Clamp it
  if (Parameters::Problem_case == Parameters::Clamped_validation)
  {
    // Set the boundary conditions
    unsigned nbound = 2;
    for (unsigned b = 0; b < nbound; b++)
    {
      const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
      for (unsigned e = 0; e < nb_element; e++)
      {
        // Get pointer to bulk element adjacent to b
        ELEMENT* el_pt =
          dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b, e));

        // A true clamp, so we set everything except the second normal to zero
        for (unsigned idof = 0; idof < 6; ++idof)
        {
          // Cannot set second normal derivative
          if (idof != 3)
          {
            el_pt->fix_out_of_plane_displacement_dof(
              idof, b, Parameters::get_null_fct);
          }
        }
      }
    }
  }
  // All other cases: simply pin and stop rotation via the centre
  else
  {
    pin_all_displacements_and_rotation_at_centre_node();
  }

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
  // Parametric curve describing the boundary
  // hierher code share with thing that defines the
  // boundaries of the mesh!
  CurvilineGeomObject* parametric_curve_pt = 0;
  switch (ibound)
  {
    case 0:
      parametric_curve_pt = &Parameters::Parametric_curve_top;
      break;

    case 1:
      parametric_curve_pt = &Parameters::Parametric_curve_bottom;
      break;

    default:
      throw OomphLibError("Unexpected boundary number.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
      break;
  } // end parametric curve switch


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
    // Fill in vertices' positions (this step should be moved inside the
    // curveable Bell element)
    Vector<Vector<double>> xn(nnode, Vector<double>(2, 0.0));
    for (unsigned n = 0; n < nnode; ++n)
    {
      Node* nod_pt = bulk_el_pt->node_pt(n);
      xn[n][0] = nod_pt->x(0);
      xn[n][1] = nod_pt->x(1);

      // Check if it is on the outer boundaries
      if (!(nod_pt->is_on_boundary(Outer_boundary0) ||
            nod_pt->is_on_boundary(Outer_boundary1)))
      {
        index_of_interior_node = n;
        ++nnode_on_neither_boundary;
      }
    } // end record boundary nodes


    // hierher: ouch! This seems to map (x,y) to zeta! This is at best possible
    // to within a tolerance. Needs a redesign!

    // aidan: Inverse exists for (x,y) in map image (the boundary), so shouldn't
    // this be fine for boundary nodes?

    // s at the next (cyclic) node after interior
    const double s_ubar =
      parametric_curve_pt->get_zeta(xn[(index_of_interior_node + 1) % 3]);

    // s at the previous (cyclic) node before interior
    const double s_obar =
      parametric_curve_pt->get_zeta(xn[(index_of_interior_node + 2) % 3]);

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

    // Upgrade it // hierher what is "3"?
    bulk_el_pt->upgrade_element_to_curved(
      edge, s_ubar, s_obar, parametric_curve_pt, Parameters::Boundary_order);
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
  // std::ofstream debugstream;
  // debugstream.open("test_file");
  // debugstream << "x y nx ny tx ty nx_x nx_y ny_x ny_y tx_x tx_y ty_x ty_y"
  //             << std::endl;
  // debugstream.close();

  // Get the number of boundaries
  unsigned n_bound = 2;
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
          double coord =
            Parametric_curve_pt[b]->get_zeta(el_pt->node_pt(n)->position());
          boundary_coordinate_of_node.push_back(coord);
        }
      }

      // If the element has nodes on the boundary, rotate the Hermite dofs
      if (!boundary_node.empty())
      {
        // Rotate the nodes by passing the index of the nodes and the
        // normal / tangent vectors to the element
        el_pt->rotated_boundary_helper_pt()->set_nodal_boundary_parametrisation(
          boundary_node, boundary_coordinate_of_node, Parametric_curve_pt[b]);
      }
    }
  }
} // end rotate_edge_degrees_of_freedom


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const std::string& comment)
{
  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 30;

  sprintf(
    filename, "%s/soln%i.dat", Doc_info.directory().c_str(), Doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file, npts);
  some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" << comment << "\"\n";
  some_file.close();


  // Find the solution at r=0
  // ----------------------

  // hierher precompute
  MeshAsGeomObject Mesh_as_geom_obj(Bulk_mesh_pt);
  Vector<double> s(2);
  GeomObject* geom_obj_pt = 0;
  Vector<double> r(2, 0.0);
  Mesh_as_geom_obj.locate_zeta(r, geom_obj_pt, s);

  // Compute the interpolated displacement vector
  Vector<double> u_0(3, 0.0);
  u_0 = dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_fvk_disp(s);

  oomph_info << "w in the middle: " << std::setprecision(15) << u_0[2]
             << std::endl;

  Trace_file << Parameters::P_mag << " " << u_0[0] << '\n';

  // Increment the doc_info number
  Doc_info.number()++;

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


  // Clamped boundary conditions?
  CommandLineArgs::specify_command_line_flag("--use_clamped_bc");

  // Poisson Ratio
  CommandLineArgs::specify_command_line_flag("--nu", &Parameters::Nu);

  // Applied Pressure
  CommandLineArgs::specify_command_line_flag("--p", &Parameters::P_mag);

  // FvK prameter
  CommandLineArgs::specify_command_line_flag("--eta", &Parameters::Eta);

  // Element Area
  double element_area = 0.09;
  CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

  // Order of the boundary interpolation
  CommandLineArgs::specify_command_line_flag("--boundary_order",
                                             &Parameters::Boundary_order);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();


  // Constant pressure for validation case
  if (CommandLineArgs::command_line_flag_has_been_set("--use_clamped_bc"))
  {
    Parameters::Problem_case = Parameters::Clamped_validation;
  }

  // Build problem
  // UnstructuredFvKProblem<NON_WRAPPED_ELEMENT
  // UnstructuredFvKProblem<FvKPointForceAndSourceElement<NON_WRAPPED_ELEMENT>>
  UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>> problem(
    element_area);


  // Control parameters
  double dp_mag = 0.000001;
  double dt_mag = 0.000001;
  unsigned nstep = 1000;

  // Which case are we doing
  if (Parameters::Problem_case == Parameters::Clamped_validation)
  {
    nstep = 1;
    dp_mag = 0.01;
    dt_mag = 0.000001;
    Parameters::P_mag = 0.01;
    Parameters::T_mag = 0.0;
  }
  else if (Parameters::Problem_case == Parameters::Axisymmetric_shear_buckling)
  {
    nstep = 100;
    dp_mag = 0.0;
    dt_mag = 0.000001;
    Parameters::P_mag = 0.001;
    Parameters::T_mag = 0.0;
  }
  else if (Parameters::Problem_case ==
           Parameters::Nonaxisymmetric_shear_buckling)
  {
    nstep = 100;
    dp_mag = 0.0;
    dt_mag = 0.000001;
    Parameters::P_mag = 0.001;
    Parameters::T_mag = 0.0;
  }

  // Loop
  for (unsigned i = 0; i < nstep; i++)
  {
    oomph_info << "Solving for P = " << Parameters::P_mag
               << " ; Tau = " << Parameters::T_mag << "\n";

    // Do it
    problem.newton_solve();

    // Document it
    problem.doc_solution();

    // Bump it
    Parameters::T_mag += dt_mag;
    Parameters::P_mag += dp_mag;
  }

} // End of main
