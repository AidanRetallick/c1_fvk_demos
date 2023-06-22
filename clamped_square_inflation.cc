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


namespace Parameters
{
 // Plate length
 double L = 1.0;
  // The plate thickness
 double Thickness = 0.01;
 // Poisson ratio
 double Nu = 0.5;
 // Membrane coupling coefficient
 double Eta = 12.0*(1.0-Nu*Nu)*Thickness*Thickness;
 
 // Control parameters
 double P_mag = 0.0;
 double T_mag = 0.0;
 
 // Mesh parameters
 double Element_area=0.1;
 
 // Output directory
 string Output_dir="RESLT";
 
 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& x, double& pressure)
 {
  // Parabolic forcing growing outwards
  pressure = P_mag * (x[0]*x[0]+x[1]*x[1]);
 }
 
 // Pressure wrapper so we can output the pressure function
 void get_pressure(const Vector<double>& x, Vector<double>& pressure)
 {
  pressure.resize(1);
  get_pressure(x,pressure[0]);
 }
 
 // Assigns the value of in-plane forcing depending on the position (x,y)
 void get_in_plane_force(const Vector<double>& x, Vector<double>& grad)
 {
  // Parabolic forcing decaying outwards
  grad[0]= T_mag * (1.0-x[0]*x[0])*(1.0-x[1]*x[1]);
  grad[1]= T_mag * (1.0-x[0]*x[0])*(1.0-x[1]*x[1]);
 }

 //--------------------------------------------------------------------- 
 // Functions to assign boundary conditions
 
 // Null function for any zero (homogenous) BCs
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
  delete (Surface_mesh_pt);
  delete (Bulk_mesh_pt);
  delete Boundary_pt;
  delete Boundary0_pt;
  delete Boundary1_pt;
  delete Boundary2_pt;
  delete Boundary3_pt;
 };

 /// Setup and build the mesh
 void build_mesh();

 /// Helper function to (re-)set boundary condition
 /// and complete the build of all elements
 void complete_problem_setup();

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// Print information about the parameters we are trying to solve for.
 void actions_before_newton_solve()
 {
  // For what control parameter are we about to solve?
  oomph_info << "-------------------------------------------------------" << std::endl;
  oomph_info << "Solving for p = " << Parameters::P_mag << std::endl;
  oomph_info << "     Doc_info = " << Doc_info.number() << std::endl;
  oomph_info << "-------------------------------------------------------" << std::endl;
 }

 /// Update after solve (empty)
 void actions_after_newton_solve() {}

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

 // Triangle Mesh Parameter Data
 // This is the data used to set-up the mesh, we need to store the pointers
 // HERE otherwise we will not be able to clean up the memory once we have
 // finished the problem.
 // The initial (and maximum) element area
 double Element_area;
 // The mesh parameters
 TriangleMeshParameters* Triangle_mesh_parameters_pt;
 TriangleMeshClosedCurve* Boundary_pt;
 TriangleMeshPolyLine* Boundary0_pt;
 TriangleMeshPolyLine* Boundary1_pt;
 TriangleMeshPolyLine* Boundary2_pt;
 TriangleMeshPolyLine* Boundary3_pt;

 /// Trace file to document norm of solution
 ofstream Trace_file;

 /// Delete traction elements and wipe the surface mesh
 void delete_traction_elements(Mesh* const &surface_mesh_pt);

 /// Pointer to "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to "surface" mesh
 Mesh* Surface_mesh_pt;

}; // end_of_problem_class

//==start_of_problem_constructor===========================================
/// Constructs problem, by calling the build_mesh() and then the
/// complete_problem_setup() helpers. Finally assign_eqn_numbers() sets up
/// the local-global dof/equation numbering scheme.
//=========================================================================
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem()
 :
 Element_area(Parameters::Element_area)
{
 // Build the mesh
 build_mesh();

 // Store number of bulk elements
 complete_problem_setup();

 char filename[100];
 ofstream Param_file;
 strcpy(filename, (Parameters::Output_dir + "/parameters.dat").c_str());
 Param_file.open(filename);

 // Output plate parameters
 Param_file << "L            " << Parameters::L         << std::endl
	    << "thickness    " << Parameters::Thickness << std::endl
	    << "nu           " << Parameters::Nu        << std::endl
	    << "eta          " << Parameters::Eta       << std::endl
	    << "Element area " << Parameters::Element_area << std::endl;
 
 strcpy(filename, (Parameters::Output_dir + "/trace.dat").c_str());
 Trace_file.open(filename);
 
 oomph_info << "Number of equations: "
	    << assign_eqn_numbers() << '\n';
 
} // end Constructor


//==start_of_build_mesh====================================================
/// Build the rectangular mesh by specifying the location of the vertices
/// before assigning them to be the endpoints of the boundary edges.
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
 
 double L = Parameters::L;

 // Declare the vertices..
 Vector<double>
  V0(2,0.0),
  V1(2,0.0),
  V2(2,0.0),
  V3(2,0.0);
 // ..and place them at the corners of the rectangle
 V0[0] = -0.5*L;
 V0[1] = -0.5;
 V1[0] =  0.5*L;
 V1[1] = -0.5;
 V2[0] =  0.5*L;
 V2[1] =  0.5;
 V3[0] = -0.5*L;
 V3[1] =  0.5;
  
 // Declare the edges..
 Vector<Vector<double>>
  E0(2,Vector<double>(2,0.0)),
  E1(2,Vector<double>(2,0.0)),
  E2(2,Vector<double>(2,0.0)),
  E3(2,Vector<double>(2,0.0));
 // ..and assign them endpoints
 E0[0] = V0;
 E0[1] = V1;
 E1[0] = V1;
 E1[1] = V2;
 E2[0] = V2;
 E2[1] = V3;
 E3[0] = V3;
 E3[1] = V0;
  
 // Define boundaries from edges
 Boundary0_pt = new TriangleMeshPolyLine(E0, 0);
 Boundary1_pt = new TriangleMeshPolyLine(E1, 1);
 Boundary2_pt = new TriangleMeshPolyLine(E2, 2);
 Boundary3_pt = new TriangleMeshPolyLine(E3, 3);
 
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
 boundary_polyline_pt[0] = Boundary0_pt;
 boundary_polyline_pt[1] = Boundary1_pt;
 boundary_polyline_pt[2] = Boundary2_pt;
 boundary_polyline_pt[3] = Boundary3_pt;
 
 Boundary_pt = new TriangleMeshClosedCurve(boundary_polyline_pt);
 
 TriangleMeshParameters Triangle_mesh_parameters(Boundary_pt);
 
 // Set the maximum element area
 Triangle_mesh_parameters.element_area()=Element_area;
 
 // Build an assign bulk mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(Triangle_mesh_parameters);
 
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



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of
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
 // The free boundary condition is completely unpinned
 static const Vector<unsigned> free{};

 // In-plane dofs:
 // |  0  |  1  |
 // |  ux |  uy |
 // Possible boundary conditions on the in-plane displacement
 static const Vector<unsigned> pin_ux_dofs{0};
 static const Vector<unsigned> pin_uy_dofs{1};
 static const Vector<unsigned> pin_inplane_dofs{0,1};
 
 // Out-of-plane dofs:
 // |  0  |  1  |  2  |  3  |  4  |  5  |
 // |  w  | w_x | w_y | w_xx| w_xy| w_yy|
 // Possible boundary conditions depending on whether an edge is x-normal or
 // y-normal
 static const Vector<unsigned> resting_pin_xn_dofs{0,2,5};
 static const Vector<unsigned> resting_pin_yn_dofs{0,1,3};
 static const Vector<unsigned> sliding_clamp_xn_dofs{1,4};
 static const Vector<unsigned> sliding_clamp_yn_dofs{2,4};
 static const Vector<unsigned> true_clamp_xn_dofs{0,1,2,4,5};
 static const Vector<unsigned> true_clamp_yn_dofs{0,1,2,3,4};

 // Vectors to store which boundary conditions we are applying to each edge.
 Vector<Vector<unsigned>> pinned_u_dofs(4);
 Vector<Vector<unsigned>> pinned_w_dofs(4);
 
 // Choose BCs for each boundary
 pinned_u_dofs[0] = pin_inplane_dofs;
 pinned_u_dofs[1] = pin_inplane_dofs;
 pinned_u_dofs[2] = pin_inplane_dofs;
 pinned_u_dofs[3] = pin_inplane_dofs;
 pinned_w_dofs[0] = resting_pin_yn_dofs;
 pinned_w_dofs[1] = resting_pin_xn_dofs;
 pinned_w_dofs[2] = resting_pin_yn_dofs;
 pinned_w_dofs[3] = resting_pin_xn_dofs;
 
 // Loop over all the boundaries in our bulk mesh
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   // Number of elements on b
   const unsigned nb_element = Bulk_mesh_pt->nboundary_element(b);
   // Number of dofs we are pinning on b
   const unsigned n_pinned_u_dofs = pinned_u_dofs[b].size();
   const unsigned n_pinned_w_dofs = pinned_w_dofs[b].size();

   // Loop over the elements on boundary b
   for(unsigned e=0; e<nb_element; e++)
    {
     // Get pointer to bulk element adjacent to b
     ELEMENT* el_pt =
      dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(b,e));

     // Pin in-plane dofs at edge b - always null because we are applying
     // homogenous bcs
     for(unsigned i=0; i<n_pinned_u_dofs; i++)
      {
       unsigned idof=pinned_u_dofs[b][i];
       el_pt->fix_in_plane_displacement_dof(idof, b, Parameters::get_null_fct);
      }
     // Pin out-of-plane dofs at edge b - always null because we are applying
     // homogenous bcs
     for(unsigned i=0; i<n_pinned_w_dofs; i++)
      {
       unsigned idof=pinned_w_dofs[b][i];
       el_pt->fix_out_of_plane_displacement_dof(idof, b, Parameters::get_null_fct);
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
 
 unsigned npts;
 // Number of plot points for coarse output
 npts = 2;
 sprintf(filename, "%s/coarse_soln_%i.dat",
	 Parameters::Output_dir.c_str(),
	 Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();

 // Number of plot points for fine output
 npts = 10;
 sprintf(filename, "%s/soln_%i.dat",
	 Parameters::Output_dir.c_str(),
	 Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();


 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \""
	   << comment << "\"\n";
 some_file.close();

 // Write the pressure, degree of swelling and
 // deflection to the trace file
 //-----------------------------------------------------
 // Get the centre deflection first
 Vector<double> origin(2,0.0), s(2,0.0);
 GeomObject* centre_element_pt;
 Vector<double> w_centre(1,0.0);
  
 unsigned n_element=Bulk_mesh_pt->nelement();
 for(unsigned i=0; i<n_element; i++)
  {
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i))
    ->locate_zeta(origin, centre_element_pt, s);
   if(centre_element_pt!=NULL)
    {
     w_centre = dynamic_cast<ELEMENT*>(centre_element_pt)
      ->interpolated_u_foeppl_von_karman(s);
     break;
    }
  }
 Trace_file << Doc_info.number()          << " "
	    << Parameters::P_mag                     << " "
	    << w_centre[0]                       << endl;

 // Increment the doc_info numbers
 Doc_info.number()++;

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
/// Driver code for demo of unstructured C1 Foeppl-von Karman
/// elements
//============================================================
int main(int argc, char **argv)
{
 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
 
 // Create the problem
 UnstructuredFvKProblem<FoepplVonKarmanC1CurvableBellElement<4>>
  problem;
 
 // Set solver paramters
 problem.newton_solver_tolerance()=1e-8;
 problem.max_residuals()=1e4;
 problem.max_newton_iterations()=10;
 
 // Increase the pressure under the membrane
 Parameters::P_mag+=10;
 // Solve the system
 problem.newton_solve();
 // Document the current solution
 problem.doc_solution();
 
 // Print success
 oomph_info<<"Exiting Normally\n";
} //End of main
