/*
 * Adapted from:
 * This file is part of TinyAD and released under the MIT license.
 * Author: Patrick Schmidt
 */





// #include <igl/opengl/glfw/Viewer.h>
#include "polyscope/polyscope.h"

#include "VizHook.h"



#include <thread>

// /**
//  * Compute tutte embedding with boundary on circle.
//  * Per-vertex 2D coordinates returned as n_vertices-by-2 matrix.
//  */
// inline Eigen::MatrixXd tutte_embedding(
//     const Eigen::MatrixXd& _V,
//     const Eigen::MatrixXi& _F)
// {
//   Eigen::VectorXi b; // #constr boundary constraint indices
//   Eigen::MatrixXd bc; // #constr-by-2 2D boundary constraint positions
//   Eigen::MatrixXd P; // #V-by-2 2D vertex positions
//   igl::boundary_loop(_F, b); // Identify boundary vertices
//   igl::map_vertices_to_circle(_V, b, bc); // Set boundary vertex positions
//   igl::harmonic(_F, b, bc, 1, P); // Compute interior vertex positions

//   return P;
// }


static PhysicsHook *hook = NULL;


void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
    {
        hook->run();
    }
    else
        hook->pause();
}

void resetSimulation()
{
    
    if (!hook)
        return;
    std::cout << "try to reset" << std::endl;
    hook->reset();
}

void drawGUICallback()
{
	ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
							   // instead of full width. Must have 
							   // matching PopItemWidth() below.

    if(hook->showSimButtons() || true)
    {
        if (ImGui::CollapsingHeader("Start Simulation.", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Run/Pause Sim"))
            {
                toggleSimulation();
            }
            if (ImGui::Button("Reset Sim"))
            {
                resetSimulation();
            }
        }
    }
    hook->drawGUI();
	ImGui::PopItemWidth();
    hook->render();
}

int main(int argc, char **argv) {
  // Configure the argument parser

// Comment out parser.  Should be easy to add in commandline scripting with this though.

  // args::ArgumentParser parser("A simple demo of Polyscope with libIGL.\nBy "
  //                             "Nick Sharp (nsharp@cs.cmu.edu)",
  //                             "");
  // args::Positional<std::string> inFile(parser, "mesh", "input mesh");

  // // Parse args
  // try {
  //   parser.ParseCLI(argc, argv);
  // } catch (args::Help) {
  //   std::cout << parser;
  //   return 0;
  // } catch (args::ParseError e) {
  //   std::cerr << e.what() << std::endl;

  //   std::cerr << parser;
  //   return 1;
  // }


  // Options
  polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  polyscope::options::autocenterStructures = true;
// polyscope::options::autoscaleStructures = true;

  // Initialize polyscope
  polyscope::init();

  hook = new VizHook();
  hook->reset();

  polyscope::state::userCallback = drawGUICallback;



  polyscope::show();

  return 0;
}



/**
 * Injectively map a disk-topology triangle mesh to the plane
 * and optimize the symmetric Dirichlet energy via projected Newton.
 */
//   // View resulting parametrization
//   igl::opengl::glfw::Viewer viewer;
//   viewer.core().is_animating = true;
//   viewer.data().set_mesh(P, F);
//   viewer.core().camera_zoom = 2;
//   viewer.data().show_lines = false;
//   Eigen::MatrixXd N;
//   igl::per_vertex_normals(V,F,N);
//   viewer.data().set_colors( ((N.array()*0.5)+0.5).eval());
//   viewer.callback_pre_draw = [&] (igl::opengl::glfw::Viewer& viewer)
//   {
//     if(redraw)
//     {
//       viewer.data().set_vertices(P);
//       viewer.core().align_camera_center(P);
//       viewer.core().camera_zoom = 2;
//       {
//         std::lock_guard<std::mutex> lock(m);
//         redraw = false;
//       }
//     }
//     return false;
//   };
//   viewer.launch();