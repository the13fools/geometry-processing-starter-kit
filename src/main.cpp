#include "polyscope/polyscope.h"

#include <thread>

#include "VizHook.h"


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


// void callback() {

//   static int numPoints = 2000;
//   static float param = 3.14;

//   ImGui::PushItemWidth(100);

//   // Curvature
//   if (ImGui::Button("add curvature")) {
//  //   addCurvatureScalar();
//   }
  
//   // Normals 
//   if (ImGui::Button("add normals")) {
// //    computeNormals();
//   }

//   // Param
//   if (ImGui::Button("add parameterization")) {
//  //   computeParameterization();
//   }

//   // Geodesics
//   if (ImGui::Button("compute distance")) {
//  //   computeDistanceFrom();
//   }
//   ImGui::SameLine();
// //  ImGui::InputInt("source vertex", &iVertexSource);

//   ImGui::PopItemWidth();
// }

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

  // Initialize polyscope
  polyscope::init();

  hook = new VizHook();
  hook->reset();

  polyscope::state::userCallback = drawGUICallback;

  polyscope::show();

  return 0;
}

