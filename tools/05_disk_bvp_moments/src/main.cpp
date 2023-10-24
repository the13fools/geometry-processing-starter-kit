/*
 * Adapted from:
 * This file is part of TinyAD and released under the MIT license.
 * Author: Patrick Schmidt
 */



#include "polyscope/polyscope.h"

#include "VizHook.h"

#include <thread>



static VizHook *hook = NULL;


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

  double w_bound_prev = hook->w_bound;
  double w_smooth_prev = hook->w_smooth;
  double w_curl_prev = hook->w_curl;
  double w_s_perp_prev = hook->w_s_perp;

    std::cout << "try to reset" << std::endl;
    hook->reset();

  hook->w_bound = w_bound_prev;
  hook->w_smooth = w_smooth_prev;
  hook->w_curl = w_curl_prev;
  hook->w_s_perp = w_s_perp_prev;

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

//   polyscope::GroundPlaneMode::None;
// polyscope::options::groundPlaneEnabled = false;

// polyscope::screenshotExtension = ".jpg";


  polyscope::show();

  return 0;
}