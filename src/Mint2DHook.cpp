    
    
// #include <mutex>
// #include <thread>

// #include "PhysicsHook.h"
#include "Mint2DHook.h"


#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

// #include <fstream>
#include <sys/stat.h>
#include <iostream>


// #include "polyscope/polyscope.h"

// #include "polyscope/surface_mesh.h"

// #include <Eigen/Core>


// #include "VizHelper.h"

// #include "UtilsMisc.h"

// // #include <TinyAD/ScalarFunction.hh>
// // #include <TinyAD/Utils/NewtonDirection.hh>
// // #include <TinyAD/Utils/NewtonDecrement.hh>
// // #include <TinyAD/Utils/LineSearch.hh>

// // #include <fstream>
// #include <sys/stat.h>
    
    
    void Mint2DHook::drawGUI()
    {
        //	ImGui::SliderFloat("k scale", &k_scale, 0.0001f, 2.0f, "k * %.3f");
        //	ImGui::SliderFloat("dt scale", &dt_scale, 0.0001f, 2.0f, "dt * %.3f");
            // ImGui::InputFloat("k scale", &k_scale);
            // ImGui::InputFloat("dt scale", &dt_scale);

        ImGui::InputDouble("Smoothness Weight", &w_smooth);
        ImGui::InputDouble("S Perp Weight", &w_s_perp);
        ImGui::InputDouble("Curl Weight", &w_curl);
        ImGui::InputDouble("Bound Weight", &w_bound);


/// Whatever maybe make this a dropdown eventually 
// From line 556 of imgui demo: https://skia.googlesource.com/external/github.com/ocornut/imgui/+/refs/tags/v1.73/imgui_demo.cpp
            const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Vector Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
            const char* current_element_name = (current_element >= 0 && current_element < Field_View::Element_COUNT) ? element_names[current_element] : "Unknown";
            ImGui::PushItemWidth(300);
            ImGui::SliderInt("Shading Mode", (int *) &current_element, 0, Element_COUNT - 1, current_element_name);
            ImGui::PopItemWidth();
          
          
          
            // ImGui::SameLine(); // ImGui::HelpMarker("Using the format string parameter to display a name instead of the underlying integer.");

    }

    void Mint2DHook::updateRenderGeometry()
    {



      vc.d().frames.resize(frames.rows(), 3);
      vc.d().frames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);
      vc.d().deltas = deltas;

      vc.updateVizState();

      vc.d().vec_curl = curls_primal;
      vc.d().sym_curl = curls_sym;

      vc.d().frame_smoothness = smoothness_primal;
      vc.d().moment_smoothness = smoothness_sym;



    }




    // Need to fill out viewer for each of: Field_View { vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual,
    void Mint2DHook::renderRenderGeometry()
    {
		polyscope::getSurfaceMesh()->updateVertexPositions(renderP);
        
        polyscope::getSurfaceMesh()->centerBoundingBox();
        polyscope::getSurfaceMesh()->resetTransform();

        switch (current_element) 
        { 
            case Field_View::vec_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", vc.d().frame_norms)->setEnabled(true);
                } 
                break;

            case Field_View::delta_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("delta_norms", vc.d().delta_norms)->setEnabled(true);
                } 
                break; 
              
            case Field_View::primal_curl_residual: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("primal_curl_residual", vc.d().vec_curl)->setEnabled(true);
                } 
                break;             
            case Field_View::sym_curl_residual: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("sym_curl_residual", vc.d().sym_curl)->setEnabled(true);
                } 
                break;           
            case Field_View::vec_dirch: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_dirch", vc.d().frame_smoothness)->setEnabled(true);
                } 
                break; 
            case Field_View::moment_dirch: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("moment_dirch", vc.d().moment_smoothness)->setEnabled(true);
                } 
                break; 

            default: 
                { 
                    // std::cout << "Unknown color!"; 
                } 
                break; 
        } 

        auto vectors = polyscope::getSurfaceMesh()->addFaceVectorQuantity("frames", vc.getFrames()); //   ( ((N.array()*0.5)+0.5).eval());
        vectors->vectorColor = glm::vec3(.7,.7,.7);
        
        // vectors->setVectorLengthScale(1., true);
        // vectors->setEnabled(true);
        // vectors->setVectorColor(glm::vec3(.7,.7,.7));
        
        polyscope::requestRedraw();   

        // std::ifstream file;
        std::string cur_log_file =  cur_log_folder + "/cur_file_iter_" + std::to_string(10000 + cur_iter) + ".png";

        struct stat buffer;   
        bool exists = (stat(cur_log_file.c_str(), &buffer) == 0); 

        

        if (!exists)
        {
            polyscope::screenshot(cur_log_file, true);
            std::cout << cur_log_file << std::endl;
            std::cout << "Current File Path to log" << std::endl;
        }


        // // opening the file
        // file.open(cur_log_file.c_str(), 'r');

        // if (file == NULL)
              




    }