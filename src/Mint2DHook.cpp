    
    
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

#include <igl/readOBJ.h>

#include "date.h"

#include <igl/on_boundary.h>


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


void Mint2DHook::initSimulation()
{

    // igl::readOBJ(std::string(SOURCE_PATH) + "/../shared/" + cur_mesh_name + ".obj", V, F);
    igl::readOBJ("/home/josh/Documents/mint_redux/geometry-processing-starter-kit/tools/shared/" + cur_mesh_name + ".obj", V, F);

      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_subdiv.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/../shared/circle_1000.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole2.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_little_hole.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole_descimate.obj", V, F);

      // std::chrono::time_point now_clock = std::chrono::system_clock::now();
      // std::chrono::year_month_day now_time = std::chrono::floor<std::chrono::day>(now_clock);

// namespace C = std::chrono;
    // namespace D = date;
    // namespace S = std;

    auto tp = std::chrono::system_clock::now(); // tp is a C::system_clock::time_point
    auto dp = date::floor<date::days>(tp);  // dp is a sys_days, which is a
                                      // type alias for a C::time_point
    auto ymd = date::year_month_day{dp};
    auto time = date::make_time(std::chrono::duration_cast<std::chrono::milliseconds>(tp-dp));


      int month = static_cast<unsigned>( ymd.month() ); 
      int day = static_cast<unsigned>( ymd.day() );
      int hour = time.hours().count();
      int minute = time.minutes().count();
      cur_log_folder = "../../results/" + cur_mesh_name + "_" + std::to_string(month) + "_" + std::to_string(day) + "_" + std::to_string(hour); 
      int mk1succeeded = mkdir(cur_log_folder.c_str(), 0777);

      
      cur_log_folder = cur_log_folder + "/" + std::to_string(minute); // + std::to_string(now_time.month()) + "_" + std::to_string(now_time.day());
      std::cout << "log folder path: " << cur_log_folder << std::endl;
      int mk2succeeded = mkdir(cur_log_folder.c_str(), 0777);
      if (!mk2succeeded)
      {
        std::cout << "WARNING did not succeed in creating logging directory" << std::endl;
      }
// 

      cur_surf = Surface(V, F);


     
      // P = tutte_embedding(V, F); 

      // std::cout << "v size " <<  V.size() << " f size " << F.size() << " p size " << P.size() <<std::endl;

      buffer = 0;
      prev_energy = -1;


      cur_iter = 0; 
      inner_loop_iter = 0;

      // w_bound = 1e3; 
      // w_smooth = 1e-5; 
      // w_s_perp = 0;
      // w_curl = 1e5;

      w_bound = 1e6; 
      w_smooth = 10; // 1e4; // 1e3; 
      w_smooth_vector = 1;
      w_s_perp = 0; // 1e1 
      w_curl = 1e2; // 1e3;
      w_attenuate = 1; // 1e2;

      identity_weight = 1e-6;

      useProjHessian = true;

      // current_element = Field_View::vec_dirch;
 
      polyscope::removeAllStructures();
      // renderP.resize(P.rows(), 3);
      // renderP << P, Eigen::MatrixXd::Zero(P.rows(), 1);
      renderP = V;
      renderF = F; 

      polyscope::registerSurfaceMesh("cur state", renderP, renderF);
      // polyscope::getSurfaceMesh()->setEdgeWidth(.6);
      polyscope::getSurfaceMesh()->edgeWidth = .6;
      polyscope::view::resetCameraToHomeView();

      frames = Eigen::MatrixXd::Zero(F.rows(), 2);
      deltas = Eigen::MatrixXd::Zero(F.rows(), 4);
    //   metadata = Eigen::MatrixXd::Zero(F.rows(), 2);
      curls_sym = Eigen::VectorXd::Zero(F.rows());
      curls_primal = Eigen::VectorXd::Zero(F.rows());
      smoothness_primal = Eigen::VectorXd::Zero(F.rows());
      smoothness_sym = Eigen::VectorXd::Zero(F.rows());


      // frames = Eigen::MatrixXd::Random(F.rows(), 2);


      // bound_edges.resize(F.rows(),2);
      // const Eigen::MatrixXi blah = F;

      Eigen::MatrixXi bound_edges;

      Eigen::MatrixXi K;

      igl::on_boundary(F,bound_face_idx, K);
      // igl::boundary_facets(F, bound_edges, bound_face_idx, K);

      int nbf = bound_face_idx.size();
      for(int i = 0; i < nbf; i++)
      {

        if (bound_face_idx(i) == 1)
        {
          Eigen::VectorXd v0 = V.row(F(i,0));
          Eigen::VectorXd v1 = V.row(F(i,1));
          Eigen::VectorXd v2 = V.row(F(i,2));
          Eigen::VectorXd c = (( v0 + v1 + v2 ) / 3);

          if(std::sqrt(c.squaredNorm()) < .45)
          {
            bound_face_idx(i) = -1;
          }
          else
          {
            c.normalize();
            frames.row(i) = Eigen::Vector2d(c(1),-c(0)); // circulation 
            // frames.row(i) = Eigen::Vector2d(c(0),c(1)); // diverging

          }

          
          // std::cout << "i" << i << "bound_face_idx(i)" << bound_face_idx(i) << std::endl;
        }

      }
      }











    bool Mint2DHook::simulateOneStep()
    {
        // auto func = func_wrapper->get();
        if (cur_iter < max_iters)
        {
            cur_iter++;
            inner_loop_iter++;

 auto t1 = std::chrono::high_resolution_clock::now();

             TINYAD_DEBUG_OUT("Energy in iteration " << cur_iter << ": " << f);

// func.eval_with_hessian_proj

            // if( dec < convergence_eps )
            // {
            //   buffer -= 1;
            //   identity_weight = identity_weight * 10.;
            // }

            if ( dec < convergence_eps || (inner_loop_iter > 300 && dec / f < 1e-5))
            {
              std::cout << "***** current decrement: " << dec << std::endl;
              buffer = 5;
              identity_weight = 1e-6;
              if (w_smooth_vector > 0)
              {
                w_smooth_vector = 0;
                
              }
              else {
                w_attenuate = w_attenuate / 10.;
                std::cout << "New attenuation value is set to: " << w_attenuate << std::endl;
                // inner_loop_iter = 0;
                if (w_attenuate < 1e-12)
                  cur_iter = max_iters;
              }

              useProjHessian = true;
                 

                // Eigen::MatrixXd tmp =TinyAD::to_passive(H_proj);
                //  igl::writeDMAT("converged_hessian.dmat",tmp,true);
            }

            // Eigen::MatrixXd tmp =TinyAD::to_passive(H_proj);
            // igl::writeDMAT("curr_hessian.dmat",tmp,true);

            // cur_iter = max_iters; // break
            x = TinyAD::line_search(x, d, f, g, func, 1., .8, 512, 1e-3);

            auto t2 = std::chrono::high_resolution_clock::now();
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            std::cout << ms_int.count() << "ms\n";



            ///// Move this out 
            func.x_to_data(x, [&] (int f_idx, const Eigen::VectorXd& v) {
                frames.row(f_idx) = v.head<2>();
                // metadata.row(f_idx) = v.segment(2, 2);
                deltas.row(f_idx) = v.tail<4>();
                // if (bound_face_idx(f_idx) == 1)
                // {
                //   frames.row(f_idx) = frames_orig.row(f_idx);
                // }
                });

    


        }
        else if (cur_iter == max_iters) 
        {
            TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));
            cur_iter++;

             // FINAL LOGGING.  
        }
        else{
            this->pause();
        }


        
        return false;
    }
