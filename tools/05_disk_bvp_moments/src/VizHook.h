#include "PhysicsHook.h"
#include "Mint2DHook.h"

#include "Surface.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


#include <Eigen/IterativeLinearSolvers>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include <igl/readOBJ.h>

#include <igl/on_boundary.h>

#include <igl/writeDMAT.h>

#include <sys/stat.h>

#include <igl/map_vertices_to_circle.h>

#include "date.h"


#include <chrono>

// #include <fstream>
#include <sys/stat.h>

#include "UtilsMisc.h"



class VizHook : public Mint2DHook
{
public:
    VizHook() : Mint2DHook() {
      current_element = Field_View::vec_norms;
    }

//     virtual void drawGUI()
//     {
// 	//	ImGui::SliderFloat("k scale", &k_scale, 0.0001f, 2.0f, "k * %.3f");
// 	//	ImGui::SliderFloat("dt scale", &dt_scale, 0.0001f, 2.0f, "dt * %.3f");
// 		// ImGui::InputFloat("k scale", &k_scale);
// 		// ImGui::InputFloat("dt scale", &dt_scale);

//     ImGui::InputDouble("Smoothness Weight", &w_smooth);
//     ImGui::InputDouble("S Perp Weight", &w_s_perp);
//     ImGui::InputDouble("Curl Weight", &w_curl);
//     ImGui::InputDouble("Bound Weight", &w_bound);


// /// Whatever maybe make this a dropdown eventually 
// // From line 556 of imgui demo: https://skia.googlesource.com/external/github.com/ocornut/imgui/+/refs/tags/v1.73/imgui_demo.cpp
//             const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Vector Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
//             const char* current_element_name = (current_element >= 0 && current_element < Field_View::Element_COUNT) ? element_names[current_element] : "Unknown";
//             ImGui::PushItemWidth(300);
//             ImGui::SliderInt("Shading Mode", (int *) &current_element, 0, Element_COUNT - 1, current_element_name);
//             ImGui::PopItemWidth();
          
          
          
//             // ImGui::SameLine(); // ImGui::HelpMarker("Using the format string parameter to display a name instead of the underlying integer.");

//     }


    virtual void drawGUI()
    {
      Mint2DHook::drawGUI();

    }




    virtual void initSimulation()
    {

      // cur_mesh_name = "circle_subdiv";

      // cur_mesh_name = "circle";
      cur_mesh_name = "circle_1000";

      igl::readOBJ(std::string(SOURCE_PATH) + "/../shared/" + cur_mesh_name + ".obj", V, F);

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
      metadata = Eigen::MatrixXd::Zero(F.rows(), 2);
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

      int nedges = cur_surf.nEdges();

      rots.clear();// 
      rstars.clear();
      e_projs.clear();
      e_projs_primal.clear();

      e_projs2.resize(nedges,4); 

 

      for (int i = 0; i < nedges; i++)
      {

        Eigen::Vector3d estart = V.row(cur_surf.data().edgeVerts(i,0));
        Eigen::Vector3d eend = V.row(cur_surf.data().edgeVerts(i,1));
        Eigen::Vector3d edge_dir = (eend - estart).normalized();
        Eigen::Matrix2d e_to_x;
        e_to_x << edge_dir(0),edge_dir(1),-edge_dir(1),edge_dir(0); // Note this rotates the edge into [1,0]
        // std::cout << e_to_x * edge_dir.head(2) << std::endl<< std::endl; // sanity check.
// std::cout << flatten(edge_dir.head(2) * edge_dir.head(2).transpose()) << std::endl<< std::endl;

        rots.push_back(e_to_x);

        Eigen::Vector4d e_proj = rstar_xcomp_from_r(e_to_x);
        e_projs.push_back(e_proj);
        e_projs_primal.push_back(edge_dir.head(2));
        e_projs2.row(i) = e_proj;

      }
// std::cout << e_projs2 << std::endl;


    std::cout << e_projs2.rows() << std::endl;
    std::cout << "blah" << std::endl;

      frames_orig = frames;

      // std::cout << frames << std::endl;



      vc.d().frames.resize(frames.rows(), 3);
      vc.d().frames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);

      polyscope::getSurfaceMesh()->addFaceVectorQuantity("orig normals", vc.d().frames); //   ( ((N.array()*0.5)+0.5).eval());
      // polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", frames.rowwise().squaredNorm())->setEnabled(true); //   ( ((N.array()*0.5)+0.5).eval());



      // Set up function with 2D vertex positions as variables.
      func = TinyAD::scalar_function<6>(TinyAD::range(F.rows()));

      // Add objective term per face. Each connecting 3 vertices.
      func.add_elements<4>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
          {
          // Evaluate element using either double or TinyAD::Double
          using T = TINYAD_SCALAR_TYPE(element);



          // Get variable 2D vertex positions
          Eigen::Index f_idx = element.handle;
          Eigen::VectorX<T> s_curr = element.variables(f_idx);
          Eigen::Vector2<T> curr =  s_curr.head(2);
          Eigen::Vector4<T> delta = s_curr.tail(4);

          Eigen::Matrix2<T> currcurr = curr*curr.transpose();
          Eigen::Vector4<T> currcurrt = flatten(currcurr);
          // metadata field
          // Eigen::Vector2<T> metadata = s_curr.segment(2, 2);


          if (bound_face_idx(f_idx) == 1)
          {

            Eigen::Vector2<T> targ = frames_orig.row(f_idx);
            return w_bound*(curr-targ).squaredNorm() + w_bound*delta.squaredNorm();
          }

          if (bound_face_idx(f_idx) == -1)
          {
            T ret = w_bound*delta.squaredNorm();
            for(int i = 0; i < 3; i++)
            {
              int neighbor_edge_idx = cur_surf.data().faceNeighbors(f_idx, i);
              if(neighbor_edge_idx > -1)
              {
                Eigen::VectorX<T> s_n = element.variables(neighbor_edge_idx);
                Eigen::Vector2<T> n_i = s_n.head(2);
                // ret = ret + (n_i-curr).squaredNorm() * w_smooth;
                Eigen::Matrix2<T> nini = n_i*n_i.transpose();
                Eigen::Vector4<T> ninit = flatten(nini);
                ret = ret + (ninit-currcurrt).squaredNorm() * w_smooth * w_attenuate;
                // ret = ret + (n_i*n_i.transpose()-currcurr).norm() * w_smooth;
              }
            }
            return ret;
          }



          Eigen::VectorX<T> s_a = element.variables(cur_surf.data().faceNeighbors(f_idx, 0));
          Eigen::VectorX<T> s_b = element.variables(cur_surf.data().faceNeighbors(f_idx, 1));
          Eigen::VectorX<T> s_c = element.variables(cur_surf.data().faceNeighbors(f_idx, 2));



          Eigen::Vector2<T> a = s_a.head(2);
          Eigen::Vector2<T> b = s_b.head(2);
          Eigen::Vector2<T> c = s_c.head(2);

          Eigen::Matrix2<T> aa = a*a.transpose();
          Eigen::Matrix2<T> bb = b*b.transpose();
          Eigen::Matrix2<T> cc = c*c.transpose();


          Eigen::Vector4<T> a_delta = s_a.tail(4);
          Eigen::Vector4<T> b_delta = s_b.tail(4);
          Eigen::Vector4<T> c_delta = s_c.tail(4);

          Eigen::Vector4<T> aat = flatten(aa);
          Eigen::Vector4<T> bbt = flatten(bb);
          Eigen::Vector4<T> cct = flatten(cc);

          aat = aat + a_delta;
          bbt = bbt + b_delta; 
          cct = cct + c_delta;



                    // return (T) 0.;
         






          Eigen::Vector2<T> curr_normalized = curr.normalized();
          Eigen::Vector2<T> curr_perp; // = curr_normalized;
          curr_perp(0) = curr_normalized(1);
          curr_perp(1) = -curr_normalized(0);

          // T s_perp_term = pow(a.dot(curr_perp),2) + pow(b.dot(curr_perp),2) + pow(c.dot(curr_perp), 2);

          // T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (curr_perp * curr_perp.transpose())).norm();
          T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (currcurrt)).squaredNorm();

          // T primal_dirichlet_term = (a + b + c - 3*curr).squaredNorm();
          // T dirichlet_term = (aat+bbt+cct-3*currcurrt).squaredNorm();

          T primal_dirichlet_term = (a - curr).squaredNorm() + (b - curr).squaredNorm() + (c - curr).squaredNorm();
          T dirichlet_term = (aat-currcurrt).squaredNorm() + (bbt-currcurrt).squaredNorm() + (cct-currcurrt).squaredNorm();



          // T dirichlet_term = (aa + bb + cc - 3*currcurr).norm();
  // dirichlet_term += 1e-5*abs(dirichlet_term - metadata(0));
          // T delta_rescale = std::max(frames.row(f_idx).squaredNorm(), 1e-8);
          // delta_rescale = (.0001 + 1./delta_rescale);
          T delta_rescale = 1.;
          // std::cout << delta_rescale << std::endl;

          smoothness_primal(f_idx) = TinyAD::to_passive(primal_dirichlet_term);
          smoothness_sym(f_idx) = TinyAD::to_passive(dirichlet_term);

          // T delta_dirichlet = (a_delta+b_delta+c_delta-3*delta).squaredNorm()*delta_rescale;

          T delta_norm_term = delta_rescale * delta.squaredNorm();// + delta_dirichlet;

          Eigen::Vector4d ea = e_projs2.row(cur_surf.data().faceEdges(f_idx, 0));
          Eigen::Vector4d eb = e_projs2.row(cur_surf.data().faceEdges(f_idx, 1));
          Eigen::Vector4d ec = e_projs2.row(cur_surf.data().faceEdges(f_idx, 2));

          // Eigen::Vector4<T> ea = e_projs.at(cur_surf.data().faceEdges(f_idx, 0));
          // Eigen::Vector4<T> eb = e_projs.at(cur_surf.data().faceEdges(f_idx, 1));
          // Eigen::Vector4<T> ec = e_projs.at(cur_surf.data().faceEdges(f_idx, 2));

          T curl_term = pow(ea.dot(aat + a_delta) - ea.dot(currcurrt + delta),2);
          curl_term +=  pow(eb.dot(bbt + b_delta) - eb.dot(currcurrt + delta),2);
          curl_term +=  pow(ec.dot(cct + c_delta) - ec.dot(currcurrt + delta),2);

          curls_sym(f_idx) = TinyAD::to_passive(curl_term);


          Eigen::Vector2<T> ea_primal = e_projs_primal.at(cur_surf.data().faceEdges(f_idx, 0));
          Eigen::Vector2<T> eb_primal = e_projs_primal.at(cur_surf.data().faceEdges(f_idx, 1));
          Eigen::Vector2<T> ec_primal = e_projs_primal.at(cur_surf.data().faceEdges(f_idx, 2));

          T curl_term_primal = pow(ea_primal.dot(a) - ea_primal.dot(curr),2);
          curl_term_primal +=  pow(eb_primal.dot(b) - eb_primal.dot(curr),2);
          curl_term_primal +=  pow(ec_primal.dot(c) - ec_primal.dot(curr),2);

          curls_primal(f_idx) = TinyAD::to_passive(curl_term_primal);

          // curl_term += 1e-5*abs(curl_term - metadata(1));


          // T atten = 1./(cur_iter*cur_iter + 1);
          // T atten = 1./(cur_iter + 1);

          // T atten = 1.;

          T delta_weight = .1; // std::min(w_curl/100., 1./w_attenuate);
          T w_curl_new = std::min(1e4, 1./w_attenuate) * w_curl;

          T ret = delta_norm_term * delta_weight;
          if (w_smooth_vector > 0)
            return w_smooth_vector * primal_dirichlet_term + ret;
          if (w_smooth > 0)
            ret = ret + w_attenuate * w_smooth * dirichlet_term;
          if (w_s_perp > 0)
            ret = ret + w_attenuate * w_s_perp * s_perp_term;
          if (w_curl_new > 0)
            ret = ret + w_curl_new * curl_term;

          return ret;

// (w_smooth * dirichlet_term + 
//                   w_s_perp * s_perp_term) *  + 
//                   w_curl*curl_term  + 
//                  delta_norm_term * delta_weight;

          // return (w_smooth * dirichlet_term + 
          //         w_s_perp * s_perp_term) * atten + 
          //        w_curl*curl_term  + 
          //        delta_norm_term;

                 
          });

      // Assemble inital x vector from P matrix.
      // x_from_data(...) takes a lambda function that maps
      // each variable handle (vertex index) to its initial 2D value (Eigen::Vector2d).
        x = func.x_from_data([&] (int f_idx) {
          Eigen::VectorXd ret;
          ret = Eigen::VectorXd::Zero(6); // resize(10);
          ret.head(2) = Eigen::VectorXd::Random(2) * 0.0001;
          // ret << frames.row(f_idx), deltas.row(f_idx);
          return ret;
          });

    }


    virtual void updateRenderGeometry()
    {
      Mint2DHook::updateRenderGeometry();

    }


    virtual bool simulateOneStep()
    {
        if (cur_iter < max_iters)
        {
            cur_iter++;
            inner_loop_iter++;


            auto t1 = std::chrono::high_resolution_clock::now();


            // auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
            auto [f, g, H_proj] = func.eval_with_derivatives(x);
            TINYAD_DEBUG_OUT("Energy in iteration " << cur_iter << ": " << f);
            // std::cout<<"the number of nonzeros "<<H_proj.nonZeros() << "number of non-zeros per dof " << H_proj.nonZeros() / (6*F.rows()) << " # rows " << H_proj.rows() << " faces " << F.rows() <<std::endl;

            // std::cout<<"the number of nonzeros "<<H_proj.nonZeros()<<std::endl;
            Eigen::VectorXd d;
            double dec;
            // d = TinyAD::newton_direction(g, H_proj, solver);
             // = TinyAD::newton_decrement(d, g);

            if (prev_energy < 0)
            {
              prev_energy = f + 100 * convergence_eps;
            }
            
            try
            {
              if (w_smooth_vector > 0 || useProjHessian)
              {
                auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
                f = f_h;
                g = g_h;
                H_proj = H_proj_h;
                d = TinyAD::newton_direction(g, H_proj, solver, 0.);
                dec = TinyAD::newton_decrement(d, g);

                if ( dec / f < 1e-3)
                {
                  useProjHessian = false;
                  std::cout << "switch off projected hessian to fine-tune result" << std::endl;
                }


              }
              else
              {
                d = TinyAD::newton_direction(g, H_proj, solver, identity_weight);
                dec = TinyAD::newton_decrement(d, g);
                identity_weight = identity_weight / 2.;
              }
              
            }
            catch(const std::exception& e)
            {
              auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
              f = f_h;
              g = g_h;
              H_proj = H_proj_h;
              d = TinyAD::newton_direction(g, H_proj, solver);
              dec = TinyAD::newton_decrement(d, g);
              if ( !useProjHessian )
                identity_weight = identity_weight * 10.;
            }
            
            // 
            std::cout << "current decrement: " << dec << std::endl;
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


    virtual void renderRenderGeometry()
    {
      Mint2DHook::renderRenderGeometry();
    }



protected:
  // Read mesh and compute Tutte embedding



  Eigen::MatrixXd metadata;


  Surface cur_surf;


      std::vector<Eigen::Matrix2d> rots;// 
      std::vector<Eigen::Matrix4d> rstars;
      std::vector<Eigen::Vector4d> e_projs;
      std::vector<Eigen::Vector2d> e_projs_primal;

      Eigen::MatrixXd e_projs2;



  
  Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 

  // Eigen::MatrixXd renderFrames;
  // Eigen::MatrixXd renderDeltas;





  bool useProjHessian = true;


  
  std::vector<Eigen::Matrix2d> rest_shapes;
// %%% 2 + 4 + 4
  decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;
  Eigen::VectorXd x;


  double prev_energy; 


  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;


    
};