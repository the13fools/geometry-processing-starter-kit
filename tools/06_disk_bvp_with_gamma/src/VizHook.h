#include "PhysicsHook.h"
#include "Surface.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <Eigen/Core>

#include <Eigen/IterativeLinearSolvers>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include <igl/readOBJ.h>

#include <igl/on_boundary.h>


#include <igl/map_vertices_to_circle.h>

#include <UtilsMisc.h>

#include <VizHelper.h>

enum Field_View { vec_norms, delta_norms, gamma_norms, vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual, gui_free, Element_COUNT };


class VizHook : public PhysicsHook
{
public:
    VizHook() : PhysicsHook() {
      current_element = Field_View::vec_norms;
    }

    virtual void drawGUI()
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
            const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Gamma Norms", "Vector Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
            const char* current_element_name = (current_element >= 0 && current_element < Field_View::Element_COUNT) ? element_names[current_element] : "Unknown";
            ImGui::PushItemWidth(300);
            ImGui::SliderInt("Shading Mode", (int *) &current_element, 0, Element_COUNT - 1, current_element_name);
            ImGui::PopItemWidth();
          
          
          
            // ImGui::SameLine(); // ImGui::HelpMarker("Using the format string parameter to display a name instead of the underlying integer.");

    }


template <typename ScalarType, int Rows, int Cols>
Eigen::Matrix<ScalarType, Rows * Cols, 1> flatten(const Eigen::Matrix<ScalarType, Rows, Cols>& matrix) {
    Eigen::Matrix<ScalarType, Rows * Cols, 1> flattened;
    flattened << matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1);
    return flattened;
}

template <typename ScalarType>
Eigen::Matrix<ScalarType, 2,2> fold(const Eigen::Matrix<ScalarType, 4, 1>& matrix) {
    Eigen::Matrix<ScalarType, 2,2> folded;
    folded << matrix(0), matrix(1), matrix(2), matrix(3);
    return folded;
}


    Eigen::Matrix4d rstar_from_r(Eigen::Matrix2d r)
    {
      Eigen::Matrix4d ret;
      double a = r(0,0);
      double b = r(0,1);
      double c = r(1,0);
      double d = r(1,1);
      ret << a*a, a*b, a*b, b*b,
             a*c, b*c, a*d, b*d,
             a*c, b*c, a*d, b*d,
             c*c, c*d, c*d, d*d;

      // ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
      return ret;
    }


    void updateJacobians(const Eigen::MatrixXd& frames, std::vector<Eigen::MatrixXd>& jacobians)
    {
      jacobians.clear();
      int nrows = frames.rows();
      for(int i = 0; i < nrows; i++)
      {
        Eigen::VectorXd cv = frames.row(i);
        Eigen::MatrixXd jac(4,2);
        jac.col(0) << 2*cv(0), cv(1), cv(1), 0;
        jac.col(1) << 0, cv(0), cv(0), 2*cv(1);
        jacobians.push_back(jac);
      }
    }


    Eigen::Vector4d rstar_xcomp_from_r(Eigen::Matrix2d r)
    {
      Eigen::Vector4d ret;
      double a = r(0,0);
      double b = r(0,1);
      double c = r(1,0);
      double d = r(1,1);
      ret << a*a, a*b, a*b, b*b;

      // ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
      return ret;
    }



    virtual void initSimulation()
    {

      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle.obj", V, F);

      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_subdiv.obj", V, F);
      igl::readOBJ(std::string(SOURCE_PATH) + "/../shared/circle_1000.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole2.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_little_hole.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole_descimate.obj", V, F);
// 

      cur_surf = Surface(V, F);


     
      // P = tutte_embedding(V, F); 

      // std::cout << "v size " <<  V.size() << " f size " << F.size() << " p size " << P.size() <<std::endl;

      cur_iter = 0; 


      // w_bound = 1e3; 
      // w_smooth = 1e-5; 
      // w_s_perp = 0;
      // w_curl = 1e5;

      w_bound = 1e5; 
      w_smooth = 1; // 10000; // 1e4; // 1e3; 
      w_s_perp = 0; // 1e1 
      w_curl = 0; // 1e3;
      w_attenuate = 1.;

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
      gammas = Eigen::MatrixXd::Zero(F.rows(), 4);
      metadata = Eigen::MatrixXd::Zero(F.rows(), 2);
      curls = Eigen::VectorXd::Zero(F.rows());

      // frames = Eigen::MatrixXd::Random(F.rows(), 2);
      updateJacobians(frames, jacobians);


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
        e_projs2.row(i) = e_proj;

      }
// std::cout << e_projs2 << std::endl;


    std::cout << e_projs2.rows() << std::endl;
    std::cout << "blah" << std::endl;

      frames_orig = frames;

      // std::cout << frames << std::endl;



      renderFrames.resize(frames.rows(), 3);
      renderFrames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);

      polyscope::getSurfaceMesh()->addFaceVectorQuantity("orig normals", renderFrames); //   ( ((N.array()*0.5)+0.5).eval());
      polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", frames.rowwise().squaredNorm())->setEnabled(true); //   ( ((N.array()*0.5)+0.5).eval());



      // Set up function with 2D vertex positions as variables.
      func = TinyAD::scalar_function<8>(TinyAD::range(F.rows()));

      // Add objective term per face. Each connecting 3 vertices.
      func.add_elements<4>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
          {
          // Evaluate element using either double or TinyAD::Double
          using T = TINYAD_SCALAR_TYPE(element);



          // Get variable 2D vertex positions
          Eigen::Index f_idx = element.handle;
          Eigen::VectorX<T> s_curr = element.variables(f_idx);
          // Eigen::Vector2<T> curr =  s_curr.head(2);
          // Eigen::Vector2<T> curr = frames.row(f_idx);
          Eigen::Vector4<T> delta = s_curr.tail(4);
          // Eigen::Vector2<T> gamma = s_curr.segment(2, 2);

          Eigen::Vector4<T> gamma = s_curr.head(4);
          Eigen::Matrix2<T> gamma_folded = fold(gamma);

          Eigen::Vector2<T> primal = frames.row(f_idx);
          Eigen::Vector2<T> curr = primal;// + s_curr.segment(2, 2);
          
          Eigen::Matrix2<T> currcurr = curr*curr.transpose();
          Eigen::Vector4<T> currcurrt = flatten(currcurr) + gamma;
          
          
          // metadata field
          // Eigen::Vector2<T> metadata = s_curr.segment(2, 2);


          if (bound_face_idx(f_idx) == 1)
          {

            Eigen::Vector2<T> targ = frames_orig.row(f_idx);
            Eigen::Matrix2<T> targtarg = targ*targ.transpose();
            Eigen::Vector4<T> targtargt = flatten(targtarg);
            return w_bound*(currcurrt-targtargt).squaredNorm() + w_bound*delta.squaredNorm();
          }

          if (bound_face_idx(f_idx) == -1)
          {
            T ret = w_bound*(delta.squaredNorm());
            for(int i = 0; i < 3; i++)
            {
              int neighbor_edge_idx = cur_surf.data().faceNeighbors(f_idx, i);
              if(neighbor_edge_idx > -1)
              {
                // Eigen::VectorX<T> s_n = element.variables(neighbor_edge_idx);
                // Eigen::Vector2<T> n_i = s_n.head(2);
                // ret = ret + (n_i-curr).squaredNorm() * w_smooth;
                Eigen::Vector2<T> n_i = frames.row(neighbor_edge_idx);
                Eigen::Matrix2<T> nini = n_i*n_i.transpose();
                Eigen::Vector4<T> ninit = flatten(nini);
                ret = ret + (currcurrt-ninit).squaredNorm() * w_smooth * w_attenuate;
                // ret = ret + (n_i*n_i.transpose()-currcurr).norm() * w_smooth;
              }
            }
            return ret;
          }



          Eigen::VectorX<T> s_a = element.variables(cur_surf.data().faceNeighbors(f_idx, 0));
          Eigen::VectorX<T> s_b = element.variables(cur_surf.data().faceNeighbors(f_idx, 1));
          Eigen::VectorX<T> s_c = element.variables(cur_surf.data().faceNeighbors(f_idx, 2));

          // Eigen::Vector2<T> s_a_gamma = s_a.segment(2, 2);
          // Eigen::Vector2<T> s_b_gamma = s_b.segment(2, 2);
          // Eigen::Vector2<T> s_c_gamma = s_c.segment(2, 2);

          Eigen::Vector2<T> a = frames.row(cur_surf.data().faceNeighbors(f_idx, 0));
          Eigen::Vector2<T> b = frames.row(cur_surf.data().faceNeighbors(f_idx, 1));
          Eigen::Vector2<T> c = frames.row(cur_surf.data().faceNeighbors(f_idx, 2));

          // a = a + s_a_gamma;
          // b = b + s_b_gamma;
          // c = c + s_c_gamma;

          Eigen::Matrix2<T> aa = a*a.transpose();
          Eigen::Matrix2<T> bb = b*b.transpose();
          Eigen::Matrix2<T> cc = c*c.transpose();




          Eigen::Vector4<T> aat = flatten(aa);
          Eigen::Vector4<T> bbt = flatten(bb);
          Eigen::Vector4<T> cct = flatten(cc);

          Eigen::Vector4<T> s_a_gamma = s_a.head(4);
          Eigen::Vector4<T> s_b_gamma = s_b.head(4);
          Eigen::Vector4<T> s_c_gamma = s_c.head(4);

          aat = aat + s_a_gamma;
          bbt = bbt + s_b_gamma;
          cct = cct + s_c_gamma;



                    // return (T) 0.;
         


          Eigen::Vector4<T> a_delta = s_a.tail(4);
          Eigen::Vector4<T> b_delta = s_b.tail(4);
          Eigen::Vector4<T> c_delta = s_c.tail(4);



          Eigen::Vector2<T> curr_normalized = curr.normalized();
          Eigen::Vector2<T> curr_perp; // = curr_normalized;
          curr_perp(0) = curr_normalized(1);
          curr_perp(1) = -curr_normalized(0);

          // T s_perp_term = pow(a.dot(curr_perp),2) + pow(b.dot(curr_perp),2) + pow(c.dot(curr_perp), 2);

          // T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (curr_perp * curr_perp.transpose())).norm();
          T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (currcurrt)).squaredNorm();

          // T dirichlet_term = (a + b + c - 3*curr).squaredNorm();
          T dirichlet_term = (aat+bbt+cct-3*currcurrt).squaredNorm();

          // T dirichlet_term = (aa + bb + cc - 3*currcurr).norm();
  // dirichlet_term += 1e-5*abs(dirichlet_term - metadata(0));


          T delta_weight = std::min(w_curl/100., 1./w_attenuate);

          T delta_rescale = std::max(frames.row(f_idx).squaredNorm(), 1e-8);
          delta_rescale = (.0001 + 1./delta_rescale);
          T delta_norm_term = delta.squaredNorm();
          T gamma_norm_term = gamma.squaredNorm(); // * delta_rescale;

          T gamma_dirichlet_term = (s_a_gamma+s_b_gamma+s_c_gamma-3*gamma).squaredNorm()*delta_rescale;

          Eigen::Vector4d ea = e_projs2.row(cur_surf.data().faceEdges(f_idx, 0));
          Eigen::Vector4d eb = e_projs2.row(cur_surf.data().faceEdges(f_idx, 1));
          Eigen::Vector4d ec = e_projs2.row(cur_surf.data().faceEdges(f_idx, 2));

          // Eigen::Vector4<T> ea = e_projs.at(cur_surf.data().faceEdges(f_idx, 0));
          // Eigen::Vector4<T> eb = e_projs.at(cur_surf.data().faceEdges(f_idx, 1));
          // Eigen::Vector4<T> ec = e_projs.at(cur_surf.data().faceEdges(f_idx, 2));

          T curl_term = pow(ea.dot(aat + a_delta) - ea.dot(currcurrt + delta),2);
          curl_term +=  pow(eb.dot(bbt + b_delta) - eb.dot(currcurrt + delta),2);
          curl_term +=  pow(ec.dot(cct + c_delta) - ec.dot(currcurrt + delta),2);

          curls(f_idx) = TinyAD::to_passive(curl_term);

          // curl_term += 1e-5*abs(curl_term - metadata(1));


          // T atten = 1./(cur_iter*cur_iter + 1);
          // T atten = 1./(cur_iter + 1);

          // T atten = 1.;


         

          T ret = delta_norm_term * 10000.;
          ret = gamma_norm_term * 10. * w_smooth; // + gamma_dirichlet_term * .1 * w_smooth;
          if (w_smooth > 0)
            ret = ret + w_attenuate * w_smooth * dirichlet_term;
          if (w_s_perp > 0)
            ret = ret + w_attenuate * w_s_perp * s_perp_term;
          if (w_curl > 0)
            ret = ret + w_curl * curl_term;

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
          ret = Eigen::VectorXd::Zero(8); // resize(10);
          // ret.head(4) = Eigen::VectorXd::Random(4);
          // ret << frames.row(f_idx), deltas.row(f_idx);
          return ret;
          });

    }


    virtual void updateRenderGeometry()
    {

      renderFrames.resize(frames.rows(), 3);
      renderFrames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);

      renderDeltas = deltas;
      renderGammas = gammas;

      // sym_curl.resize(frames.rows());
      sym_curl = curls;

      // vec_smoothness.resize(frames.rows());
      vec_smoothness = metadata.col(1);

    }


    virtual bool simulateOneStep()
    {
        if (cur_iter < max_iters)
        {
            cur_iter++;

            auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
            TINYAD_DEBUG_OUT("Energy in iteration " << cur_iter << ": " << f);

            Eigen::VectorXd d = TinyAD::newton_direction(g, H_proj, solver);
            if (TinyAD::newton_decrement(d, g) < convergence_eps)
            {
              w_attenuate = w_attenuate / 10.;
              std::cout << "New attenuation value is set to: " << w_attenuate << std::endl;
              if (w_attenuate < 1e-12)
                 cur_iter = max_iters;
            }
            // cur_iter = max_iters; // break
            x = TinyAD::line_search(x, d, f, g, func, 1., .8, 512, 1e-3);


            ///// Move this out 
            func.x_to_data(x, [&] (int f_idx, const Eigen::VectorXd& v) {
                // frames.row(f_idx) = v.head<2>();
                // metadata.row(f_idx) = v.segment(2, 2);
                // gammas.row(f_idx) = v.segment(2, 2);
                gammas.row(f_idx) = v.head<4>();

                deltas.row(f_idx) = v.tail<4>();
                // if (bound_face_idx(f_idx) == 1)
                // {
                //   frames.row(f_idx) = frames_orig.row(f_idx);
                // }
                });


            int nrows = frames.rows();
            for(int i = 0; i < nrows; i++)
            {
              Eigen::Vector2d v_curr = frames.row(i);
              // Eigen::Vector2d v_curr = gammas.row(i);
              Eigen::Matrix2d vtv_curr = v_curr*v_curr.transpose();
              Eigen::Matrix2d p_curr = vtv_curr;
              // Eigen::VectorXd gamma_curr = jacobians.at(i)*gammas.row(i).transpose();
              Eigen::VectorXd gamma_curr = gammas.row(i);


              Eigen::VectorXd delta_curr = deltas.row(i);
              
              p_curr(0,0) = p_curr(0,0) + gamma_curr(0);// + delta_curr(0);
              p_curr(0,1) = p_curr(0,1) + gamma_curr(1);// + delta_curr(0);
              p_curr(1,0) = p_curr(1,0) + gamma_curr(2);// + delta_curr(0);
              p_curr(1,1) = p_curr(1,1) + gamma_curr(3);// + delta_curr(0);

              // std::cout << "normalized " << v_curr.normalized().transpose() << "norm " << v_curr.norm() << "v_curr " << v_curr.transpose() <<  std::endl;

              GN_proj_to_rank_1(p_curr, v_curr);
              frames.row(i) = v_curr;

              Eigen::Matrix2d vtv_after_proj = v_curr*v_curr.transpose();
              Eigen::Matrix2d gamma_after_proj = p_curr - vtv_after_proj;
              gammas.row(i) = flatten(gamma_after_proj) ;

              // std::cout << gammas.row(i) << std::endl;



              // x = x*0;
              x = func.x_from_data([&] (int f_idx) {
                Eigen::VectorXd ret;
                ret = Eigen::VectorXd::Zero(8); // resize(10);
                ret.tail(4) = deltas.row(f_idx);
                ret.head(4) = gammas.row(f_idx);
                // ret.head(4) = Eigen::VectorXd::Random(4);
                // ret << frames.row(f_idx), deltas.row(f_idx);
                return ret;
              });
              

            }
            updateJacobians(frames, jacobians);


                          // Eigen::Map<Eigen::Matrix2d> gamma_curr(gammas.row(i), 2,2);
              // Eigen::Matrix2d p_curr = v_curr*v_curr.transpose() + 

// Svd2x2Helper(p_curr);
// Eigen::JacobiSVD<Eigen::Matrix2d> svd(p_curr,Eigen::ComputeFullU | Eigen::ComputeFullV);
// std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
// // std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
// std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;

// Vector3f rhs(1, 0, 0);
// cout << "Now consider this rhs vector:" << endl << rhs << endl;
// cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;



              // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>> eigen_solver(p_curr);
              // Eigen::Vector<double, 2> eigenvalues = eigen_solver.eigenvalues().reverse(); // values are OK
              // Eigen::Matrix<double, 2, 2> eigenvectors = eigen_solver.eigenvectors().transpose().colwise().reverse();
              // std::cout << "eigenvectors:\n" << eigenvectors.matrix() << "\n";
              // std::cout << "eigenvectors:\n" << eigenvalues.matrix() << "\n";
              // std::cout << "vtv_curr:\n" << p_curr << "\n";


// Eigen::Vector2d v_curr = frames.row(i);
              // Eigen::Vector2d gamma_debug = gammas.row(i);

              // Eigen::Vector2d curr_row = frames.row(i);
              // curr_row = v_curr + gamma_debug;
              // frames.row(i) = curr_row; // eigenvectors.row(0)*eigenvalues(0);

            // updateJacobians(frames, jacobians);


          //        x = func.x_from_data([&] (int f_idx) {
          // Eigen::VectorXd ret;
          // ret = Eigen::VectorXd::Zero(10); // resize(10);
          // ret.head(2) = Eigen::VectorXd::Random(2);
          // // ret << frames.row(f_idx), deltas.row(f_idx);
          // return ret;
          // });



        }
        else if (cur_iter == max_iters) 
        {
            TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));
            cur_iter++;
        }
        else{
            this->pause();
        }


        
        return false;
    }


// Need to fill out viewer for each of: Field_View { vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual,
    virtual void renderRenderGeometry()
    {
		polyscope::getSurfaceMesh()->updateVertexPositions(renderP);
        
        polyscope::getSurfaceMesh()->centerBoundingBox();
        polyscope::getSurfaceMesh()->resetTransform();

        switch (current_element) 
        { 
            case Field_View::vec_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", renderFrames.rowwise().squaredNorm())->setEnabled(true);
                } 
                break;

            case Field_View::delta_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("delta_norms", renderDeltas.rowwise().squaredNorm())->setEnabled(true);
                } 
                break; 
            case Field_View::gamma_norms:
                { 
                  polyscope::getSurfaceMesh()->addFaceScalarQuantity("delta_norms", renderGammas.rowwise().squaredNorm())->setEnabled(true);
                } 
                break; 
              
            case Field_View::primal_curl_residual: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("primal_curl_residual", renderFrames.rowwise().norm())->setEnabled(true);
                } 
                break;             
            case Field_View::sym_curl_residual: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("sym_curl_residual", sym_curl)->setEnabled(true);
                } 
                break;           
            case Field_View::vec_dirch: 
                { 
                    polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_dirch", vec_smoothness)->setEnabled(true);
                } 
                break; 

            default: 
                { 
                    // std::cout << "Unknown color!"; 
                } 
                break; 
        } 
      
  
        auto vectors = polyscope::getSurfaceMesh()->addFaceVectorQuantity("frames", renderFrames); //   ( ((N.array()*0.5)+0.5).eval());
        vectors->vectorColor = glm::vec3(.7,.7,.7);
        
        // vectors->setVectorLengthScale(1., true);
        vectors->setEnabled(true);
        // vectors->setVectorColor(glm::vec3(.7,.7,.7));
        
        polyscope::requestRedraw();   
    }


// static auto func;
  double w_bound;
  double w_smooth; 
  double w_curl;
  double w_s_perp;

  double w_attenuate;



private:
  // Read mesh and compute Tutte embedding
  Eigen::MatrixXd V; // #V-by-3 3D vertex positions
  Eigen::MatrixXi F; // #F-by-3 indices into V
  Eigen::MatrixXd P; //  = tutte_embedding(V, F); // #V-by-2 2D vertex positions
  Eigen::MatrixXd frames;
  Eigen::MatrixXd deltas;
  Eigen::MatrixXd gammas;
  Eigen::MatrixXd metadata;
  Eigen::MatrixXd frames_orig;

  Eigen::VectorXd curls; 

  Surface cur_surf;


      std::vector<Eigen::Matrix2d> rots;// 
      std::vector<Eigen::Matrix4d> rstars;
      std::vector<Eigen::Vector4d> e_projs;

  std::vector<Eigen::MatrixXd> jacobians;

      Eigen::MatrixXd e_projs2;


  
  Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 

  Eigen::MatrixXd renderFrames;
  Eigen::MatrixXd renderDeltas;
  Eigen::MatrixXd renderGammas;
  Eigen::MatrixXd renderP;
  Eigen::MatrixXi renderF;
  Field_View current_element;

  Eigen::VectorXd vec_curl;
  Eigen::VectorXd sym_curl;
  Eigen::VectorXd vec_norms;
  Eigen::VectorXd delta_norms;
  Eigen::VectorXd vec_smoothness;
  Eigen::VectorXd moment_smoothness;
  // Eigen::VectorXd moment_smoothness;


  
  std::vector<Eigen::Matrix2d> rest_shapes;
// %%% 2 + 4 + 4
  decltype(TinyAD::scalar_function<8>(TinyAD::range(1))) func;
  Eigen::VectorXd x;

  int max_iters = 5000;
  int cur_iter = 0;
  double convergence_eps = 1e-10;

  TinyAD::LinearSolver<double> solver;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;


    
};