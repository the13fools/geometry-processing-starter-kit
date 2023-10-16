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

class VizHook : public PhysicsHook
{
public:
    VizHook() : PhysicsHook() {}

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

    }

// Helper fn
    // Eigen::Vector4<TinyAD::Scalar<24, double, true>> flatten(Eigen::Matrix2<TinyAD::Scalar<24, double, true>> r)
    // {
    //   Eigen::Vector4<TinyAD::Scalar<24, double, true>> ret;
    //   ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
    //   return ret;
    // }

template <typename ScalarType, int Rows, int Cols>
Eigen::Matrix<ScalarType, Rows * Cols, 1> flatten(const Eigen::Matrix<ScalarType, Rows, Cols>& matrix) {
    Eigen::Matrix<ScalarType, Rows * Cols, 1> flattened;
    flattened << matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1);
    return flattened;
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
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_1000.obj", V, F);
      igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole2.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_little_hole.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole_descimate.obj", V, F);


      cur_surf = Surface(V, F);


     
      // P = tutte_embedding(V, F); 

      // std::cout << "v size " <<  V.size() << " f size " << F.size() << " p size " << P.size() <<std::endl;

      cur_iter = 0; 


      // w_bound = 1e3; 
      // w_smooth = 1e-5; 
      // w_s_perp = 0;
      // w_curl = 1e5;

      w_bound = 1e4; 
      w_smooth = 0; // 1e4; // 1e3; 
      w_s_perp = 1; // 1e1 
      w_curl = 1e3;
 
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

      frames = Eigen::MatrixXd::Random(F.rows(), 2);


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
            // frames.row(i) = Eigen::Vector2d(c(1),-c(0)); // circulation 
            frames.row(i) = Eigen::Vector2d(c(0),c(1)); // diverging

          }

          
          // std::cout << "i" << i << "bound_face_idx(i)" << bound_face_idx(i) << std::endl;
        }

      }

      int nedges = cur_surf.nEdges();
      // std::vector<Eigen::Matrix2d> rots;// 
      // std::vector<Eigen::Matrix4d> rstars;
      // std::vector<Eigen::Vector4d> e_projs;

      // Eigen::MatrixXd e_projs2;
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
          
          if (bound_face_idx(f_idx) == 1)
          {

            Eigen::Vector2<T> targ = frames_orig.row(f_idx);
            return w_bound*(curr-targ).squaredNorm() + w_bound*delta.squaredNorm();
          }

          if (bound_face_idx(f_idx) == -1)
          {
            return (T) 0;
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

          Eigen::Matrix2<T> currcurr = curr*curr.transpose();


          Eigen::Vector4<T> aat = flatten(aa);
          Eigen::Vector4<T> bbt = flatten(bb);
          Eigen::Vector4<T> cct = flatten(cc);
          Eigen::Vector4<T> currcurrt = flatten(currcurr);


                    // return (T) 0.;
         


          Eigen::Vector4<T> a_delta = s_a.tail(4);
          Eigen::Vector4<T> b_delta = s_b.tail(4);
          Eigen::Vector4<T> c_delta = s_c.tail(4);



          Eigen::Vector2<T> curr_normalized = curr.normalized();
          Eigen::Vector2<T> curr_perp; // = curr_normalized;
          curr_perp(0) = curr_normalized(1);
          curr_perp(1) = -curr_normalized(0);

          // T s_perp_term = pow(a.dot(curr_perp),2) + pow(b.dot(curr_perp),2) + pow(c.dot(curr_perp), 2);

          T s_perp_term = ((a.dot(curr_perp) + b.dot(curr_perp) + c.dot(curr_perp)) * (curr_perp * curr_perp.transpose())).norm();


          T dirichlet_term = (a + b + c - 3*curr).squaredNorm();

          // T dirichlet_term = (aa + bb + cc - 3*currcurr).norm();

          T delta_norm_term = delta.squaredNorm();

          // // Eigen::Vector2i ea_idx = cur_surf.data().edgeVerts.row(cur_surf.data().faceEdges(f_idx, 0));
          // // Eigen::Vector2i eb_idx = cur_surf.data().edgeVerts.row(cur_surf.data().faceEdges(f_idx, 1));
          // // Eigen::Vector2i ec_idx = cur_surf.data().edgeVerts.row(cur_surf.data().faceEdges(f_idx, 2));

          // // Eigen::Vector2<T> ea = (V.row(ea_idx(0)) - V.row(ea_idx(1))).head<2>();
          // // Eigen::Vector2<T> eb = (V.row(eb_idx(0)) - V.row(eb_idx(1))).head<2>();
          // // Eigen::Vector2<T> ec = (V.row(ec_idx(0)) - V.row(ec_idx(1))).head<2>();

          // // T curl_term = pow(ea.dot(a + a_delta) - ea.dot(curr + delta),2);
          // // curl_term +=  pow(eb.dot(b + b_delta) - eb.dot(curr + delta),2);
          // // curl_term +=  pow(ec.dot(c + c_delta) - ec.dot(curr + delta),2);



          // // Eigen::Vector4<T> ea = e_projs2.row(cur_surf.data().faceEdges(f_idx, 0));
          // // Eigen::Vector4<T> eb = e_projs2.row(cur_surf.data().faceEdges(f_idx, 1));
          // // Eigen::Vector4<T> ec = e_projs2.row(cur_surf.data().faceEdges(f_idx, 2));
          // // int tmp = cur_surf.data().faceEdges(f_idx, 0);
          // // std::cout << "target row" << e_projs.at(tmp) << std::endl<< std::endl;

          // // std::cout << "row contents" << e_projs2.row(tmp) << std::endl;


          Eigen::Vector4d ea = e_projs2.row(cur_surf.data().faceEdges(f_idx, 0));
          Eigen::Vector4d eb = e_projs2.row(cur_surf.data().faceEdges(f_idx, 1));
          Eigen::Vector4d ec = e_projs2.row(cur_surf.data().faceEdges(f_idx, 2));

          // Eigen::Vector4<T> ea = e_projs.at(cur_surf.data().faceEdges(f_idx, 0));
          // Eigen::Vector4<T> eb = e_projs.at(cur_surf.data().faceEdges(f_idx, 1));
          // Eigen::Vector4<T> ec = e_projs.at(cur_surf.data().faceEdges(f_idx, 2));

          T curl_term = pow(ea.dot(aat + a_delta) - ea.dot(currcurrt + delta),2);
          curl_term +=  pow(eb.dot(bbt + b_delta) - eb.dot(currcurrt + delta),2);
          curl_term +=  pow(ec.dot(cct + c_delta) - ec.dot(currcurrt + delta),2);


          // T atten = 1./(cur_iter*cur_iter + 1);
          T atten = 1./(cur_iter + 1);

          // T atten = 1.;

return (w_smooth * dirichlet_term + 
                  w_s_perp * s_perp_term) * atten + 
                  w_curl*curl_term  + 
                 delta_norm_term;

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
          ret.resize(6);
          ret.head(2) = Eigen::VectorXd::Random(2);
          // ret << frames.row(f_idx), deltas.row(f_idx);
          return ret;
          });

    }


    virtual void updateRenderGeometry()
    {

      renderFrames.resize(frames.rows(), 3);
      renderFrames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);

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
            cur_iter = max_iters; // break
            x = TinyAD::line_search(x, d, f, g, func, 1., .8, 512, 1e-8);


            ///// Move this out 
            func.x_to_data(x, [&] (int f_idx, const Eigen::VectorXd& v) {
                frames.row(f_idx) = v.head<2>();
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
        }
        else{
            this->pause();
        }


        
        return false;
    }

    virtual void renderRenderGeometry()
    {
		polyscope::getSurfaceMesh()->updateVertexPositions(renderP);
        
        polyscope::getSurfaceMesh()->centerBoundingBox();
        polyscope::getSurfaceMesh()->resetTransform();
        polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", renderFrames.rowwise().squaredNorm())->setEnabled(true);
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



private:
  // Read mesh and compute Tutte embedding
  Eigen::MatrixXd V; // #V-by-3 3D vertex positions
  Eigen::MatrixXi F; // #F-by-3 indices into V
  Eigen::MatrixXd P; //  = tutte_embedding(V, F); // #V-by-2 2D vertex positions
  Eigen::MatrixXd frames;
  Eigen::MatrixXd deltas;
  Eigen::MatrixXd frames_orig;

  Surface cur_surf;


      std::vector<Eigen::Matrix2d> rots;// 
      std::vector<Eigen::Matrix4d> rstars;
      std::vector<Eigen::Vector4d> e_projs;

      Eigen::MatrixXd e_projs2;


  
  Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 

  Eigen::MatrixXd renderFrames;
  Eigen::MatrixXd renderP;
  Eigen::MatrixXi renderF;

  
  std::vector<Eigen::Matrix2d> rest_shapes;

  decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;
  Eigen::VectorXd x;

  int max_iters = 5000;
  int cur_iter = 0;
  double convergence_eps = 1e-8;

  TinyAD::LinearSolver<double> solver;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;


    
};