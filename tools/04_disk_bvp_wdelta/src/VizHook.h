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

    virtual void initSimulation()
    {

      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_1000.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole2.obj", V, F);
      // igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_little_hole.obj", V, F);
            igl::readOBJ(std::string(SOURCE_PATH) + "/circle_pent_hole_descimate.obj", V, F);



      cur_surf = Surface(V, F);


     
      // P = tutte_embedding(V, F); 

      // std::cout << "v size " <<  V.size() << " f size " << F.size() << " p size " << P.size() <<std::endl;

      cur_iter = 0; 


      w_bound = 1e3; 
      w_smooth = 1e-5; 
      w_s_perp = 0;
      w_curl = 1e5;
 
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
      deltas = Eigen::MatrixXd::Zero(F.rows(), 2);

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
            frames.row(i) = Eigen::Vector2d(c(1),-c(0));
          }

          
          // std::cout << "i" << i << "bound_face_idx(i)" << bound_face_idx(i) << std::endl;
        }

      }


      frames_orig = frames;

      // std::cout << frames << std::endl;



      renderFrames.resize(frames.rows(), 3);
      renderFrames << frames, Eigen::MatrixXd::Zero(frames.rows(), 1);

      polyscope::getSurfaceMesh()->addFaceVectorQuantity("orig normals", renderFrames); //   ( ((N.array()*0.5)+0.5).eval());
      polyscope::getSurfaceMesh()->addFaceScalarQuantity("vec_norms", frames.rowwise().squaredNorm())->setEnabled(true); //   ( ((N.array()*0.5)+0.5).eval());



      // Set up function with 2D vertex positions as variables.
      func = TinyAD::scalar_function<4>(TinyAD::range(F.rows()));

      // Add objective term per face. Each connecting 3 vertices.
      func.add_elements<4>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
          {
          // Evaluate element using either double or TinyAD::Double
          using T = TINYAD_SCALAR_TYPE(element);



          // Get variable 2D vertex positions
          Eigen::Index f_idx = element.handle;
          Eigen::Vector4<T> s_curr = element.variables(f_idx);
          Eigen::Vector2<T> curr =  s_curr.head(2);
          Eigen::Vector2<T> delta = s_curr.tail(2);
          
          if (bound_face_idx(f_idx) == 1)
          {

            Eigen::Vector2<T> targ = frames_orig.row(f_idx);
            return w_bound*(curr-targ).squaredNorm() + w_bound*delta.squaredNorm();
          }

          if (bound_face_idx(f_idx) == -1)
          {
            return (T) 0;
          }
         


          Eigen::Vector4<T> s_a = element.variables(cur_surf.data().faceNeighbors(f_idx, 0));
          Eigen::Vector4<T> s_b = element.variables(cur_surf.data().faceNeighbors(f_idx, 1));
          Eigen::Vector4<T> s_c = element.variables(cur_surf.data().faceNeighbors(f_idx, 2));

          Eigen::Vector2<T> a = s_a.head(2);
          Eigen::Vector2<T> b = s_b.head(2);
          Eigen::Vector2<T> c = s_c.head(2);


          Eigen::Vector2<T> a_delta = s_a.tail(2);
          Eigen::Vector2<T> b_delta = s_b.tail(2);
          Eigen::Vector2<T> c_delta = s_c.tail(2);



          Eigen::Vector2<T> curr_normalized = curr.normalized();
          Eigen::Vector2<T> curr_perp; // = curr_normalized;
          curr_perp(0) = curr_normalized(1);
          curr_perp(1) = -curr_normalized(0);

          T s_perp_term = pow(a.dot(curr_perp),2) + pow(b.dot(curr_perp),2) + pow(c.dot(curr_perp), 2);

          T dirichlet_term = (a + b + c - 3*curr).squaredNorm();

          T delta_norm_term = delta.squaredNorm();

          Eigen::Vector2i ea_idx = cur_surf.data().edgeVerts.row(cur_surf.data().faceEdges(f_idx, 0));
          Eigen::Vector2i eb_idx = cur_surf.data().edgeVerts.row(cur_surf.data().faceEdges(f_idx, 1));
          Eigen::Vector2i ec_idx = cur_surf.data().edgeVerts.row(cur_surf.data().faceEdges(f_idx, 2));

          Eigen::Vector2<T> ea = (V.row(ea_idx(0)) - V.row(ea_idx(1))).head<2>();
          Eigen::Vector2<T> eb = (V.row(eb_idx(0)) - V.row(eb_idx(1))).head<2>();
          Eigen::Vector2<T> ec = (V.row(ec_idx(0)) - V.row(ec_idx(1))).head<2>();

          T curl_term = pow(ea.dot(a + a_delta) - ea.dot(curr + delta),2);
          curl_term +=  pow(eb.dot(b + b_delta) - eb.dot(curr + delta),2);
          curl_term +=  pow(ec.dot(c + c_delta) - ec.dot(curr + delta),2);

          T atten = 1./(cur_iter + 1);

          return (w_smooth*dirichlet_term + 
                 w_s_perp * s_perp_term) * atten + 
                 w_curl*curl_term  + 
                 delta_norm_term;
          });

      // Assemble inital x vector from P matrix.
      // x_from_data(...) takes a lambda function that maps
      // each variable handle (vertex index) to its initial 2D value (Eigen::Vector2d).
        x = func.x_from_data([&] (int f_idx) {
          Eigen::Vector4d ret;
          // ret << frames.row(f_idx), deltas.row(f_idx);
          return ret;
          });

    }


    virtual void updateRenderGeometry()
    {

      renderFrames.resize(frames.rows(), 3);
      renderFrames << frames + deltas, Eigen::MatrixXd::Zero(frames.rows(), 1);

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
            x = TinyAD::line_search(x, d, f, g, func);


            ///// Move this out 
            func.x_to_data(x, [&] (int f_idx, const Eigen::Vector4d& v) {
                frames.row(f_idx) = v.head<2>();
                deltas.row(f_idx) = v.tail<2>();
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


  
  Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 

  Eigen::MatrixXd renderFrames;
  Eigen::MatrixXd renderP;
  Eigen::MatrixXi renderF;

  
  std::vector<Eigen::Matrix2d> rest_shapes;

  decltype(TinyAD::scalar_function<4>(TinyAD::range(1))) func;
  Eigen::VectorXd x;

  int max_iters = 5000;
  int cur_iter = 0;
  double convergence_eps = 1e-8;

  TinyAD::LinearSolver<double> solver;
  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;


    
};