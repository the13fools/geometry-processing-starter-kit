#include "PhysicsHook.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <Eigen/Core>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include <igl/readOBJ.h>

#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/per_vertex_normals.h>
#include <igl/map_vertices_to_circle.h>


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

    }

    virtual void initSimulation()
    {

      igl::readOBJ(std::string(SOURCE_PATH) + "/armadillo_cut_low.obj", V, F);
    // igl::readOBJ(std::string(SOURCE_PATH) + "/decimated-max.obj", V, F);
    //  igl::readOBJ(std::string(SOURCE_PATH) + "/inspired_mesh.obj", V, F);

     
      P = tutte_embedding(V, F); 

      std::cout << "v size " <<  V.size() << " f size " << F.size() << " p size " << P.size() <<std::endl;

      cur_iter = 0; 

      polyscope::removeAllStructures();
      renderP.resize(P.rows(), 3);
      renderP << P, Eigen::MatrixXd::Zero(P.rows(), 1);

      renderF = F; 
      polyscope::registerSurfaceMesh("cur state", renderP, renderF);
      polyscope::getSurfaceMesh()->setEdgeWidth(.6);
      polyscope::view::resetCameraToHomeView();

      Eigen::MatrixXd N;
      igl::per_vertex_normals(V,F,N);
      polyscope::getSurfaceMesh()->addVertexColorQuantity("orig normals", ((N.array()*0.5)+0.5).eval())->setEnabled(true); //   ( ((N.array()*0.5)+0.5).eval());

      // Pre-compute triangle rest shapes in local coordinate systems
      rest_shapes.clear();
      rest_shapes.resize(F.rows());
      for (int f_idx = 0; f_idx < F.rows(); ++f_idx)
      {
        // Get 3D vertex positions
        Eigen::Vector3d ar_3d = V.row(F(f_idx, 0));
        Eigen::Vector3d br_3d = V.row(F(f_idx, 1));
        Eigen::Vector3d cr_3d = V.row(F(f_idx, 2));

        // Set up local 2D coordinate system
        Eigen::Vector3d n = (br_3d - ar_3d).cross(cr_3d - ar_3d);
        Eigen::Vector3d b1 = (br_3d - ar_3d).normalized();
        Eigen::Vector3d b2 = n.cross(b1).normalized();

        // Express a, b, c in local 2D coordiante system
        Eigen::Vector2d ar_2d(0.0, 0.0);
        Eigen::Vector2d br_2d((br_3d - ar_3d).dot(b1), 0.0);
        Eigen::Vector2d cr_2d((cr_3d - ar_3d).dot(b1), (cr_3d - ar_3d).dot(b2));

        // Save 2-by-2 matrix with edge vectors as colums
        rest_shapes[f_idx] = TinyAD::col_mat(br_2d - ar_2d, cr_2d - ar_2d);
      };

      // Set up function with 2D vertex positions as variables.
      func = TinyAD::scalar_function<2>(TinyAD::range(V.rows()));

      // Add objective term per face. Each connecting 3 vertices.
      func.add_elements<3>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
          {
          // Evaluate element using either double or TinyAD::Double
          using T = TINYAD_SCALAR_TYPE(element);

          // Get variable 2D vertex positions
          Eigen::Index f_idx = element.handle;
          Eigen::Vector2<T> a = element.variables(F(f_idx, 0));
          Eigen::Vector2<T> b = element.variables(F(f_idx, 1));
          Eigen::Vector2<T> c = element.variables(F(f_idx, 2));

          // Triangle flipped?
          Eigen::Matrix2<T> M = TinyAD::col_mat(b - a, c - a);
          if (M.determinant() <= 0.0)
          return (T)INFINITY;

          // Get constant 2D rest shape of f
          Eigen::Matrix2d Mr = rest_shapes[f_idx];
          double A = 0.5 * Mr.determinant();

          // Compute symmetric Dirichlet energy
          Eigen::Matrix2<T> J = M * Mr.inverse();
          return A * (J.squaredNorm() + J.inverse().squaredNorm());
          });

      // Assemble inital x vector from P matrix.
      // x_from_data(...) takes a lambda function that maps
      // each variable handle (vertex index) to its initial 2D value (Eigen::Vector2d).
        x = func.x_from_data([&] (int v_idx) {
          return P.row(v_idx);
          });

    }

        /**
     * Compute tutte embedding with boundary on circle.
     * Per-vertex 2D coordinates returned as n_vertices-by-2 matrix.
     * 
     * helper function
     */
    Eigen::MatrixXd tutte_embedding(
        const Eigen::MatrixXd& _V,
        const Eigen::MatrixXi& _F)
    {
    Eigen::VectorXi b; // #constr boundary constraint indices
    Eigen::MatrixXd bc; // #constr-by-2 2D boundary constraint positions
    Eigen::MatrixXd P; // #V-by-2 2D vertex positions
    igl::boundary_loop(_F, b); // Identify boundary vertices
    igl::map_vertices_to_circle(_V, b, bc); // Set boundary vertex positions
    igl::harmonic(_F, b, bc, 1, P); // Compute interior vertex positions

    // P.resize(P.rows(),3);

    return P;
    }


    virtual void updateRenderGeometry()
    {
        renderP.resize(P.rows(), 3);
        renderP << P, Eigen::MatrixXd::Zero(P.rows(), 1);
        renderF = F;
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
            func.x_to_data(x, [&] (int v_idx, const Eigen::Vector2d& p) {
                P.row(v_idx) = p;
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
        
        polyscope::requestRedraw();   
    }


// static auto func;

private:
  // Read mesh and compute Tutte embedding
  Eigen::MatrixXd V; // #V-by-3 3D vertex positions
  Eigen::MatrixXi F; // #F-by-3 indices into V
//   igl::readOBJ(std::string(SOURCE_PATH) + "/armadillo_cut_low.obj", V, F);
  Eigen::MatrixXd P; //  = tutte_embedding(V, F); // #V-by-2 2D vertex positions
  
  std::vector<Eigen::Matrix2d> rest_shapes;

  decltype(TinyAD::scalar_function<2>(TinyAD::range(1))) func;
  Eigen::VectorXd x;

  int max_iters = 5000;
  int cur_iter = 0;
  double convergence_eps = 1e-2;

  TinyAD::LinearSolver<double> solver;


// using T = TINYAD_SCALAR_TYPE(element);
// int f_idx = element.handle;
// Eigen::Vector2<T> a = element.variables(F(f_idx, 0));
// Eigen::Vector2<T> b = element.variables(F(f_idx, 1));
// Eigen::Vector2<T> c = element.variables(F(f_idx, 2));
//   TinyAD::ScalarFunction<2,double, int> func;
                       
                                             //
    bool redraw = false;
    // std::mutex m;


    Eigen::MatrixXd renderP;
    Eigen::MatrixXi renderF;
};