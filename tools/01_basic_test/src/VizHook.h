#include "PhysicsHook.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#include <Eigen/LU>

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>
#include <chrono>

class VizHook : public PhysicsHook
{
public:
    VizHook() : PhysicsHook() {}

    virtual void drawGUI()
    {
	//	ImGui::SliderFloat("k scale", &k_scale, 0.0001f, 2.0f, "k * %.3f");
	//	ImGui::SliderFloat("dt scale", &dt_scale, 0.0001f, 2.0f, "dt * %.3f");
		ImGui::InputFloat("k scale", &k_scale);
		ImGui::InputFloat("dt scale", &dt_scale);

    }

    virtual void initSimulation()
    {
        origQ.resize(4, 3);
        origQ << -1, -1, 0,
            1, -1, 0,
            -1, 1, 0,
            1, 1, 0;
        Q = origQ*1.3;
        V.resize(4, 3);
        V.setZero();
        F.resize(2, 3);
        F << 0, 1, 2,
            2, 1, 3;

        dt = 1e-5;
        k = 1e-2;

        polyscope::removeAllStructures();
        renderQ = Q;
        renderF = F; 
        polyscope::registerSurfaceMesh("cur state", renderQ, renderF);
		k_scale = 1.;
		dt_scale = 1.;

    

        /* 

        std::cout << "Test linear solver performance" << std::endl;

        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(0.0,1.0);

        int rows=10000;
        int cols=10000;

        std::vector<Eigen::Triplet<double> > tripletList;
        for(int i=0;i<rows;++i)
        {
            for(int j=0;j<cols;++j)
            {
                auto v_ij=dist(gen);
                auto v_ij_val=dist(gen);     //generate random number
                if(v_ij < 0.01)
                {
                    tripletList.push_back(Eigen::Triplet<double>(i,j,v_ij_val));      //if larger than treshold, insert it
                }
            }
        }
        Eigen::SparseMatrix<double> mat(rows,cols);
        mat.setFromTriplets(tripletList.begin(), tripletList.end()); 
        Eigen::SparseMatrix<double> id(rows,cols);
        id.setIdentity();

        

        Eigen::SparseMatrix<double> test;
        test = mat.transpose(); // + 100. * id;
        test = test + mat + 1000000. * id;
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> eigllt; 
        Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double>> cholmodllt; 

        auto t1 = std::chrono::high_resolution_clock::now();
        eigllt.compute(test);
        auto t2 = std::chrono::high_resolution_clock::now();
        cholmodllt.compute(test);
        auto t3 = std::chrono::high_resolution_clock::now();
        auto eigllt_ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        auto cholmodllt_ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2);
        std::cout << "eig llt took: " << eigllt_ms_int.count() << "ms\n";
        std::cout << "cholmod llt took: " << cholmodllt_ms_int.count() << "ms\n";


        if(eigllt.info()!=Eigen::Success) {
            std::cout << "llt failed" << std::endl;
        }
        
        if(cholmodllt.info()!=Eigen::Success) {
            std::cout << "llt failed" << std::endl;
        }
        */

    }

    virtual void updateRenderGeometry()
    {
        renderQ = Q;
        renderF = F;
    }

    virtual bool simulateOneStep()
    {
		Q += V * (dt * (double)dt_scale);
		Eigen::MatrixXd Force = (origQ - Q)*(k *(double)k_scale);
        V += dt * Force;
	//	std::cout << V << std::endl;
        
        return false;
    }

    virtual void renderRenderGeometry()
    {
		polyscope::getSurfaceMesh()->updateVertexPositions(renderQ);
        // polyscope::getSurfaceMesh()->centerBoundingBox();
        // polyscope::getSurfaceMesh()->resetTransform();
        // polyscope::view::resetCameraToHomeView();

        polyscope::requestRedraw();   
    }

private:
    double k;
    double dt;
	float dt_scale;
	float k_scale;
    Eigen::MatrixXd origQ;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
};