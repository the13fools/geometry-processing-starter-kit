#include "PhysicsHook.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <Eigen/Core>

class VizHook : public PhysicsHook
{
public:
    VizHook() : PhysicsHook() {}

    virtual void drawGUI()
    {

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

        dt = 1e-2;
        k = 1e-1;

        polyscope::removeAllStructures();
        renderQ = Q;
        renderF = F; 
        polyscope::registerSurfaceMesh("cur state", renderQ, renderF);
    }

    virtual void updateRenderGeometry()
    {
        renderQ = Q;
        renderF = F;
        polyscope::getSurfaceMesh()->updateVertexPositions(renderQ);
    }

    virtual bool simulateOneStep()
    {
        Q += dt * V;
        Eigen::MatrixXd Force = k*(origQ - Q);
        V += dt * Force;
  //      std::cout << V << std::endl;
        
        return false;
    }

 //    virtual void renderRenderGeometry()
 //    {
 //    //    viewer.data().set_mesh(renderQ, renderF);
	// //	renderQ = .9 * renderQ;
	// 	polyscope::getSurfaceMesh()->updateVertexPositions<Eigen::MatrixXd>(renderQ);
	// 	std::cout << "here" << std::endl;
			
	// 	//	updateVertexPositions(renderQ);
 //        polyscope::requestRedraw();   
 //    }

private:
    double k;
    double dt;
    Eigen::MatrixXd origQ;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
};