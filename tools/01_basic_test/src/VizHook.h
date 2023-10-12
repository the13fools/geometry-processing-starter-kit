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