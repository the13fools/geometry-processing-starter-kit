
#include "VizHelper.h"

#include <Eigen/Core>

using namespace VizHelper;

void VizCache::init_matricies()
{

    int nfaces = nFaces();
    d().frame_norms = Eigen::VectorXd::Zero(nfaces);
    d().moment_norms = Eigen::VectorXd::Zero(nfaces);
    d().delta_norms = Eigen::VectorXd::Zero(nfaces);
    d().gamma_norms = Eigen::VectorXd::Zero(nfaces);

}