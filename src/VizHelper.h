#ifndef VIZHELPER_H
#define VIZHELPER_H

#include <mutex>
#include <thread>

#include "polyscope/polyscope.h"

#include <Eigen/Core>
#include <vector>

#include "Surface.h"

namespace VizHelper 
{

struct UI_State
{

    bool save_data = true;
    std::string save_path = "prev_data";
    std::string obj_path = "SET THIS IN YOUR TOOL CODE DOOD!";




};


struct VizData
{
    // Surface Data Structures;
    Surface s; // surface data

    UI_State ui;

    // Duplicated here for convenience
    const Eigen::MatrixXd V;
    const Eigen::MatrixXi F;

    // Optimization State 
    Eigen::MatrixXd frames;
    Eigen::MatrixXd moments;
    Eigen::MatrixXd deltas;
    Eigen::MatrixXd gammas;

    // Cache of computed quantites for visualization, not every simulation will use all of these.

    Eigen::VectorXd frame_norms;
    Eigen::VectorXd moment_norms;
    Eigen::VectorXd delta_norms;
    Eigen::VectorXd gamma_norms;

    Eigen::VectorXd frame_smoothness;
    Eigen::VectorXd moment_smoothness;
    Eigen::VectorXd delta_smoothness;
    Eigen::VectorXd gamma_smoothness;

    Eigen::VectorXd vec_curl;
    Eigen::VectorXd sym_curl;

    // Other data for making charts later.
    int step; 

    std::vector<double> objective_fun;
    // std::vector<double> smoothness_weight;
    // std::vector<double> delta_weight;
    // std::vector<double> smoothness_weight;
    // std::vector<double> smoothness_weight;
    

};




// Class wrapping a triangle mesh surface embedded in R^3, along with its combinatorial and geometric data structures
class VizCache
{
public:
    VizCache(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    VizCache(){};
    virtual ~VizCache() {}

    const VizData &data() const { return d; }

    int nVerts() const { return d.s.data().V.rows(); }
    int nFaces() const { return d.s.data().F.rows(); }
    int nEdges() const { return d.s.data().E.rows(); }

    void updateVizState()
    {
        updateVizState(d.frames, d.moments, d.deltas, d.gammas);
    }

    void updateVizState(Eigen::VectorXd frames,     
                        Eigen::MatrixXd moments,
                        Eigen::MatrixXd deltas)
    {
        Eigen::MatrixXd gammas;
        updateVizState(frames, moments, deltas, gammas);
    }

    void updateVizState(Eigen::VectorXd frames,     
                        Eigen::MatrixXd moments,
                        Eigen::MatrixXd deltas,
                        Eigen::MatrixXd gammas)
    {
        int nfaces = nFaces();


        d.frame_norms = frames.rowwise().squaredNorm();
        d.moment_norms = moments.rowwise().squaredNorm();
        d.delta_norms = deltas.rowwise().squaredNorm();
        d.gamma_norms = gammas.rowwise().squaredNorm();



    }
    // int numInteriorEdges() const;

    // Eigen::Vector3d faceNormal(int face) const;
    // double faceArea(int face) const;

    // // Finds shortest (combinatorial) path from start to end vertex. Each path entry is a combination of (1) the edge index along the path, and (2) the orientation: the jth path segment goes from
    // // edgeVerts(path[j].first, path[j].second) to edgeVerts(path[j].first, 1 - path[j].second).
    // // List will be empty if no path exists (vertices lie on disconnected components).
    // void shortestPath(int startVert, int endVert, std::vector<std::pair<int, int> > &path) const;

private:
    // computes E and faceedges/faceWings from V and F
    // void buildConnectivityStructures();
    // // compute the geometric data structs
    // void buildGeometricStructures();

    VizData d;
};

}

#endif
