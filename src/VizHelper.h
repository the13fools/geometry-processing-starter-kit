#ifndef VIZHELPER_H
#define VIZHELPER_H

#include <mutex>
#include <thread>

// #include "polyscope/polyscope.h"

#include <Eigen/Core>
#include <vector>

#include "Surface.h"

namespace VizHelper 
{

struct UI_State
{

    bool save_data = true;
    bool show_frames = true;
    std::string save_path = "prev_data";
    std::string obj_path = "SET THIS IN YOUR TOOL CODE DOOD!";





};


struct VizData
{
    // Surface Data Structures;
    Surface s; // surface data

    UI_State ui;

    // Duplicated here for convenience
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    int rank = 1;
    int nmoments = 0;
    int vec_dim = 2;

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

// enum Field_View { vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT };



// Class wrapping a triangle mesh surface embedded in R^3, along with its combinatorial and geometric data structures
class VizCache
{
public:
    VizCache(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
    {
        VizCache(V,F,0,1,2);
    };
    VizCache(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int nmoments, int rank, int vec_dim)
    {
        d().V = V;
        d().F = F;
        d().nmoments = nmoments;
        d().rank = rank;
        d().vec_dim = vec_dim;

        init_matricies();
    };
    VizCache(){};
    virtual ~VizCache() {}

    void init_matricies();

    VizData &d() { return data_; }

    int nVerts()  { return d().s.data().V.rows(); }
    int nFaces()  { return d().s.data().F.rows(); }
    int nEdges()  { return d().s.data().E.rows(); }

    Eigen::MatrixXd getFrames() { return d().frames; }

    void updateVizState()
    {
        // std::cout << "frames " << d().frames.row(0) << " deltas " << d().deltas.row(0) << std::endl;
        // updateVizState(d().frames, d().moments, d().deltas, d().gammas);
            updateVizState(d().frames, d().frames * 0, d().deltas, d().frames * 0);
    }

    void updateVizState(Eigen::MatrixXd frames,     
                        Eigen::MatrixXd moments,
                        Eigen::MatrixXd deltas)
    {
        Eigen::MatrixXd gammas;
        updateVizState(frames, moments, deltas, gammas);
    }

    void updateVizState(Eigen::MatrixXd frames,     
                        Eigen::MatrixXd moments,
                        Eigen::MatrixXd deltas,
                        Eigen::MatrixXd gammas)
    {
        int nfaces = nFaces();

        // std::cout << "frames " << frames.row(0) << " deltas " << deltas.row(0) << std::endl;

        d().frame_norms = frames.rowwise().squaredNorm();
        // d().moment_norms = moments.rowwise().squaredNorm();
        d().delta_norms = deltas.rowwise().squaredNorm();
        // d().gamma_norms = gammas.rowwise().squaredNorm();



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

    VizData data_;
};

}

#endif
