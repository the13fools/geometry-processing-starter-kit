#ifndef MINT2DHOOK_H
#define MINT2DHOOK_H

#include <mutex>
#include <thread>

#include "PhysicsHook.h"

#include "polyscope/polyscope.h"

#include "polyscope/surface_mesh.h"

#include <Eigen/Core>


#include "VizHelper.h"

// #include "UtilsMisc.h"

#include <TinyAD/Utils/LinearSolver.hh>

// #include <TinyAD/ScalarFunction.hh>
// #include <TinyAD/Utils/NewtonDirection.hh>
// #include <TinyAD/Utils/NewtonDecrement.hh>
// #include <TinyAD/Utils/LineSearch.hh>

// // #include <fstream>
// #include <sys/stat.h>
// #include <iostream>

enum Field_View { vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT };




class Mint2DHook : public virtual PhysicsHook
{
public:

    Mint2DHook() : PhysicsHook() {
    //   current_element = Field_View::vec_norms;
    }



    virtual void drawGUI();



    virtual void updateRenderGeometry();

    virtual void renderRenderGeometry();


    virtual ~Mint2DHook()
    {      
    }

    public: 

        /// TODO: SAVE LOAD stuff 

            // static auto func;
        double w_bound;
        double w_smooth;
        double w_smooth_vector; 
        double w_curl;
        double w_s_perp;

        double w_attenuate;

        Field_View current_element;


    protected:
        VizHelper::VizCache vc;

        Eigen::MatrixXd V; // #V-by-3 3D vertex positions
        Eigen::MatrixXi F; // #F-by-3 indices into V
        Eigen::MatrixXd P; //  = tutte_embedding(V, F); // #V-by-2 2D vertex positions


        Eigen::MatrixXd renderP;
        Eigen::MatrixXi renderF;

        Eigen::MatrixXd frames;
        Eigen::MatrixXd deltas;

        Eigen::MatrixXd frames_orig;

        Eigen::VectorXd curls_sym;
        Eigen::VectorXd curls_primal; 
        Eigen::VectorXd smoothness_primal;
        Eigen::VectorXd smoothness_sym;

        std::string cur_mesh_name;
        std::string cur_log_folder;


        int max_iters = 5000;
        int cur_iter = 0;
        int inner_loop_iter = 0;
        double convergence_eps = 1e-10;
        double identity_weight = 1e-6;

        TinyAD::LinearSolver<double> solver;

        int buffer;

};
    

#endif
