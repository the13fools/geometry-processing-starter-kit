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

#include <TinyAD/ScalarFunction.hh>
// #include <TinyAD/Utils/NewtonDirection.hh>
// #include <TinyAD/Utils/NewtonDecrement.hh>
// #include <TinyAD/Utils/LineSearch.hh>

// // #include <fstream>
// #include <sys/stat.h>
// #include <iostream>


/// A bit of hackery, curtousy of https://stackoverflow.com/a/13980645
//// Having to pass around templates is a bit ugly, and unnecessary in this case.
// decltype(TinyAD::scalar_function<6>(TinyAD::range(1)))
// decltype(TinyAD::scalar_function<6>(TinyAD::range(1)));
// class TinyADFuncBase
// {
// public:
//     virtual ~TinyADFuncBase() {}
//     virtual evaluate() {};
//     template<class T> const T& get() const; //to be implimented after Parameter
//     template<class T, class U> void setValue(const U& rhs); //to be implimented after Parameter
// };

// template <typename T>
// class TinyADFunc : public TinyADFuncBase
// {
// public:
//     TinyADFunc(const T& rhs) :value(rhs) {}
//     const T& get() const {return value;}
//     void setValue(const T& rhs) {value=rhs;}    
// private:
//     T value;
// };

// //Here's the trick: dynamic_cast rather than virtual
// template<class T> const T& TinyADFuncBase::get() const
// { return dynamic_cast<const TinyADFunc<T>&>(*this).get(); }
// template<class T, class U> void TinyADFuncBase::setValue(const U& rhs)
// { return dynamic_cast<TinyADFunc<T>&>(*this).setValue(rhs); }


// // template<typename N> 
// class dummy 




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

    virtual void initSimulation();

    virtual bool simulateOneStep();


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

        // using funcw = decltype(TinyAD::scalar_function<6>(TinyAD::range(1)));


    protected:
        VizHelper::VizCache vc;


        Surface cur_surf;

        Eigen::MatrixXd V; // #V-by-3 3D vertex positions
        Eigen::MatrixXi F; // #F-by-3 indices into V
        Eigen::MatrixXd P; //  = tutte_embedding(V, F); // #V-by-2 2D vertex positions

        Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 


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

        double prev_energy; 
        bool useProjHessian = true;
        Eigen::VectorXd x;

        // funcw func;

          decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;
        // TinyADFuncBase func_wrapper;

        int buffer;

};
    

#endif
