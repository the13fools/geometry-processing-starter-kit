#ifndef ADHOOK_H
#define ADHOOK_H

// #include "polyscope/polyscope.h"

#include <TinyAD/Utils/LinearSolver.hh>


#include <Eigen/Core>

class ADFuncRunner 
{ 
    virtual double eval_func_local(const Eigen::VectorXd &x) = 0; 
    
    virtual void eval_func_with_derivatives(const Eigen::VectorXd &x) = 0; 
    virtual void eval_func_and_proj_hess_to_psd_local(const Eigen::VectorXd &x) = 0; 

    // virtual Eigen::VectorXd eval_grad_local(const Eigen::VectorXd &x) = 0; 
    // virtual Eigen::MatrixXd eval_hess_local(const Eigen::VectorXd &x) = 0; 

    double eval_func_at(const Eigen::VectorXd &x) { return eval_func_local(x); } 

// can copy tinyad and return a tuple here for convenience.
    void eval_func_with_derivatives_at(const Eigen::VectorXd &x) { eval_func_with_derivatives(x);  } 
    void eval_func_and_proj_hess_to_psd_at(const Eigen::VectorXd &x) { eval_func_and_proj_hess_to_psd_local(x);  } 
    void reset_diag_hessian_reg() { identity_weight = 1e-8; }


    Eigen::VectorXd get_current_dofs_x() { return _cur_x; } 
    double get_fval_at_x() { return _fun_val; } 
    Eigen::VectorXd get_grad_at_x() { return _grad; } 
    Eigen::SparseMatrix<double> get_hessian_at_x() { return _hess; } 


    void take_newton_step(const Eigen::VectorXd &x);


    public: 
        bool useProjHessian = true;
        double prev_energy = -100000.;
        double identity_weight = 1e-8;  // This term controls how much identity we add to the hessian to make it psd.  



    protected: 

        // current state  
        Eigen::VectorXd _cur_x;

        // runner state 
        bool x_curr_is_new = false;
        TinyAD::LinearSolver<double> solver; // make this changable 


        // config flags 


        // cached quantities 

        double _fun_val; 
        Eigen::VectorXd _grad;
        Eigen::SparseMatrix<double> _hess;
        Eigen::SparseMatrix<double> _hess_proj;
}; 


/// This class wraps your tiny ad func in a templated way to allow for taking newton steps in sequence
class ADFunc_TinyAD_Instance : public ADFuncRunner { 
            // ScalarFunction _f; public: 
            // FuncEvaluator(f) : _f(f) {} 
            // In practice: use f to evaluate the function, gradient, and hessian 
            
    double eval_func_local(const VectorXd &x) override { return x.squaredNorm(); } 
    VectorXd eval_grad_local(const VectorXd &x) override { return Matrix::Ones(); } 
    MatrixXd eval_hess_local(const VectorXd &x) override { return Matrix::Identity(x.size(), x.size()); } 


    decltype(TinyAD::scalar_function<N>(TinyAD::range(1))) _func;


}; 
            

///// TODO:  Add enzyme connector here  
// template <size_t N> class ADFunc_enzyme_Instance : public ADFuncRunner { 
//             // ScalarFunction _f; public: 
//             // FuncEvaluator(f) : _f(f) {} 
//             // In practice: use f to evaluate the function, gradient, and hessian 
            
//     double eval_func_local(const VectorXd &x) override { return x.squaredNorm(); } 
//     VectorXd eval_grad_local(const VectorXd &x) override { return Matrix::Ones(); } 
//     MatrixXd eval_hess_local(const VectorXd &x) override { return Matrix::Identity(x.size(), x.size()); } 


//     decltype(TinyAD::scalar_function<N>(TinyAD::range(1))) _func;


// }; 




#endif


    // void clear_state()
    // {
    //     _cur_x = _cur_x * 0;
    //     _fun_val = -1;
    //     _grad.clear();
    //     _hess.clear();
    // }

    // double eval_func(const Eigen::VectorXd &x) { return eval_func_local(x); } 
    // Eigen::VectorXd eval_grad(const Eigen::VectorXd &x) { return eval_grad_local(x); } 
    // Eigen::MatrixXd eval_hess(const Eigen::VectorXd &x) { return eval_hess_local(x); } 
        
    // double eval_func(const Eigen::VectorXd &x) { return eval_func_local(x); } 
    // Eigen::VectorXd eval_grad(const Eigen::VectorXd &x) { return eval_grad_local(x); } 
    // Eigen::MatrixXd eval_hess(const Eigen::VectorXd &x) { return eval_hess_local(x); } 
