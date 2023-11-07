#ifndef UTILSMISC
#define UTILSMISC
#include <Eigen/Core>

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
//////
/// A robust 2x2 svd apparently.  Not well tested.  
/////

// https://scicomp.stackexchange.com/a/28506
template <class T>
void Rq2x2Helper(const Eigen::Matrix<T, 2, 2>& A, T& x, T& y, T& z, T& c2, T& s2) {
    T a = A(0, 0);
    T b = A(0, 1);
    T c = A(1, 0);
    T d = A(1, 1);

    if (c == 0) {
        x = a;
        y = b;
        z = d;
        c2 = 1;
        s2 = 0;
        return;
    }
    T maxden = std::max(abs(c), abs(d));

    T rcmaxden = 1/maxden;
    c *= rcmaxden;
    d *= rcmaxden;

    T den = 1/sqrt(c*c + d*d);

    T numx = (-b*c + a*d);
    T numy = (a*c + b*d);
    x = numx * den;
    y = numy * den;
    z = maxden/den;

    s2 = -c * den;
    c2 = d * den;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template <class T>
void Svd2x2Helper(const Eigen::Matrix<T, 2, 2>& A, T& c1, T& s1, T& c2, T& s2, T& d1, T& d2) {
    // Calculate RQ decomposition of A
    T x, y, z;
    Rq2x2Helper(A, x, y, z, c2, s2);

    // Calculate tangent of rotation on R[x,y;0,z] to diagonalize R^T*R
    T scaler = T(1)/std::max(abs(x), abs(y));
    T x_ = x*scaler, y_ = y*scaler, z_ = z*scaler;
    T numer = ((z_-x_)*(z_+x_)) + y_*y_;
    T gamma = x_*y_;
    gamma = numer == 0 ? 1 : gamma;
    T zeta = numer/gamma;

    T t = 2*sgn(zeta)/(abs(zeta) + sqrt(zeta*zeta+4));

    // Calculate sines and cosines
    c1 = T(1) / sqrt(T(1) + t*t);
    s1 = c1*t;

    // Calculate U*S = R*R(c1,s1)
    T usa = c1*x - s1*y; 
    T usb = s1*x + c1*y;
    T usc = -s1*z;
    T usd = c1*z;

    // Update V = R(c1,s1)^T*Q
    t = c1*c2 + s1*s2;
    s2 = c2*s1 - c1*s2;
    c2 = t;

    // Separate U and S
    d1 = std::hypot(usa, usc);
    d2 = std::hypot(usb, usd);
    T dmax = std::max(d1, d2);
    T usmax1 = d2 > d1 ? usd : usa;
    T usmax2 = d2 > d1 ? usb : -usc;

    T signd1 = sgn(x*z);
    dmax *= d2 > d1 ? signd1 : 1;
    d2 *= signd1;
    T rcpdmax = 1/dmax;

    c1 = dmax != T(0) ? usmax1 * rcpdmax : T(1);
    s1 = dmax != T(0) ? usmax2 * rcpdmax : T(0);
}

template <class T>
void Svd2x2Helper(const Eigen::Matrix<T, 2, 2>& A) {
    // Calculate RQ decomposition of A
    T c1;
    T s1; 
    T c2; 
    T s2;
    T d1; 
    T d2;

    Svd2x2Helper(A, c1, s1, c2, s2, d1, d2);

    std::cout << "vector " << c1 << " " << s1 << " sing value " << d1 << " sqrt " << std::sqrt(d1) << std::endl;
    std::cout << "vector " << c2 << " " << s2 << " sing value " << d2 << " sqrt " << std::sqrt(std::abs(d2)) << std::endl;

}


///////
/// Gauss Newton rank-1 2x2 Projection 
///////

void GN_proj_to_rank_1(Eigen::Matrix2d p, Eigen::Vector2d& v)
{        
    // std::cout << "Newtons method to do rank-1 projection test" << std::endl;

    // Eigen::Matrix2d p = targ*targ.transpose();
        // std::cout << " p " << p << std::endl;
        double a = p(0,0);
        double b = p(1,0);
        double c = p(0,1);
        double d = p(1,1);
        Eigen::Vector2d grad;
        Eigen::Matrix2d hess;
        // Eigen::Vector2d v = Eigen::Vector2d::Random()*10;
        // v = .1*v + targ;

        Eigen::Vector2d v_prev;
        int iter;
        
        for(int i = 0; i < 500; i++)
        {
            // std::cout << v.transpose() << std::endl;
            v_prev = v;

            grad(0) = -4 * v(0) * (a - v(0)*v(0)) - 2 * v(1) * (b - v(0) * v(1)) - 2 * v(1) * (c - v(0) * v(1));
            grad(1) = -2 * v(0) * (b - v(0) * v(1)) - 2 * v(0) * (c - v(0) * v(1)) - 4 * v(1) * (d - v(1)*v(1));
            hess(0,0) = 8 * v(0)*v(0) - 4 * (a - v(0)*v(0)) + 4 * v(1)*v(1); 
            hess(0,1) = 4 * v(0) * v(1) - 2 * (b - v(0) * v(1)) - 2 * (c - v(0) * v(1));
            hess(1,0) = 4 * v(0) * v(1) - 2 * (b - v(0) * v(1)) - 2 * (c - v(0) * v(1));
            hess(1,1) = 4 * v(0)*v(0) + 8 * v(1)*v(1) - 4 * (d - v(1)*v(1));
            
            hess.colPivHouseholderQr().solve(-grad);
            Eigen::Matrix2d inv = hess.inverse();
            v += hess.colPivHouseholderQr().solve(-grad); // .5 * (inv * grad);
            // v += inv*(-grad-hess*v);

            if ( ( v - v_prev ).squaredNorm() < 1e-8 )
            {
                // std::cout << "break at " << i << std::endl;
                break;
            }
      
        }

        std::cout << "p \n " << p << " v \n" << v << v.transpose() << std::endl;
        if(v.norm() < 1e-4)
        {
            v = Eigen::Vector2d::Random();
        }
}

///////////
///// 2x2 matrix ops 
///////////

template <typename ScalarType, int Rows, int Cols>
Eigen::Matrix<ScalarType, Rows * Cols, 1> flatten(const Eigen::Matrix<ScalarType, Rows, Cols>& matrix) {
    Eigen::Matrix<ScalarType, Rows * Cols, 1> flattened;
    flattened << matrix(0, 0), matrix(0, 1), matrix(1, 0), matrix(1, 1);
    return flattened;
}


    Eigen::Matrix4d rstar_from_r(Eigen::Matrix2d r)
    {
      Eigen::Matrix4d ret;
      double a = r(0,0);
      double b = r(0,1);
      double c = r(1,0);
      double d = r(1,1);
      ret << a*a, a*b, a*b, b*b,
             a*c, b*c, a*d, b*d,
             a*c, b*c, a*d, b*d,
             c*c, c*d, c*d, d*d;

      // ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
      return ret;
    }


    Eigen::Vector4d rstar_xcomp_from_r(Eigen::Matrix2d r)
    {
      Eigen::Vector4d ret;
      double a = r(0,0);
      double b = r(0,1);
      double c = r(1,0);
      double d = r(1,1);
      ret << a*a, a*b, a*b, b*b;

      // ret << r(0,0), r(0, 1), r(1,0), r(1,1); // could maybe also use v1.resize(1, 4); possibly faster
      return ret;
    }

////////////
////  Serialize/Deserialize 
////////////

template<typename Derived>
void serializeEigenMatrix(const Eigen::MatrixBase<Derived>& matrix, const std::string& filename) {
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if (!out.is_open()) {
        throw std::runtime_error("Unable to open file for writing");
    }

    // Write the matrix dimensions
    typename Derived::Index rows = matrix.rows(), cols = matrix.cols();
    out.write(reinterpret_cast<const char*>(&rows), sizeof(typename Derived::Index));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(typename Derived::Index));

    // Write the matrix data
    out.write(reinterpret_cast<const char*>(matrix.data()), rows * cols * sizeof(typename Derived::Scalar));

    out.close();
}


template<typename MatrixType>
MatrixType deserializeEigenMatrix(const std::string& filename) {
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("Unable to open file for reading");
    }

    // Read the matrix dimensions
    typename MatrixType::Index rows = 0, cols = 0;
    in.read(reinterpret_cast<char*>(&rows), sizeof(typename MatrixType::Index));
    in.read(reinterpret_cast<char*>(&cols), sizeof(typename MatrixType::Index));

    // Allocate matrix
    MatrixType matrix(rows, cols);

    // Read the matrix data
    in.read(reinterpret_cast<char*>(matrix.data()), rows * cols * sizeof(typename MatrixType::Scalar));

    in.close();

    return matrix;
}



//     /////// Simulation Step 

//     void take_one_step()
//     {

//             auto t1 = std::chrono::high_resolution_clock::now();


//             // auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
//             auto [f, g, H_proj] = func.eval_with_derivatives(x);
//             TINYAD_DEBUG_OUT("Energy in iteration " << cur_iter << ": " << f);
//             // std::cout<<"the number of nonzeros "<<H_proj.nonZeros() << "number of non-zeros per dof " << H_proj.nonZeros() / (6*F.rows()) << " # rows " << H_proj.rows() << " faces " << F.rows() <<std::endl;

//             // std::cout<<"the number of nonzeros "<<H_proj.nonZeros()<<std::endl;
//             Eigen::VectorXd d;
//             double dec;
//             // d = TinyAD::newton_direction(g, H_proj, solver);
//              // = TinyAD::newton_decrement(d, g);

//             if (prev_energy < 0)
//             {
//               prev_energy = f + 100 * convergence_eps;
//             }
            
//             try
//             {
//               if (w_smooth_vector > 0 || useProjHessian)
//               {
//                 auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
//                 f = f_h;
//                 g = g_h;
//                 H_proj = H_proj_h;
//                 d = TinyAD::newton_direction(g, H_proj, solver, 0.);
//                 dec = TinyAD::newton_decrement(d, g);

//                 if ( dec / f < 1e-3)
//                 {
//                   useProjHessian = false;
//                   std::cout << "switch off projected hessian to fine-tune result" << std::endl;
//                 }


//               }
//               else
//               {
//                 d = TinyAD::newton_direction(g, H_proj, solver, identity_weight);
//                 dec = TinyAD::newton_decrement(d, g);
//                 identity_weight = identity_weight / 2.;
//               }
              
//             }
//             catch(const std::exception& e)
//             {
//               auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
//               f = f_h;
//               g = g_h;
//               H_proj = H_proj_h;
//               d = TinyAD::newton_direction(g, H_proj, solver);
//               dec = TinyAD::newton_decrement(d, g);
//               if ( !useProjHessian )
//                 identity_weight = identity_weight * 10.;
//             }
            
//             // 
//             std::cout << "current decrement: " << dec << std::endl;
//             // if( dec < convergence_eps )
//             // {
//             //   buffer -= 1;
//             //   identity_weight = identity_weight * 10.;
//             // }

//             if ( dec < convergence_eps || (inner_loop_iter > 300 && dec / f < 1e-5))
//             {
//               std::cout << "***** current decrement: " << dec << std::endl;
//               buffer = 5;
//               identity_weight = 1e-6;
//               if (w_smooth_vector > 0)
//               {
//                 w_smooth_vector = 0;
                
//               }
//               else {
//                 w_attenuate = w_attenuate / 10.;
//                 std::cout << "New attenuation value is set to: " << w_attenuate << std::endl;
//                 // inner_loop_iter = 0;
//                 if (w_attenuate < 1e-12)
//                   cur_iter = max_iters;
//               }

//               useProjHessian = true;
                 

//                 // Eigen::MatrixXd tmp =TinyAD::to_passive(H_proj);
//                 //  igl::writeDMAT("converged_hessian.dmat",tmp,true);
//             }

//             // Eigen::MatrixXd tmp =TinyAD::to_passive(H_proj);
//             // igl::writeDMAT("curr_hessian.dmat",tmp,true);

//             // cur_iter = max_iters; // break
//             x = TinyAD::line_search(x, d, f, g, func, 1., .8, 512, 1e-3);

//             auto t2 = std::chrono::high_resolution_clock::now();
//             auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
//             std::cout << ms_int.count() << "ms\n";


// }






#endif