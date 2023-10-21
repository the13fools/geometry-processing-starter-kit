#ifndef UTILSMISC
#define UTILSMISC
#include <Eigen/Core>

#include <Eigen/Dense>

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
        
        for(int i = 0; i < 5; i++)
        {
            // std::cout << v.transpose() << std::endl;

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
      
        }

        std::cout << v.transpose() << std::endl;
        if(v.norm() < 1e-4)
        {
            v = Eigen::Vector2d::Random();
        }
}


#endif