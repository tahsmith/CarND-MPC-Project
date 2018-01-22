#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

class PID {
public:
    /*
    * Errors
    */
    double p_error;
    double i_error;
    double d_error;

    /*
    * Coefficients
    */
    double Kp;
    double Ki;
    double Kd;

    /*
    * Constructor
    */
    PID();

    /*
    * Initialize PID.
    */
    void Init(double Kp, double Ki, double Kd);

    /*
    * Update the PID error variables given cross track error.
    */
    void UpdateError(double cte);

    /*
    * Calculate the total PID error.
    */
    double TotalError();
};

/// Measures acceleration from sequential velocity measurements.
/// Produces a
class AccelerationController {
    PID pid;
public:
    double v, a;
    double t;


    AccelerationController(double Kp, double Ki, double Kd);

    void Update(double v);

    double Throttle(double accelerationSetPoint);
};

class MPC {
public:
    MPC();

    struct Solution {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> psi;
        std::vector<double> v;
        std::vector<double> a;
        std::vector<double> delta;
    };

    // Solve the model given an initial state and polynomial coefficients.
    Solution Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
