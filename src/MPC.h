#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

class MPC {
public:
    double v, a;
    double t;

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
    // Return the first actuatotions.
    Solution Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

    void Update(double v);
};

#endif /* MPC_H */
