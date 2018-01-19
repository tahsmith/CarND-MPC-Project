#include <algorithm>

#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

using CppAD::AD;

// TODO: Set the timestep length and duration
const size_t N = 10;
const double dt = 0.1;
const size_t n_state = 4 * N;
const size_t n_control = 2 * (N - 1);
const size_t n_vars = n_state + n_control;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

AD<double> polyEval(Eigen::VectorXd coeffs, AD<double> x)
{
    AD<double> result = 0.0;
    for (int i = 0; i < coeffs.size(); i++)
    {
        result += coeffs[i] * CppAD::pow(x, i);
    }
    return result;
}


Eigen::VectorXd globalKinematic(Eigen::VectorXd state,
                                Eigen::VectorXd actuators, double dt)
{
    Eigen::VectorXd next_state(state.size());

    // NOTE: state is [x, y, psi, v]
    // NOTE: actuators is [delta, a]
    double r = state(3) * dt;
    next_state(0) = state(0) + r * cos(state(2));
    next_state(1) = state(1) + r * sin(state(2));
    next_state(2) = state(2) + r * actuators(0) / Lf;
    next_state(3) = state(3) + actuators(1) * dt;
    return next_state;
}

namespace state_indices
{
    size_t x_i(size_t i)
    {
        return i;
    }

    size_t y_i(size_t i)
    {
        return i + x_i(N);
    }

    size_t psi_i(size_t i)
    {
        return i + y_i(N);
    }

    size_t v_i(size_t i)
    {
        return i + psi_i(N);
    }

    size_t delta_i(size_t i)
    {
        return i + v_i(N);
    }

    size_t a_i(size_t i)
    {
        return i + delta_i(N - 1);
    }
}

class FG_eval
{
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;

    FG_eval(Eigen::VectorXd coeffs)
    { this->coeffs = coeffs; }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    AD<double> crossTrackError(const ADvector& vars) {
        using namespace state_indices;
        AD<double> cte = 0.0;
        for (size_t i = 0; i < N; ++i) {
            auto x = vars[x_i(i)];
            auto err = polyEval(coeffs, x) - x;
            cte += CppAD::pow(err, 2);
        }

        return cte;
    }

    AD<double> distance(const ADvector& vars) {
        using namespace state_indices;
        AD<double> distance = 0.0;
        for (size_t i = 1; i < N; ++i) {
            auto x0 = vars[x_i(i - 1)];
            auto y0 = vars[y_i(i - 1)];
            auto x1 = vars[x_i(i)];
            auto y1 = vars[y_i(i)];
            auto dx = x1 - x0;
            auto dy = y1 - y0;
            distance += CppAD::sqrt(dx * dx + dy * dy);
        }

        return distance;
    }

    void applyProcessModelConstraint(const ADvector& vars, size_t i, ADvector& fg) {
        using namespace state_indices;
        auto x0 = vars[x_i(i - 1)];
        auto y0 = vars[y_i(i - 1)];
        auto psi0 = vars[psi_i(i - 1)];
        auto v0 = vars[v_i(i - 1)];
        auto delta = vars[delta_i(i - 1)];
        auto a = vars[a_i(i - 1)];

        auto r = v0 * dt;
        auto x1_pred = x0 + r * CppAD::cos(psi0);
        auto y1_pred = y0 + r * CppAD::sin(psi0);
        auto psi1_pred = psi0 + r * delta / Lf;
        auto v1_pred = v0 + a * dt;

        auto x1 = vars[x_i(i)];
        auto y1 = vars[y_i(i)];
        auto psi1 = vars[psi_i(i)];
        auto v1 = vars[v_i(i)];

        auto x_err = x1_pred - x1;
        auto y_err = y1_pred - y1;
        auto psi_err = psi1_pred - psi1;
        auto v_err = v1_pred - v1;
        fg[i + 1] = x_err * x_err;
        fg[i + N + 1] = y_err * y_err;
        fg[i + 2 * N + 1] = psi_err * psi_err;
        fg[i + 3 * N + 1] = v_err * v_err;
    }

    void operator()(ADvector& fg, const ADvector& vars)
    {
        // TODO: implement MPC
        // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
        // NOTE: You'll probably go back and forth between this function and
        // the Solver function below.

        fg[0] = crossTrackError(vars) - distance(vars);
        for(size_t i = 1; i < N; ++i) {
            applyProcessModelConstraint(vars, i, fg);
        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC()
{}

MPC::~MPC()
{}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs)
{
    using namespace state_indices;
    bool ok = true;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    // TODO: Set the number of model variables (includes both states and inputs).
    // For example: If the state is a 4 element vector, the actuators is a 2
    // element vector and there are 10 timesteps. The number of variables is:
    //
    // 4 * 10 + 2 * 9
    // TODO: Set the number of constraints
    size_t n_constraints = 4 * N;

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++)
    {
        vars[i] = 0;
    }
    vars[x_i(0)] = state(0);
    vars[y_i(0)] = state(1);
    vars[psi_i(0)] = state(2);
    vars[v_i(0)] = state(3);

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    // TODO: Set lower and upper limits for variables.
    for (size_t i = 0; i < N - 1; ++i) {
        vars_lowerbound[a_i(i)] = -1.0;
        vars_upperbound[a_i(i)] = 1.0;
    }

    // Lower and upper limits for the constraints
    // Should be 0 besides initial state.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (int i = 0; i < n_constraints; i++)
    {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[x_i(0)] = state(0);
    constraints_lowerbound[y_i(0)] = state(1);
    constraints_lowerbound[psi_i(0)] = state(2);
    constraints_lowerbound[v_i(0)] = state(3);

    constraints_upperbound[x_i(0)] = state(0);
    constraints_upperbound[y_i(0)] = state(1);
    constraints_upperbound[psi_i(0)] = state(2);
    constraints_upperbound[v_i(0)] = state(3);

    // object that computes objective and constraints
    FG_eval fg_eval(move(coeffs));

    //
    // NOTE: You don't have to worry about these options
    //
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
        options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
        constraints_upperbound, fg_eval, solution);

    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    // TODO: Return the first actuator values. The variables can be accessed with
    // `solution.x[i]`.
    //
    // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
    // creates a 2 element double vector.
    vector<double> output;
    std::copy(solution.x.data(), solution.x.data() + solution.x.size(), std::back_inserter(output));
    return output;
}
