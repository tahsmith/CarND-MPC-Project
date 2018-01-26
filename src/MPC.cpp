#include <algorithm>
#include <chrono>


#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

using CppAD::AD;

// Config for the MPC system
const size_t N = 20;
const double dt = 0.05;
const size_t n_state = 4 * N;
const size_t n_control = 2 * (N - 1);
const size_t n_vars = n_state + n_control;
const size_t n_constraints = 4 * (N - 1);
const double target_v = 120;
const double max_a = 100;
const double max_steer = 25.0 * M_PI / 180.0;
const double steering_aggression = 1.0 * max_steer;
const double acceleration_aggression = 100;


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

// Utility functions to keep sanity
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

    /// Calculate the cost associated with a trajectory that does not follow the waypoint curve.
    AD<double> crossTrackError(const ADvector& vars) {
        using namespace state_indices;
        AD<double> cte = 0.0;
        for (size_t i = 0; i < N; ++i) {
            auto x = vars[x_i(i)];
            auto y = vars[y_i(i)];
            auto err = polyEval(coeffs, x) - y;
            cte += CppAD::pow(err, 2);
        }

        return cte;
    }

    /// Calculate the cost associated with a trajectory that is not tangential to the of the waypoint curve
    AD<double> tangentialError(const ADvector& vars) {
        using namespace state_indices;
        AD<double> tangentialError = 0.0;
        for (size_t i = 1; i < N; ++i) {
            auto psi = vars[psi_i(i)];
            auto x0 = vars[x_i(i - 1)];
            auto x1 = vars[x_i(i)];
            auto y0 = polyEval(coeffs, x0);
            auto y1 = polyEval(coeffs, x1);
            if (CppAD::abs(x0 - x1) > 0.001) {
                tangentialError += CppAD::pow(psi - atan2(y1 - y0, x1 - x0), 2);
            }
        }
        return tangentialError;
    }

    /// Calculate the cost associated with a trajectory that does not keep up with the speed set point (target_v)
    AD<double> speedPenalty(const ADvector& vars) {
        using namespace state_indices;
        AD<double> speedError = 0.0;
        for (size_t i = 0; i < N; ++i) {
            speedError += CppAD::pow(vars[v_i(i)] - target_v, 2);
        }
        return speedError;
    }

    /// Calculate the cost associated a control trajectory that changes rapidly. This attempt to smooth the controls by
    /// minimising the difference between successive control values.
    AD<double> smooth(const ADvector& vars) {
        using namespace state_indices;
        AD<double> smooth_cost = 0;
        for (size_t i = 0; i < N - 2; i++) {
            smooth_cost += CppAD::pow((vars[delta_i(i + 1)] - vars[delta_i(i)]) / steering_aggression, 2);
            smooth_cost += CppAD::pow((vars[a_i(i + 1)] - vars[a_i(i)] ) / acceleration_aggression, 2);
        }
        return smooth_cost;
    }

    /// Calculate the cost associated with a trajectory that corners too hard. The cost per time step here is
    /// proportional to the centripetal force.
    AD<double> corneringPenalty(const ADvector& vars) {
        using namespace state_indices;
        AD<double> cornering_penalty = 0;
        for (size_t i = 0; i < N - 2; i++) {
            cornering_penalty += CppAD::pow(vars[v_i(i)] * vars[v_i(i)] * vars[delta_i(i)]/ max_steer, 2);
        }
        return cornering_penalty;
    }

    /// Set up the constrains associated with the process model, i.e., all pairs successive time points in the state
    /// space must obey the dynamic model. Here that model is the bicycle model.
    void applyProcessModelConstraint(const ADvector& vars, size_t i, ADvector& fg) {
        using namespace state_indices;
        auto x0 = vars[x_i(i)];
        auto y0 = vars[y_i(i)];
        auto psi0 = vars[psi_i(i)];
        auto v0 = vars[v_i(i)];
        auto delta = vars[delta_i(i)];
        auto a = vars[a_i(i)];

        auto r = v0 * dt;
        auto x1_pred = x0 + r * CppAD::cos(psi0);
        auto y1_pred = y0 + r * CppAD::sin(psi0);
        auto psi1_pred = psi0 + r * delta / Lf;
        auto v1_pred = v0 + a * dt;

        auto x1 = vars[x_i(i + 1)];
        auto y1 = vars[y_i(i + 1)];
        auto psi1 = vars[psi_i(i + 1)];
        auto v1 = vars[v_i(i + 1)];

        auto x_err = x1_pred - x1;
        auto y_err = y1_pred - y1;
        auto psi_err = psi1_pred - psi1;
        auto v_err = v1_pred - v1;
        fg[i + 0 * (N - 1) + 1] = x_err;
        fg[i + 1 * (N - 1) + 1] = y_err;
        fg[i + 2 * (N - 1) + 1] = psi_err;
        fg[i + 3 * (N - 1) + 1] = v_err;
    }

    /// Apply all of our costs constraints.
    void operator()(ADvector& fg, const ADvector& vars)
    {
        fg[0] = 0;
        fg[0] += 100 * crossTrackError(vars);
        fg[0] += 100 * tangentialError(vars);
        fg[0] += 100 * smooth(vars);
        fg[0] += 2 * speedPenalty(vars);
        fg[0] += 0.005 * corneringPenalty(vars);
        for(size_t i = 0; i < N - 1; ++i) {
            applyProcessModelConstraint(vars, i, fg);
        }

        std::cout << fg << std::endl;
    }
};

auto MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) -> Solution
{
    using namespace state_indices;
    bool ok = true;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    Dvector vars(n_vars);
    for (int i = 0; i < n_vars; i++)
    {
        vars[i] = 0;
    }

    // Set the initial state for the solver. We are going to constrain these points to be equal to these values later.
    vars[x_i(0)] = state(0);
    vars[y_i(0)] = state(1);
    vars[psi_i(0)] = state(2);
    vars[v_i(0)] = state(3);
    vars[a_i(0)] = state(4);
    vars[delta_i(0)] = state(5);

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);

    // N.B. there is one less control point that state points.
    for (size_t i = 0; i < N - 1; ++i) {
        vars_lowerbound[a_i(i)] = -max_a;
        vars_upperbound[a_i(i)] = max_a;
        vars_lowerbound[delta_i(i)] = -max_steer;
        vars_upperbound[delta_i(i)] = max_steer;
    }

    // Set a big number for all the state variables. We are not interested in constraining the problem this way.
    for (size_t i = 1; i < N; ++i) {
        vars_lowerbound[x_i(i)] = -10000;
        vars_upperbound[x_i(i)] = 10000;
        vars_lowerbound[y_i(i)] = -10000;
        vars_upperbound[y_i(i)] = 10000;
        vars_lowerbound[psi_i(i)] = -10000;
        vars_upperbound[psi_i(i)] = 10000;
        vars_lowerbound[v_i(i)] = -10000;
        vars_upperbound[v_i(i)] = 10000;
    }

    // Constrain the initial state + controls to be the ones we have measured and passed in via `state`.
    vars_lowerbound[x_i(0)] = state(0);
    vars_upperbound[x_i(0)] = state(0);
    vars_lowerbound[y_i(0)] = state(1);
    vars_upperbound[y_i(0)] = state(1);
    vars_lowerbound[psi_i(0)] = state(2);
    vars_upperbound[psi_i(0)] = state(2);
    vars_lowerbound[v_i(0)] = state(3);
    vars_upperbound[v_i(0)] = state(3);
    vars_lowerbound[a_i(0)] = state(4);
    vars_upperbound[a_i(0)] = state(4);
    vars_lowerbound[delta_i(0)] = state(5);
    vars_upperbound[delta_i(0)] = state(5);

    // All of our constrains are strict equality constraints, so set upper and lower bound to 0
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (int i = 0; i < n_constraints; i++)
    {
        constraints_lowerbound[i] = 0.0;
        constraints_upperbound[i] = 0.0;
    }

    // object that computes objective and constraints
    FG_eval fg_eval(std::move(coeffs));

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
    switch (solution.status) {
        case CppAD::ipopt::solve_result<Dvector>::not_defined:
            std::cout << "not_defined\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::success:
            std::cout << "success\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::feasible_point_found:
            std::cout << "feasible_point_found\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::unknown:
            std::cout << "unknown\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::maxiter_exceeded:
            std::cout << "maxiter_exceeded\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::stop_at_tiny_step:
            std::cout << "stop_at_tiny_step\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::stop_at_acceptable_point:
            std::cout << "stop_at_acceptable_point\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::local_infeasibility:
            std::cout << "local_infeasibility\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::user_requested_stop:
            std::cout << "user_requested_stop\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::diverging_iterates:
            std::cout << "diverging_iterates\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::restoration_failure:
            std::cout << "restoration_failure\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::error_in_step_computation:
            std::cout << "error_in_step_computation\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::invalid_number_detected:
            std::cout << "invalid_number_detected\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::too_few_degrees_of_freedom:
            std::cout << "too_few_degrees_of_freedom\n";
            break;
        case CppAD::ipopt::solve_result<Dvector>::internal_error:
            std::cout << "internal_error\n";
            break;
    }
    if (!ok)
    {
        exit(1);
    }

    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    std::vector<double> x;
    std::copy(solution.x.data() + x_i(0), solution.x.data() + x_i(N), std::back_inserter(x));
    std::vector<double> y;
    std::copy(solution.x.data() + y_i(0), solution.x.data() + y_i(N), std::back_inserter(y));
    std::vector<double> psi;
    std::copy(solution.x.data() + psi_i(0), solution.x.data() + psi_i(N), std::back_inserter(psi));
    std::vector<double> v;
    std::copy(solution.x.data() + v_i(0), solution.x.data() + v_i(N), std::back_inserter(v));
    std::vector<double> a;
    std::copy(solution.x.data() + a_i(0), solution.x.data() + a_i(N - 1), std::back_inserter(a));
    std::vector<double> delta;
    std::copy(solution.x.data() + delta_i(0), solution.x.data() + delta_i(N - 1), std::back_inserter(delta));

    return {
            x,
            y,
            psi,
            v,
            a,
            delta
    };
}

MPC::MPC(){

}

PID::PID():
        p_error(0),
        i_error(0),
        d_error(0)
{}

void PID::Init(double Kp, double Ki, double Kd) {
    this->Kp = Kp;
    this->Ki = Ki;
    this->Kd = Kd;
}

void PID::UpdateError(double cte) {
    d_error = cte - p_error;
    p_error = cte;
    i_error += cte;
}

double PID::TotalError() {
    return -(Kp * p_error + Kd * d_error + Ki * i_error);
}

AccelerationController::AccelerationController(double Kp, double Ki, double Kd) : v(0), a(0), t(0.0)  {
    pid.Init(Kp, Ki, Kd);
}

void AccelerationController::Update(double v) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    if (t == 0.0) {
        this->v = v;
        t = duration_cast<milliseconds>(high_resolution_clock::now().time_since_epoch()).count() * 1e-3;
    }
    else
    {
        double now = duration_cast<milliseconds>(high_resolution_clock::now().time_since_epoch()).count() * 1e-3;
        double dt = now - t;
        a = (v - this->v) / dt;
        this->v = v;
        t = now;
    }
}

double AccelerationController::Throttle(double accelerationSetPoint) {
    pid.UpdateError(a - accelerationSetPoint);
    return pid.TotalError();
}
