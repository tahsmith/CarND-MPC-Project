#include <math.h>
#include <uWS/uWS.h>
#include <thread>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
using std::string;
using std::cout;
using std::endl;
using std::vector;
namespace this_thread = std::this_thread;
namespace chrono = std::chrono;


// For converting back and forth between radians and degrees.
constexpr double pi()
{ return M_PI; }

double deg2rad(double x)
{ return x * pi() / 180; }

double rad2deg(double x)
{ return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s)
{
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos)
    {
        return "";
    }
    else if (b1 != string::npos && b2 != string::npos)
    {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x)
{
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++)
    {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order)
{
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++)
    {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++)
    {
        for (int i = 0; i < order; i++)
        {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}

void transformGlobalToLocal(vector<double>& ptsx, vector<double>& ptsy, double px, double py, double psi)
{
    for(int i = 0; i < ptsx.size(); i++) {
        double shift_x = ptsx[i] - px;
        double shift_y = ptsy[i] - py;

        ptsx[i] = shift_x * cos(-psi) - shift_y*sin(-psi);
        ptsy[i] = shift_x * sin(-psi) + shift_y*cos(-psi);
    }
}

void transformLocalToGlobal(vector<double>& ptsx, vector<double>& ptsy, double px, double py, double psi)
{
    for(int i = 0; i < ptsx.size(); i++) {
        double x = ptsx[i] * cos(psi) - ptsy[i] * sin(psi);
        double y = ptsy[i] * sin(psi) + ptsy[i] * cos(psi);
        ptsx[i] = px + x;
        ptsy[i] = py + y;
    }
}

int main()
{
    uWS::Hub h;

    h.onMessage(
        [](uWS::WebSocket<uWS::SERVER> ws, char* data, size_t length,
               uWS::OpCode opCode) {
            // "42" at the start of the message means there's a websocket message event.
            // The 4 signifies a websocket message
            // The 2 signifies a websocket event
            string sdata = string(data).substr(0, length);
            cout << sdata << endl;
            if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2')
            {
                string s = hasData(sdata);
                if (s != "")
                {
                    auto j = json::parse(s);
                    string event = j[0].get<string>();
                    if (event == "telemetry")
                    {
                        // j[1] is the data JSON object
                        vector<double> ptsx = j[1]["ptsx"];
                        vector<double> ptsy = j[1]["ptsy"];
                        double px = j[1]["x"];
                        double py = j[1]["y"];
                        double psi = j[1]["psi"];
                        double v = j[1]["speed"];
                        double a = j[1]["throttle"];
                        double delta = j[1]["steering_angle"];
                        delta *= -deg2rad(25);

                        auto waypoints_x = ptsx;
                        auto waypoints_y = ptsy;
                        transformGlobalToLocal(waypoints_x, waypoints_y, px, py, psi);

                        /*
                        * TODO: Calculate steering angle and throttle using MPC.
                        *
                        * Both are in between [-1, 1].
                        *
                        */

                        Eigen::VectorXd state(6);
                        state << 0, 0, 0, v, a, delta;

                        auto coeffs = polyfit(
                            Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(waypoints_x.data(), waypoints_x.size()),
                            Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(waypoints_y.data(), waypoints_y.size()),
                            4);
                        auto soln = MPC::Solve(state, coeffs);

                        double future_i = 2;
                        double future_x = soln.x[future_i];
                        double future_y = soln.y[future_i];
                        double future_psi = soln.psi[future_i];
                        double steer_value = -soln.delta[future_i] / deg2rad(25);
                        double throttle_value = soln.a[future_i];

                        json msgJson;
                        // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
                        // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
                        msgJson["steering_angle"] = steer_value;
                        msgJson["throttle"] = throttle_value;


                        //Display the MPC predicted trajectory
                        vector<double> mpc_x_vals = soln.x;
                        vector<double> mpc_y_vals = soln.y;
                        transformGlobalToLocal(mpc_x_vals, mpc_y_vals, future_x, future_y, future_psi);

                        //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                        // the points in the simulator are connected by a Green line

                        msgJson["mpc_x"] = mpc_x_vals;
                        msgJson["mpc_y"] = mpc_y_vals;

                        //Display the waypoints/reference line
                        vector<double> next_x_vals = waypoints_x;
                        vector<double> next_y_vals = waypoints_y;
                        for (size_t i = 0; i < waypoints_x.size(); ++i) {
                            next_y_vals.push_back(polyeval(coeffs, waypoints_x[i]));
                        }
                        transformGlobalToLocal(next_x_vals, next_y_vals, future_x, future_y, future_psi);

                        //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                        // the points in the simulator are connected by a Yellow line

                        msgJson["next_x"] = next_x_vals;
                        msgJson["next_y"] = next_y_vals;


                        auto msg = "42[\"steer\"," + msgJson.dump() + "]";
                        std::cout << msg << std::endl;
                        // Latency
                        // The purpose is to mimic real driving conditions where
                        // the car does actuate the commands instantly.
                        //
                        // Feel free to play around with this value but should be to drive
                        // around the track with 100ms latency.
                        //
                        // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
                        // SUBMITTING.
                        this_thread::sleep_for(chrono::milliseconds(100));
                        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                    }
                }
                else
                {
                    // Manual driving
                    std::string msg = "42[\"manual\",{}]";
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            }
        });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse* res, uWS::HttpRequest req, char* data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1)
        {
            res->end(s.data(), s.length());
        }
        else
        {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char* message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port))
    {
        std::cout << "Listening to port " << port << std::endl;
    }
    else
    {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}
