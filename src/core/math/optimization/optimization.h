#pragma once
#include <iostream>
#include <vector>
#include <chrono>
#include <set>
#include <iterator>
#include "../../../../src/errors/core/math.h"

struct OptimizerResult
{
    const std::vector<double> x_; 
    const double function_value_; 
    const int number_iterations_; 
    const double tolerance_threshold_; 
    const int maximum_iterations_;
    const double time_taken_; 
    OptimizerResult(
        const std::vector<double> x, const double function_value, const int number_iterations,
        const double tolerance_threshold, const int maximum_iterations, const double time_taken): 
        x_(x), function_value_(function_value), number_iterations_(abs(number_iterations)),
        tolerance_threshold_(abs(tolerance_threshold)), maximum_iterations_(abs(maximum_iterations)),
        time_taken_(time_taken){};
    ~OptimizerResult(){}; 
        
};

class Optimizer
{
    public: 
        virtual std::shared_ptr<OptimizerResult> optimize() const = 0;
        Optimizer(){};
        virtual ~Optimizer(){};
};

class NewtonRaphson : public Optimizer
{
    public: 
        void set_tolerance_rate(double tol) {tol_ = abs(tol);}; 
        void set_max_iterations(int max_iter){max_iter_ = abs(max_iter);};
        std::shared_ptr<OptimizerResult> optimize() const override;
        NewtonRaphson(double x0,std::function<double(double)> f, std::function<double(double)> df):
        x0_(x0), f_(f), df_(df)
        {
            set_max_iterations(100); 
            set_tolerance_rate(1e-6);
        };
        ~NewtonRaphson(){};
    private: 
        const double x0_;
        const std::function<double(double)> f_;
        const std::function<double(double)> df_;
        int max_iter_;
        double tol_; 
};

struct NelderMeadVertex
{
    std::vector<double> points_; 
    double value_; 
    NelderMeadVertex(const std::vector<double> points, const double value): points_(points), value_(value){}; 
    ~NelderMeadVertex(){};

};

struct NelderMeadSimplex
{
    std::vector<std::shared_ptr<NelderMeadVertex>> vertices_; 
    void sort_vertices();
    void set_vertex(std::shared_ptr<NelderMeadVertex> vertex, int i){vertices_[i] = vertex; sort_vertices();};
    std::shared_ptr<NelderMeadVertex> get_best_vertex() const {return vertices_[0];};
    std::shared_ptr<NelderMeadVertex> get_worst_vertex() const {return vertices_.back();};
    std::shared_ptr<NelderMeadVertex> get_second_worst_vertex() const {int n = vertices_.size(); return vertices_[n-2];};
    NelderMeadSimplex(std::vector<std::shared_ptr<NelderMeadVertex>> vertices): vertices_(vertices) {sort_vertices();};
    ~NelderMeadSimplex(){};
};

class NelderMead : public Optimizer
{
    public: 
        void set_initial_simplex_method(std::string init_method);
        void set_max_iterations(int max_iter){max_iter_ = abs(max_iter);};
        void set_perturbation_parameter(double epsilon){epsilon_ = abs(epsilon);};
        void set_reflection_parameter(double alpha){alpha_ = abs(alpha);};
        void set_expansion_parameter(double beta) {if (beta>1) beta_ = beta;};
        void set_contraction_parameter(double gamma) {if (gamma<1 and gamma>0) gamma_ = gamma;};
        void set_shrink_parameter(double delta) {if (delta<1 and delta>0) delta_ = delta;};
        void set_tolerance_rate(double tol) {tol_ = abs(tol);};
        std::shared_ptr<OptimizerResult> optimize() const override;
        NelderMead(std::vector<double> x0,std::function<double(std::vector<double>)> f): 
        x0_(x0), n(x0.size()), f_(f) 
        {
            set_perturbation_parameter(0.05); 
            set_reflection_parameter(1.0); 
            set_expansion_parameter(2.0); 
            set_contraction_parameter(0.5); 
            set_shrink_parameter(0.5); 
            set_max_iterations(100); 
            set_tolerance_rate(1e-6);
            set_initial_simplex_method("basic");
        }
        ~NelderMead(){};
    private: 
        const std::vector<double> x0_;
        const int n; 
        const std::function<double(std::vector<double>)> f_;
        int max_iter_;
        double tol_; 
        double epsilon_; 
        double alpha_; 
        double beta_; 
        double gamma_; 
        double delta_;
        std::string init_simplex_method;
        double get_symalpha() const {return epsilon_/(2*sqrt(n));};
        std::shared_ptr<NelderMeadVertex> get_vertex(std::vector<double> x_input) const {return std::make_shared<NelderMeadVertex>(x_input, f_(x_input));};
        std::shared_ptr<NelderMeadSimplex> get_initial_simplex() const; 
        std::shared_ptr<NelderMeadVertex> get_centroid(std::shared_ptr<NelderMeadSimplex> simplex) const;
        std::shared_ptr<NelderMeadVertex> get_reflection(std::shared_ptr<NelderMeadVertex> centroid, std::shared_ptr<NelderMeadVertex> worst) const; 
        std::shared_ptr<NelderMeadVertex> get_expansion(std::shared_ptr<NelderMeadVertex> centroid, std::shared_ptr<NelderMeadVertex> reflection) const; 
        std::shared_ptr<NelderMeadVertex> get_contraction(std::shared_ptr<NelderMeadVertex> centroid, std::shared_ptr<NelderMeadVertex> worst) const; 
        std::shared_ptr<NelderMeadSimplex> get_shrink(std::shared_ptr<NelderMeadSimplex> simplex) const; 
        bool check_convergence(std::shared_ptr<NelderMeadSimplex> simplex) const;
};

