#include "optimization.h"

std::shared_ptr<OptimizerResult> NewtonRaphson::optimize() const
{
    auto start = std::chrono::high_resolution_clock::now();
    double x = x0_;
    double fx; 
    double dfx; 
    double x_step = 0;
    std::string error_message = "None";
    int iter = 0;
    for (int i = 1; i <= max_iter_; ++i) {

        iter += 1;
        fx = f_(x);
        dfx = df_(x);
        if (std::abs(fx) < tol_){
            break;
        }
        
        if (std::abs(dfx) < 1e-12) {error_message = MathError(MathErrorCode::NEWTON_RAPHSON_DERIVATIVE_CLOSE_TO_ZERO).what();break;}
        
        double x_new = x - fx / dfx;
        x_step = std::abs(x_new - x);
        if (x_step < tol_) {
            x = x_new;
            break;
        }

        if (i==max_iter_){break;}
        
        x = x_new;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::vector<double> x_out = {x};
    return std::make_shared<OptimizerResult>(x_out,fx,iter,tol_,max_iter_,elapsed.count());
};

void NelderMeadSimplex::sort_vertices()
{
    std::sort(vertices_.begin(), vertices_.end(), 
    [](const std::shared_ptr<NelderMeadVertex>& a, const std::shared_ptr<NelderMeadVertex>& b) {
        return a->value_ < b->value_;  // Compare values directly
    });
};

void NelderMead::set_initial_simplex_method(std::string init_method)
{
    if (init_method == "basic"){init_simplex_method = init_method;}
    else if (init_method == "scaled"){init_simplex_method = init_method;}
    else if (init_method == "symmetric"){init_simplex_method = init_method;}
    
};

std::shared_ptr<NelderMeadSimplex> NelderMead::get_initial_simplex() const
{
    std::vector<std::shared_ptr<NelderMeadVertex>> simplex;
    double a = get_symalpha();
    simplex.push_back(get_vertex(x0_));
    for (int i = 0; i < n; ++i) {
        std::vector<double> x = x0_; 
        if (init_simplex_method=="basic"){
            x[i] += epsilon_;
        }
        else if (init_simplex_method=="scaled"){
            x[i] += epsilon_ * (1 + fabs(x0_[i]));
        }
        else if (init_simplex_method=="symmetric"){
            for (int j = 0; j < n; ++j) {
                x[j] += (i == j) ? a : (-a / (n - 1));
            };
        }
        simplex.push_back(get_vertex(x));
    };
    return std::make_shared<NelderMeadSimplex>(simplex);
}

std::shared_ptr<NelderMeadVertex> NelderMead::get_centroid(std::shared_ptr<NelderMeadSimplex> simplex) const
{
    std::vector<std::shared_ptr<NelderMeadVertex>> simplex_vec = simplex->vertices_;
    std::vector<double> x_c(n, 0); 
    for (int j = 0; j < n; ++j){
        for (int i = 0; i < n; ++i){
            x_c[i] +=  simplex_vec[j]->points_[i]/n;
        }
    }
    return get_vertex(x_c);
}

std::shared_ptr<NelderMeadVertex> NelderMead::get_reflection(std::shared_ptr<NelderMeadVertex> centroid, std::shared_ptr<NelderMeadVertex> worst) const
{

    std::vector<double> x(n); 
    for (int j = 0; j < n; ++j){
        double xw = worst->points_[j];
        double xc = centroid->points_[j]; 
        x[j] = xc + alpha_*(xc-xw);
    }
    return get_vertex(x);
};

std::shared_ptr<NelderMeadVertex> NelderMead::get_expansion(std::shared_ptr<NelderMeadVertex> centroid, std::shared_ptr<NelderMeadVertex> reflection) const
{
    std::vector<double> x(n); 
    for (int j = 0; j < n; ++j){
        double xr = reflection->points_[j];
        double xc = centroid->points_[j]; 
        x[j] = xc + beta_*(xr-xc);
    }
    return get_vertex(x);
};

std::shared_ptr<NelderMeadVertex> NelderMead::get_contraction(std::shared_ptr<NelderMeadVertex> centroid, std::shared_ptr<NelderMeadVertex> worst) const
{
    std::vector<double> x(n); 
    for (int j = 0; j < n; ++j){
        double xw = worst->points_[j];
        double xc = centroid->points_[j]; 
        x[j] = xc + gamma_*(xw-xc);
    }
    return get_vertex(x);
};

std::shared_ptr<NelderMeadSimplex> NelderMead::get_shrink(std::shared_ptr<NelderMeadSimplex> simplex) const
{
    std::vector<std::shared_ptr<NelderMeadVertex>> simplex_out;
    std::vector<std::shared_ptr<NelderMeadVertex>> simplex_in = simplex->vertices_;
    simplex_out.push_back(simplex_in[0]); 
    for (int j = 1; j < n; ++j){
        std::vector<double> x(n); 
        for (int i = 0; i < n; ++i){
            x[i] = simplex_in[0]->points_[i] + delta_*(simplex_in[j]->points_[i] - simplex_in[0]->points_[i]);
        }
        simplex_out.push_back(get_vertex(x)); 
    } 
    return std::make_shared<NelderMeadSimplex>(simplex_out);
};

bool NelderMead::check_convergence(std::shared_ptr<NelderMeadSimplex> simplex) const
{
    double max_distance = 0.0;
    double max_value = 0.0;
    for (int i = 1; i < n; i++) {
        double distance = 0.0;
        for (int j = 0; j < n; j++) {
            distance += std::pow(simplex->vertices_[i]->points_[j] - simplex->vertices_[0]->points_[j], 2);
        }
        max_distance = std::max(max_distance, std::sqrt(distance));
        max_value = std::max(max_value, std::abs(simplex->vertices_[i]->value_ - simplex->vertices_[0]->value_));
    }
    return (max_distance < tol_ || max_value < tol_);
};

std::shared_ptr<OptimizerResult> NelderMead::optimize() const
{
    auto start = std::chrono::high_resolution_clock::now();
    std::shared_ptr<NelderMeadSimplex> simplex = get_initial_simplex(); 
    int i = 1;
    for (int iter = 1; iter <= max_iter_; iter++) {
        if (check_convergence(simplex)){break;}
        std::shared_ptr<NelderMeadVertex> centroid = get_centroid(simplex); 
        std::shared_ptr<NelderMeadVertex> worst = simplex->get_worst_vertex(); 
        std::shared_ptr<NelderMeadVertex> reflection = get_reflection(centroid,worst);
        std::shared_ptr<NelderMeadVertex> best = simplex->get_best_vertex(); 
        std::shared_ptr<NelderMeadVertex> second_worst = simplex->get_second_worst_vertex(); 

        if (reflection->value_<best->value_){
            std::shared_ptr<NelderMeadVertex> expansion = get_expansion(centroid, reflection); 
            if (expansion->value_<reflection->value_){
                simplex->set_vertex(expansion, simplex->vertices_.size()-1);
            }else{
                simplex->set_vertex(reflection, simplex->vertices_.size()-1);
            }
        }
        else if (reflection->value_>=best->value_ and reflection->value_<second_worst->value_){
            simplex->set_vertex(reflection, simplex->vertices_.size()-1);
        }
        else {
            
            if (reflection->value_<worst->value_){
                std::shared_ptr<NelderMeadVertex> contraction = get_contraction(centroid, reflection);
                if (contraction->value_<reflection->value_){
                    simplex->set_vertex(contraction, simplex->vertices_.size()-1);
                }else{
                    simplex = get_shrink(simplex);
                }
            }else{
                std::shared_ptr<NelderMeadVertex> contraction = get_contraction(centroid, worst);
                if (contraction->value_<worst->value_){
                    simplex->set_vertex(contraction, simplex->vertices_.size()-1);
                }else{
                    simplex = get_shrink(simplex);
                }
            }
        }
        if (i==max_iter_){break;}
        
        i++;

    };

    std::shared_ptr<NelderMeadVertex> final_vertex = simplex->vertices_[0]; 
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return std::make_shared<OptimizerResult>(final_vertex->points_,final_vertex->value_,i,tol_,max_iter_,elapsed.count());
    
};