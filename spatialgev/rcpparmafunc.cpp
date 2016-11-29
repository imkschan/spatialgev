// # include <Rcpp.h>
# include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    
    for (int i=0; i < n; i++) {
        arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
        out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
    }  
      
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat Sigma) {
   int ncols = Sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}


// function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE) 
// {
//     if (min(scale) <= 0) 
//         stop("invalid scale")
//     if (length(shape) != 1) 
//         stop("invalid shape")
//     q <- (q - loc)/scale
//     if (shape == 0) 
//         p <- exp(-exp(-q))
//     else p <- exp(-pmax(1 + shape * q, 0)^(-1/shape))
//     if (!lower.tail) 
//         p <- 1 - p
//     p
// }

// [[Rcpp::export]]
double pgevcpp(double q, double loc = 0, double sc = 1, double sh = 0){
    q = (q - loc) / sc;
    if (sh == 0){
        return exp(-exp(-q));
    } else {
        double tmp = 1 + sh * q;
        if (tmp < 0) tmp = 0;
        return exp(-pow(tmp, -1 / sh));
    }
}

// [[Rcpp::export]]
NumericMatrix pgev_mat(NumericMatrix data, NumericMatrix loc, NumericVector sc, NumericVector sh){
    int data_nrow = data.nrow();
    int data_ncol = data.ncol();
    NumericMatrix mat(data_nrow,data_ncol);
    for (int i = 0; i < data_nrow; i++){
        for (int k = 0; k < data_ncol; k++)
            mat(i, k) = pgevcpp(data(i, k), loc(i, k), sc(k), sh(k));
    }
    return mat;
}

// [[Rcpp::export]]
double dgevcpp(double x, double loc = 0, double sc = 1, double sh = 0, bool takelog = false){
    double tmp = 0.0;
    x = (x - loc) / sc;
    double xx = 1 + sh * x;
    if (sh == 0){
        tmp = log(1/sc) - x - exp(-x);
    } else {
        if (xx > 0) {
            tmp = log(1/sc) - pow(xx,-1/sh) - (1/sh + 1) * log(xx);
        } else {
            tmp = -INFINITY;
        }
    }
    if (takelog == true) {
        return tmp;
    } else {
        return exp(tmp);
    }
}

// [[Rcpp::export]]
double llgev1(NumericMatrix data, NumericMatrix loc, double sc, double sh, bool neg = false){
    double sum = 0.0;
    int data_nrow = data.nrow();
    int data_ncol = data.ncol();
    for (int n = 0; n < data_nrow; n++){
        for (int s = 0; s < data_ncol; s++){
            double tmp = dgevcpp(data(n, s), loc(n, s), sc, sh, true);
            sum += tmp;
        }
    }
    if (neg == true) sum = -sum;
    return sum;
}

// [[Rcpp::export]]
NumericMatrix llgev(NumericMatrix data, NumericMatrix loc, double sc, double sh){
    int data_nrow = data.nrow();
    int data_ncol = data.ncol();
    NumericMatrix mat(data_nrow,data_ncol);
    for (int n = 0; n < data_nrow; n++){
        for (int s = 0; s < data_ncol; s++){
            mat(n,s) = dgevcpp(data(n, s), loc(n, s), sc, sh, true);
            
        }
    }
    return mat;
}


// [[Rcpp::export]]
double get_ed(NumericVector s1, NumericVector s2){
    double t1 = s1(0) - s2(0);
    double t2 = s1(1) - s2(1);
    double tmp = sqrt(t1*t1 + t2*t2);
    return tmp;
}

// [[Rcpp::export]]
NumericVector get_sitecor(NumericVector s1, NumericMatrix vec, double gvar, double gsc){
    int nrow = vec.nrow();
    NumericVector tmp(nrow);
    for (int i = 0; i < nrow; i++){
        NumericVector site_j(2); site_j(0) = vec(i, 0); site_j(1) = vec(i, 1);
        double tmped = get_ed(s1, site_j);
        tmp(i) = exp(-tmped / gsc);
    }
    return tmp;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
NumericMatrix getdistcpp(NumericMatrix site){
    int nrow = site.nrow();
    int ncol = site.nrow();
    NumericVector site_i(nrow), site_j(nrow);
    NumericMatrix mat(nrow,ncol);
    for (int i = 0; i < nrow - 1; i++){
        for (int j = i + 1; j < ncol; j++){
            NumericVector site_i(2); site_i(0) = site(i, 0); site_i(1) = site(i, 1);
            NumericVector site_j(2); site_j(0) = site(j, 0); site_j(1) = site(j, 1);
            // Rcout << i << " " << site_i(0) << " " << site_i(1) << endl;
            // Rcout << j << " " << site_j(0) << " " << site_j(1) << endl;
            mat(i, j) = get_ed(site_i, site_j);
            mat(j, i) = mat(i, j);
        }
    }
    return mat;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
NumericMatrix getexpfromdistcpp(NumericMatrix sitedist, double var, double scale){
    int nrow = sitedist.nrow();
    int ncol = sitedist.ncol();
    NumericMatrix mat(nrow,ncol);
    for (int i = 0; i < nrow ; i++){
        for (int j = i ; j < ncol; j++){
            // Rcout << i << " " << site_i(0) << " " << site_i(1) << endl;
            // Rcout << j << " " << site_j(0) << " " << site_j(1) << endl;
            mat(i, j) = var * exp(-sitedist(i, j) / scale);
            mat(j, i) = mat(i, j);
        }
    }
    return mat;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
NumericMatrix getexpcpp(NumericMatrix site, double var, double scale){
    int nrow = site.nrow();
    int ncol = site.nrow();
    NumericVector site_i(nrow), site_j(nrow);
    NumericMatrix mat(nrow,ncol);
    for (int i = 0; i < nrow - 1; i++){
        for (int j = i + 1; j < ncol; j++){
            NumericVector site_i(2); site_i(0) = site(i, 0); site_i(1) = site(i, 1);
            NumericVector site_j(2); site_j(0) = site(j, 0); site_j(1) = site(j, 1);
            // Rcout << i << " " << site_i(0) << " " << site_i(1) << endl;
            // Rcout << j << " " << site_j(0) << " " << site_j(1) << endl;
            mat(i, j) = var * exp(-get_ed(site_i, site_j) / scale);
            mat(j, i) = mat(i, j);
        }
    }
    for (int i = 0; i < nrow; i++){
        mat(i, i) = var;
    }
    return mat;
}


// ratiodiff <- function(up,down,takeexp=TRUE){
//     tmp <- up - down
//     if (takeexp==TRUE) tmp <- exp(tmp)
//     if (is.nan(tmp)==TRUE) tmp <- 0
//     return(tmp) 
// }

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double rdiff(double up, double down, bool takeexp = true){
    double tmp = up - down;
    if (takeexp == true) tmp = exp(tmp);
    if (isnan(tmp) == true){
        if (takeexp == true) tmp = 0; else tmp = -INFINITY;
    }
    if (takeexp == true) tmp = std::min(tmp, 1.0); else tmp = std::min(tmp, 0.0);
    return tmp;
}

