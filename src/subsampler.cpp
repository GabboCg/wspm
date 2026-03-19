// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Tent kernel: min(x, 1-x)
inline double tent(double x) {
    
    return std::min(x, 1.0 - x);
    
}

// Pre-averaged return series (tent kernel hard-coded)
// logp: log-price vector, K: bandwidth
// [[Rcpp::export]]
arma::vec aux_preavgk(const arma::vec& logp, int K) {
    
    arma::vec r = arma::diff(logp);
    
    int n = r.n_elem;
    int out_n = n - K + 2;
    arma::vec r_pa(out_n, arma::fill::zeros);
    
    for (int i = 1; i <= K - 1; i++) {
        
        double w = tent((double)i / K);
        r_pa += w * r.subvec(i - 1, i - 1 + out_n - 1);
        
    }
    
    return r_pa;
    
}

// E[|Z|^r], Z ~ N(0,1), product over elements of vpow
// [[Rcpp::export]]
double f_mu(const arma::vec& vpow) {
    
    double out = 1.0;
    
    for (int i = 0; i < (int)vpow.n_elem; i++) {
        
        double ri = vpow[i];
        out *= std::pow(2.0, 0.5 * ri) * R::gammafn(0.5 * (ri + 1.0)) / R::gammafn(0.5);
        
    }
    
    return out;
    
}

// Kernel constants psi1, psi2 for bias correction (tent kernel)
// [[Rcpp::export]]
Rcpp::List f_psi(int K) {
    
    arma::vec g(K + 1);
    
    for (int i = 0; i <= K; i++) {
        
        g[i] = tent((double)i / K);
        
    }
    
    arma::vec dg = arma::diff(g);
    
    double psi1 = arma::sum(dg % dg) * K;
    double psi2 = arma::sum(g % g) / K;
    
    return Rcpp::List::create(Rcpp::Named("psi1") = psi1, Rcpp::Named("psi2") = psi2);
    
}

// Long-run noise variance estimator (Jacod, Li, Zheng 2019)
// [[Rcpp::export]]
double omega2(const arma::vec& logp) {
    
    int n = logp.n_elem;      // total number of prices
    int nd = n - 1;           // number of returns (n_ in R)
    int ell1 = (int)std::round(std::pow((double)nd, 1.0 / 5.0));
    int ell2 = (int)std::round(std::pow((double)nd, 1.0 / 8.0));

    // Forward moving average: logp_avg[i] = mean(logp[i..i+ell1-1])
    arma::vec logp_avg(n, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
        
        int end = std::min(i + ell1 - 1, n - 1);
        logp_avg[i] = arma::mean(logp.subvec(i, end));
        
    }

    // Um[m] for m = 0..ell1
    arma::vec Um(ell1 + 1, arma::fill::zeros);
    
    // Indices (0-based):
    //   a = logp[0..n-5*ell1-1] - logp_avg[2*ell1..n-3*ell1-1]
    //   b = logp[m..n-5*ell1-1+m] - logp_avg[m+4*ell1..n-ell1-1+m]
    int len_a = n - 5 * ell1; // number of elements in a/b
    if (len_a <= 0) {
        
        // Fallback: return a small positive value
        return 1e-10;
        
    }
    
    arma::vec base_a = logp.subvec(0, len_a - 1) - logp_avg.subvec(2 * ell1, n - 3 * ell1 - 1);

    for (int m = 0; m <= ell1; m++) {
        
        arma::vec b = logp.subvec(m, len_a - 1 + m) - logp_avg.subvec(m + 4 * ell1, n - ell1 - 1 + m);
        Um[m] = arma::mean(base_a % b);
        
    }

    // Bartlett-weighted sum: omega2 = Um[0] + 2 * sum_{j=1}^{ell2-1} w[j]*Um[j]
    double omega2 = Um[0];
    for (int j = 1; j <= ell2 - 1; j++) {
        
        double w = 1.0 - (double)j / ell2;
        omega2 += 2.0 * w * Um[j];
        
    }
    
    return omega2;
    
}

// Subsampled power-variation + covariance Sigma*
// [[Rcpp::export]]
Rcpp::List f_subsampler(const arma::vec& logp, int K, int p, int L, const arma::mat& mQ, double q, double varpi) {
    
    arma::vec r = arma::diff(logp);
    
    int nObs = (int)r.n_elem;
    double dt = 1.0 / nObs;

    // Adaptive p/L reduction
    int c = (int)std::floor((double)std::floor((double)nObs / (p * K)) / L);
    
    while (c == 0) {
        
        if (p > L) {
            
            p--;
            
        } else {
            
            L--;
            
        }
        
        c = (int)std::floor((double)std::floor((double)nObs / (p * K)) / L);
        
    }

    // Full-sample pre-averaged returns (scaled)
    arma::vec r_pa = std::pow((double)nObs, 0.25) * aux_preavgk(logp, K);
    r_pa = arma::abs(r_pa);

    // Truncation threshold
    int out_n_full = nObs - 2 * K + 2;  // length of r_pa
    arma::uvec idx_pa = arma::regspace<arma::uvec>(0, out_n_full - 2 * K + 1 - 1);
    
    // idx_pa corresponds to R's idx_pa = 1:(nObs - 2*K + 2) -> 0-based: 0..(out_n_full - 2*K + 1)
    // Actually: out_n_full = nObs - K + 2 - 1 + 1... 
    // r_pa has length = (nObs - K + 2) elements: diff(logp) has nObs elements,
    // aux_preavgk returns out_n = nObs - K + 2 elements.
    // idx_pa in R: 1:(nObs - 2*K + 2)  -> 0-based: 0..(nObs - 2*K + 1)
    int idx_pa_len = nObs - 2 * K + 2;
    if (idx_pa_len < 1) idx_pa_len = 1;
    arma::uvec idx_pa_vec = arma::regspace<arma::uvec>(0, idx_pa_len - 1);
    arma::uvec idx_pa_K_vec = idx_pa_vec + K;  // idx_pa + K (0-based)

    arma::vec mu1_vec = {1.0, 1.0};
    double mu1 = f_mu(mu1_vec);
    arma::vec y_bv = r_pa.elem(idx_pa_vec) % r_pa.elem(idx_pa_K_vec);
    double pbv = arma::mean(y_bv) / mu1;

    double cutoff = q * std::sqrt(pbv) * std::pow((double)nObs, 0.25 - varpi);
    arma::vec r_pa_t = r_pa;
    arma::uvec above = arma::find(r_pa_t > cutoff);
    r_pa_t.elem(above).zeros();

    int m = (int)mQ.n_rows;
    arma::rowvec vn(m, arma::fill::zeros);
    arma::mat vln(L, m, arma::fill::zeros);

    // Block matrix rb: (K*p) x (c*L) columns
    int maxr = K * p * c * L;
    arma::vec r_block = r.subvec(0, maxr - 1);
    arma::mat rb(K * p, c * L);
    
    for (int j = 0; j < c * L; j++) {
        
        rb.col(j) = r_block.subvec(j * K * p, j * K * p + K * p - 1);
        
    }

    // idx_pa_b for subsampled: 1:(K*p - 2*K + 2) -> 0-based: 0..(K*p - 2*K + 1)
    int idx_pa_b_len = K * p - 2 * K + 2;
    if (idx_pa_b_len < 1) idx_pa_b_len = 1;
    arma::uvec idx_pa_b_vec = arma::regspace<arma::uvec>(0, idx_pa_b_len - 1);
    arma::uvec idx_pa_b_K_vec = idx_pa_b_vec + K;

    for (int k = 0; k < m; k++) {
        
        arma::vec mQ_row = mQ.row(k).t();
        double mu_k = f_mu(mQ_row);

        arma::vec y_full = arma::pow(r_pa_t.elem(idx_pa_vec), mQ(k, 0)) % arma::pow(r_pa_t.elem(idx_pa_K_vec), mQ(k, 1));
        vn[k] = arma::mean(y_full) / mu_k;

        for (int l = 0; l < L; l++) {
            
            // Gather columns l, l+L, l+2L, ..., l+(c-1)*L from rb
            arma::mat rb_sub(K * p, c);
            
            for (int ci = 0; ci < c; ci++) {
                
                rb_sub.col(ci) = rb.col(l + ci * L);
                
            }

            // logp_ss: prepend 0, cumsum per column -> (K*p+1) x c
            arma::mat logp_ss(K * p + 1, c, arma::fill::zeros);
            
            for (int ci = 0; ci < c; ci++) {
                
                logp_ss.col(ci).subvec(1, K * p) = arma::cumsum(rb_sub.col(ci));
                
            }

            // r_pa_ss: apply aux_preavgk_rcpp to each column
            int col_out_n = K * p - K + 2;  // output length per column
            arma::mat r_pa_ss_mat(col_out_n, c);
            
            for (int ci = 0; ci < c; ci++) {
                
                r_pa_ss_mat.col(ci) = std::pow((double)nObs, 0.25) * aux_preavgk(logp_ss.col(ci), K);
                
            }
            
            r_pa_ss_mat = arma::abs(r_pa_ss_mat);

            // Truncate
            arma::mat r_pa_t_ss = r_pa_ss_mat;
            arma::uvec above_ss = arma::find(r_pa_t_ss > cutoff);
            r_pa_t_ss.elem(above_ss).zeros();

            // Power variation for this subsample
            arma::mat y_ss = arma::pow(r_pa_t_ss.rows(idx_pa_b_vec), mQ(k, 0)) % arma::pow(r_pa_t_ss.rows(idx_pa_b_K_vec), mQ(k, 1));
           
            // mean over columns, then mean of those means
            double col_mean_sum = 0.0;
            
            for (int ci = 0; ci < c; ci++) {
                
                col_mean_sum += arma::mean(y_ss.col(ci));
                
            }
            
            vln(l, k) = (col_mean_sum / c) / mu_k;
            
        }
        
    }

    // d_vln = (vln - vn) * sqrt(c * (K*p - 2*K + 2) * sqrt(dt))
    arma::mat d_vln = vln;
    
    for (int k = 0; k < m; k++) {
        
        d_vln.col(k) -= vn[k];
        
    }
    
    double scale = std::sqrt((double)c * (double)(K * p - 2 * K + 2) * std::sqrt(dt));
    d_vln *= scale;

    arma::mat sigma = (d_vln.t() * d_vln) / L;

    // Recompute RV (mQ row with c(2,0)) without truncation
    for (int k = 0; k < m; k++) {
        
        if (mQ(k, 0) == 2.0 && mQ(k, 1) == 0.0) {
            
            arma::vec y_rv = arma::pow(r_pa.elem(idx_pa_vec), 2.0);
            vn[k] = arma::mean(y_rv);
            
        }
        
    }

    return Rcpp::List::create(Rcpp::Named("vn") = vn, Rcpp::Named("sigma") = sigma);
    
}

// Main intraday jump test
// [[Rcpp::export]]
Rcpp::List intradayJumpTest_cpp(const arma::vec& price, double theta, int p, int L, const arma::mat& mQ, double q, double varpi) {
    
    int nObs = (int)price.n_elem;
    arma::vec logp = arma::log(price);

    // kn: bandwidth, forced even
    int kn = (int)std::round(theta * std::sqrt((double)nObs));
    kn = kn + kn % 2;

    // Long-run noise variance
    double omega2 = omega2(logp);

    // Subsampled covariance
    Rcpp::List sub = f_subsampler(logp, kn, p, L, mQ, q, varpi);
    arma::rowvec vn = sub["vn"];
    arma::mat Sigma = sub["sigma"];

    // Locate BV (first non-(2,0)) and RV (2,0) columns
    int bv_col = -1, rv_col = -1;
    int m = (int)mQ.n_rows;
    
    for (int k = 0; k < m; k++) {
        
        if (mQ(k, 0) == 2.0 && mQ(k, 1) == 0.0) {
            
            rv_col = k;
            
        } else if (bv_col == -1) {
            
            bv_col = k;
            
        }
        
    }
    
    // Default fallback: BV = col 0, RV = col 1
    if (bv_col == -1) bv_col = 0;
    if (rv_col == -1) rv_col = 1;

    double BV_raw = vn[bv_col];
    double RV_raw = vn[rv_col];

    double JV = 100.0 * (1.0 - BV_raw / RV_raw);

    double avar = Sigma(rv_col, rv_col) + Sigma(bv_col, bv_col) - 2.0 * Sigma(rv_col, bv_col);

    double theta_o = (double)kn / std::sqrt((double)nObs);
    double q_critical = R::qnorm(1.0 - 0.01 / nObs, 0.0, 1.0, 1, 0);
    double t_stat = std::pow((double)(nObs - 2 * kn + 2), 0.25) * (RV_raw - BV_raw) / std::sqrt(avar);

    double JF = (t_stat > q_critical) ? 1.0 : 0.0;

    // Bias correction
    Rcpp::List psi_list = f_psi(kn);
    
    double psi1 = psi_list["psi1"];
    double psi2 = psi_list["psi2"];
    double bias = (psi1 / (psi2 * theta_o * theta_o)) * omega2;

    double RV = std::max(((RV_raw / psi2) / theta_o) - bias, 0.0);
    RV = 100.0 * std::sqrt(250.0 * RV);

    double BV = std::max(((BV_raw / psi2) / theta_o) - bias, 0.0);
    BV = 100.0 * std::sqrt(250.0 * BV);

    return Rcpp::List::create(Rcpp::Named("RV") = RV,
                              Rcpp::Named("BV") = BV,
                              Rcpp::Named("JV") = JV,
                              Rcpp::Named("JF") = JF);
    
}
