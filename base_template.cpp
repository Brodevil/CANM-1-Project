#include <bits/stdc++.h>

using namespace std;

vector<vector<double>> matrix_mult(const vector<vector<double>>& a, 
                                         const vector<vector<double>>& b) {
    int p = a.size();       // rows in a
    int q = a[0].size();    // cols in a
    int r = b.size();       // rows in b
    int s = b[0].size();    // cols in b

    if (q != r) return {};

    vector<vector<double>> c(p, vector<double>(s, 0));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < s; j++) {
            for (int k = 0; k < q; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

vector<vector<double>> row_switch(int p, int q, vector<vector<double>> a) {
    swap(a[p-1], a[q-1]);     
    return a;
}

vector<vector<double>> row_subtract(int p, int q, double scalar, vector<vector<double>> a) {
    if (p <= 0 || q <= 0 || p > a.size() || q > a.size()) {
        return a;
    }
    for (int j = 0; j < a[p-1].size(); j++) {
        a[p-1][j] -= (scalar * a[q-1][j]);
    }
    return a;
}

vector<vector<double>> gaussian_elimination(vector<vector<double>> a) {
    int rows = a.size();
    if (rows == 0) return a;
    int cols = a[0].size();

    int pivot_row = 0;
    for (int j = 0; j < cols && pivot_row < rows; j++) {
        int sel = pivot_row;
        for (int i = pivot_row + 1; i < rows; i++) {
            if (abs(a[i][j]) > abs(a[sel][j])) sel = i;
        }

        if (abs(a[sel][j]) < 1e-9) continue;

        a = row_switch(pivot_row + 1, sel + 1, a);

        for (int i = pivot_row + 1; i < rows; i++) {
            double scalar = a[i][j] / a[pivot_row][j];
            a = row_subtract(i + 1, pivot_row + 1, scalar, a);
        }
        pivot_row++;
    }
    return a;
}

vector<vector<double>> make_augmented(vector<vector<double>> A, vector<double> b) {
    for (int i = 0; i < A.size(); i++) {
        A[i].push_back(b[i]); // Add b[i] as the last column
    }
    return A;
}

vector<double> solve_system(vector<vector<double>> augmented) {
    vector<vector<double>> upper = gaussian_elimination(augmented);
    
    int n = upper.size();
    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += upper[i][j] * x[j];
        }
        x[i] = (upper[i][n] - sum) / upper[i][i];
    }
    return x;
}

vector<vector<double>> transpose(const vector<vector<double>>& a) {
    if (a.empty()) return {};

    int rows = a.size();
    int cols = a[0].size();

    vector<vector<double>> result(cols, vector<double>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = a[i][j];
        }
    }
    return result;
}

void qr_graham(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int m = A.size();       
    int n = A[0].size();   

    Q = A;
    R = vector<vector<double>>(n, vector<double>(n, 0.0));

    for (int k = 0; k < n; k++) {
        double norm = 0;
        for (int i = 0; i < m; i++) {
            norm += Q[i][k] * Q[i][k];
        }
        R[k][k] = sqrt(norm);

        for (int i = 0; i < m; i++) {
            if (R[k][k] > 1e-15) { 
                Q[i][k] /= R[k][k];
            }
        }
        for (int j = k + 1; j < n; j++) {
            double dot = 0;
            for (int i = 0; i < m; i++) {
                dot += Q[i][k] * Q[i][j];
            }
            R[k][j] = dot;
            for (int i = 0; i < m; i++) {
                Q[i][j] -= R[k][j] * Q[i][k];
            }
        }
    }
}

double vector_norm(const vector<double>& v, int start_index) {
    double sum = 0;
    for (int i = start_index; i < v.size(); i++) sum += v[i] * v[i];
    return sqrt(sum);
}


double dot(const vector<double>& v1, const vector<double>& v2) {
    double res = 0;
    for(int i = 0; i < v1.size(); i++) res += v1[i] * v2[i];
    return res;
}

vector<double> mat_vec_mult(const vector<vector<double>>& A, const vector<double>& x) {
    int n = A.size();
    vector<double> res(n, 0);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < A[i].size(); j++) {
            res[i] += A[i][j] * x[j];
        }
    }
    return res;
}

// SYMMETRIC A for steepest descent and conjugate gradient

vector<double> stepest_descent(const vector<vector<double>>& A, const vector<double>& b, double tol = 1e-9) {
    int n = A.size();
    vector<double> x(n, 0.0); // Start with guess x = 0
    vector<double> r = b;     // r_0 = b - A(0) = b

    double r_old_dot = dot(r, r);

    while (sqrt(r_old_dot)>=tol) {
        if (sqrt(r_old_dot) < tol) break;

        vector<double> Ar = mat_vec_mult(A, r);

        double alpha = r_old_dot / dot(r, Ar);

        for (int i = 0; i < n; i++) {
            x[i] += alpha * r[i];
        }

        for (int i = 0; i < n; i++) {
            r[i] -= alpha * Ar[i];
        }

        r_old_dot = dot(r, r);
    }
    return x;
}


vector<double> conjugate_gradiant(const vector<vector<double>>& A, const vector<double>& b, double tol = 1e-10) {
    int n = A.size();
    vector<double> x(n, 0.0); 
    vector<double> r = b; // Assuming x0 is 0, so r0 = b - A*0
    vector<double> d = r; // r0 = d0
    
    int iterations = 0;

    while (sqrt(dot(r, r)) >= tol) {
        vector<double> Ad = mat_vec_mult(A, d); 
        double dAd = dot(d, Ad);
        
        if (abs(dAd) < 1e-18) break; 
        
        double alpha = dot(r, d) / dAd;

        for (int i = 0; i < n; i++) {
            x[i] = x[i] + alpha * d[i];
        }

        for (int i = 0; i < n; i++) {
            r[i] = r[i] - alpha * Ad[i];
        }

        if (sqrt(dot(r, r)) < tol) break;


        double beta = -dot(r, Ad) / dAd;

        for (int i = 0; i < n; i++) {
            d[i] = r[i] + beta * d[i];
        }

        iterations++;
    }

    // cout << "Converged in " << iterations << " iterations." << endl;
    return x;
}

int main(){

    /*
      Uncomment to test
    */
    // vector<vector<double>> a = {{1, 2, 3}, {4, 5, 6}};
    // vector<vector<double>> b = {{7, 8}, {9, 10}, {11, 12}};

    // vector<vector<double>> c = matrix_mult(a,b);
    // // c = row_switch(1,2,c);
    // // c = row_subtract(1,2, (double)139/58, c);
    // c = gaussian_elimination(c);
    // for(int i = 0; i<c.size(); i++){
    //     for(int j = 0; j<c[0].size(); j++){
    //         cout << c[i][j] << ", ";
    //     }
    //     cout << endl;
    // }

    // vector<vector<double>> a = {{1, 1, 1}, {0, 2, 5}, {2, 5, -1}};
    // vector<double> b = {6, -4, 27};
    // vector<vector<double>> aug = make_augmented(a,b);
    // vector<double> soln = solve_system(aug);
    // for(int i= 0; i<soln.size(); i++){
    //     cout << soln[i] << " ";
    // }
  
  return 0;

}
