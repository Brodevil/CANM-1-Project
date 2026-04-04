#include <bits/stdc++.h>

using namespace std;

constexpr double PRINT_TOLERANCE = 1e-9;

struct NewtonInterpolation {
    vector<double> x_values;
    vector<vector<double>> table;
};

double test_function(double x) {
    return (x * x - 1.0) * (x - 1.0);
}

bool is_near_zero(double value, double tolerance = PRINT_TOLERANCE) {
    return abs(value) < tolerance;
}

string format_number(double value) {
    if (is_near_zero(value)) {
        value = 0.0;
    }

    ostringstream stream;
    stream << fixed << setprecision(10) << value;

    string formatted = stream.str();
    while (!formatted.empty() && formatted.back() == '0') {
        formatted.pop_back();
    }
    if (!formatted.empty() && formatted.back() == '.') {
        formatted.pop_back();
    }
    if (formatted == "-0") {
        formatted = "0";
    }

    return formatted;
}

string format_factor(double node) {
    if (is_near_zero(node)) {
        return "x";
    }
    if (node < 0.0) {
        return "(x + " + format_number(abs(node)) + ")";
    }
    return "(x - " + format_number(node) + ")";
}

NewtonInterpolation build_divided_difference_table(
    const vector<double>& x_values,
    const vector<double>& y_values
) {
    int n = x_values.size();

    if (n == 0 || y_values.size() != x_values.size()) {
        throw invalid_argument("x and y must be non-empty and of the same size.");
    }

    vector<vector<double>> table(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        table[i][0] = y_values[i];
    }

    for (int order = 1; order < n; order++) {
        for (int row = 0; row < n - order; row++) {
            double denominator = x_values[row + order] - x_values[row];
            if (abs(denominator) < 1e-12) {
                throw invalid_argument("x values must be distinct.");
            }

            table[row][order] =
                (table[row + 1][order - 1] - table[row][order - 1]) / denominator;
        }
    }

    return {x_values, table};
}

vector<double> extract_newton_coefficients(const NewtonInterpolation& interpolation) {
    int n = interpolation.x_values.size();
    vector<double> coefficients(n);

    for (int i = 0; i < n; i++) {
        coefficients[i] = interpolation.table[0][i];
    }

    return coefficients;
}

double evaluate_newton_polynomial(
    const NewtonInterpolation& interpolation,
    double value
) {
    vector<double> coefficients = extract_newton_coefficients(interpolation);
    int n = coefficients.size();

    double result = coefficients[0];
    double product_term = 1.0;

    for (int i = 1; i < n; i++) {
        product_term *= (value - interpolation.x_values[i - 1]);
        result += coefficients[i] * product_term;
    }

    return result;
}

void print_divided_difference_table(const NewtonInterpolation& interpolation) {
    int n = interpolation.x_values.size();

    cout << "\nDivided Difference Table\n";
    for (int i = 0; i < n; i++) {
        cout << "x" << i << ":";
        if(i<10) cout << " ";
        for (int j = 0; j < n - i; j++) {
            
            cout << setw(14) << interpolation.table[i][j] << ' ';
        }
        cout << '\n';
    }
}

void print_newton_polynomial(const NewtonInterpolation& interpolation) {
    vector<double> coefficients = extract_newton_coefficients(interpolation);

    cout << "\nNewton Interpolating Polynomial:\nP(x) = ";

    bool first_term = true;
    for (int i = 0; i < coefficients.size(); i++) {
        double coefficient = coefficients[i];
        if (is_near_zero(coefficient)) {
            continue;
        }

        double magnitude = abs(coefficient);
        if (!first_term) {
            cout << (coefficient < 0.0 ? " - " : " + ");
        } else if (coefficient < 0.0) {
            cout << '-';
        }

        bool printed_coefficient = false;
        if (i == 0 || !is_near_zero(magnitude - 1.0)) {
            cout << format_number(magnitude);
            printed_coefficient = true;
        }

        for (int j = 0; j < i; j++) {
            if (printed_coefficient || j > 0) {
                cout << ' ';
            }
            cout << format_factor(interpolation.x_values[j]);
            printed_coefficient = true;
        }

        first_term = false;
    }

    if (first_term) {
        cout << '0';
    }
    cout << '\n';
}

vector<double> expand_newton_polynomial(const NewtonInterpolation& interpolation) {
    vector<double> coefficients = extract_newton_coefficients(interpolation);
    vector<double> expanded(1, 0.0);
    vector<double> basis(1, 1.0);

    for (int i = 0; i < coefficients.size(); i++) {
        double coefficient = is_near_zero(coefficients[i]) ? 0.0 : coefficients[i];

        if (!is_near_zero(coefficient)) {
            if (expanded.size() < basis.size()) {
                expanded.resize(basis.size(), 0.0);
            }

            for (int j = 0; j < basis.size(); j++) {
                expanded[j] += coefficient * basis[j];
            }
        }

        if (i + 1 < coefficients.size()) {
            vector<double> next_basis(basis.size() + 1, 0.0);
            for (int j = 0; j < basis.size(); j++) {
                next_basis[j] -= interpolation.x_values[i] * basis[j];
                next_basis[j + 1] += basis[j];
            }
            basis = move(next_basis);
        }
    }

    for (double& coefficient : expanded) {
        if (is_near_zero(coefficient)) {
            coefficient = 0.0;
        }
    }

    while (expanded.size() > 1 && is_near_zero(expanded.back())) {
        expanded.pop_back();
    }

    return expanded;
}

void print_expanded_polynomial(const vector<double>& coefficients) {
    

    bool first_term = true;
    for (int degree = coefficients.size() - 1; degree >= 0; degree--) {
        double coefficient = coefficients[degree];
        if (is_near_zero(coefficient)) {
            continue;
        }

        double magnitude = abs(coefficient);
        if (!first_term) {
            cout << (coefficient < 0.0 ? " - " : " + ");
        } else if (coefficient < 0.0) {
            cout << '-';
        }

        bool show_coefficient = (degree == 0) || !is_near_zero(magnitude - 1.0);
        if (show_coefficient) {
            cout << format_number(magnitude);
        }

        if (degree >= 1) {
            if (show_coefficient) {
                cout << ' ';
            }
            cout << 'x';
            if (degree >= 2) {
                cout << '^' << degree;
            }
        }

        first_term = false;
    }

    if (first_term) {
        cout << '0';
    }
    cout << '\n';
}

vector<double> differentiate_polynomial(const vector<double>& poly) {
    int n = poly.size();

    // If constant or empty → derivative is 0
    if (n <= 1) {
        return {0.0};
    }

    vector<double> derivative(n - 1);

    for (int i = 1; i < n; i++) {
        derivative[i - 1] = i * poly[i];
    }

    return derivative;
}

int main() {
    cout << fixed << setprecision(10);

    vector<double> x_values = {
        -2.0, -1.8333333333333333, -1.6666666666666667, -1.5,
        -1.3333333333333335, -1.1666666666666667, -1.0, -0.8333333333333335,
        -0.6666666666666667, -0.5, -0.3333333333333335, -0.16666666666666674,
        0.0, 0.16666666666666652, 0.33333333333333304, 0.5,
        0.6666666666666665, 0.833333333333333, 1.0, 1.1666666666666665,
        1.333333333333333, 1.5, 1.6666666666666665, 1.833333333333333, 2.0
    };

    vector<double> y_values = {
        -9.0, -6.689814814814813, -4.740740740740742, -3.125,
        -1.814814814814816, -0.7824074074074079, 0.0, 0.5601851851851848,
        0.9259259259259258, 1.125, 1.1851851851851851, 1.1342592592592593,
        1.0, 0.8101851851851853, 0.592592592592593, 0.375,
        0.18518518518518534, 0.0509259259259261, 0.0, 0.06018518518518507,
        0.25925925925925875, 0.625, 1.1851851851851847, 1.967592592592591, 3.0
    };

    double query_point = 0.3;

    try {
        NewtonInterpolation interpolation =
            build_divided_difference_table(x_values, y_values);

        cout << "Testing Newton Divided Difference Interpolation\n";
        cout << "Using hardcoded data points from group6.csv\n";
        cout << "Reference function: f(x) = (x^2 - 1)(x - 1)\n";
        cout << "Evaluation point: x = " << query_point << '\n';
        cout << "Data points used:\n";
        for (int i = 0; i < x_values.size(); i++) {
            cout << "x" << i << " = " << x_values[i] << ", y"<< i << " = " << y_values[i] << '\n';
        }

        print_divided_difference_table(interpolation);
        print_newton_polynomial(interpolation);
        cout << "\nStandard Polynomial Form:\nP(x) = ";
        print_expanded_polynomial(expand_newton_polynomial(interpolation));

        cout << "\nDifferentiation of P(x):\nP'(x) = ";
        print_expanded_polynomial(differentiate_polynomial(expand_newton_polynomial(interpolation)));

        double interpolated_value =
            evaluate_newton_polynomial(interpolation, query_point);
        double exact_value = test_function(query_point);

        cout << "\nInterpolated value at x = " << query_point
             << " is " << interpolated_value << '\n';
        cout << "Exact value at x = " << query_point
             << " is " << exact_value << '\n';
        cout << "Absolute error = " << abs(interpolated_value - exact_value) << '\n';
    } catch (const exception& error) {
        cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}

