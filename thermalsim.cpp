#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <fstream>
#include <string>

using namespace std;

#define PI 3.14159265359

double get_dt(double alpha, double h) {
    return (h * h) / (4.0 * alpha);
}

int main() {
    int row = 100, col = 100;
    double L = 0.1; // 10 cm plate
    double h = L / (col - 1);
    double tp, btm, lf, rg, T0;

    vector<vector<double>> T(row, vector<double>(col, 0.0));
    vector<vector<double>> S(row, vector<double>(col, 0.0));

    cout << "Boundary conditions (TOP BOTTOM LEFT RIGHT): ";
    cin >> tp >> btm >> lf >> rg;
    cout << "Initial plate temperature: ";
    cin >> T0;

    cout << "\nPlate material:\n"
         << "1. Aluminium\n"
         << "2. Copper\n"
         << "3. Iron\n"
         << "4. Stainless Steel\n"
         << "5. Brass\n"
         << "6. Custom\n"
         << "Choice: ";
    int choice;
    cin >> choice;

    double kc, rho_m, cp;
    string mat_name;

    switch (choice) {
        case 1:
            mat_name = "Aluminium";
            kc = 237.0; rho_m = 2700.0; cp = 897.0;
            break;
        case 2:
            mat_name = "Copper";
            kc = 401.0; rho_m = 8960.0; cp = 385.0;
            break;
        case 3:
            mat_name = "Iron";
            kc = 80.0; rho_m = 7870.0; cp = 447.0;
            break;
        case 4:
            mat_name = "Stainless Steel";
            kc = 16.0; rho_m = 8000.0; cp = 500.0;
            break;
        case 5:
            mat_name = "Brass";
            kc = 109.0; rho_m = 8500.0; cp = 380.0;
            break;
        default:
            mat_name = "Custom";
            cout << "k (W/m.K): "; cin >> kc;
            cout << "Density (kg/m^3): "; cin >> rho_m;
            cout << "Specific heat (J/kg.K): "; cin >> cp;
            break;
    }

    double alpha = kc / (rho_m * cp);
    cout << "\nMaterial: " << mat_name << " (alpha = " << alpha << " m^2/s)\n";

    // Set initial temp
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            T[i][j] = T0;
            S[i][j] = T0;
        }
    }

    // Outer boundaries
    for (int j = 0; j < col; j++) {
        T[0][j] = btm;
        S[0][j] = btm;
        T[row - 1][j] = tp;
        S[row - 1][j] = tp;
    }
    for (int i = 1; i < row - 1; i++) {
        T[i][0] = lf;
        S[i][0] = lf;
        T[i][col - 1] = rg;
        S[i][col - 1] = rg;
    }

    double tol = 1e-5;
    double dt = 0.9 * get_dt(alpha, h);
    double r = (alpha * dt) / (h * h);
    long long steps = 0;
    double max_rate;
    vector<vector<double>> T_next = T;

    // Snapshot setup
    double tau = (L * L) / (2.0 * PI * PI * alpha);
    double mult[6] = {0.0, 0.1, 0.5, 1.0, 2.0, 4.0};

    vector<double> snap_t;
    vector<vector<vector<double>>> snaps;
    snap_t.push_back(0.0);
    snaps.push_back(T);
    int idx = 1;

    // Transient FTCS loop
    do {
        steps++;
        max_rate = 0.0;

        for (int i = 1; i < row - 1; i++) {
            for (int j = 1; j < col - 1; j++) {
                T_next[i][j] = T[i][j] + r * (T[i - 1][j] + T[i + 1][j] + T[i][j - 1] + T[i][j + 1] - 4.0 * T[i][j]);

                double rate = fabs(T_next[i][j] - T[i][j]) / dt;
                if (rate > max_rate) max_rate = rate;
            }
        }
        T = T_next;

        double curr_t = steps * dt;
        if (idx < 6 && curr_t >= mult[idx] * tau) {
            snap_t.push_back(curr_t);
            snaps.push_back(T);
            idx++;
        }
    } while (max_rate > tol && steps < 20000000);

    double t_final = steps * dt;
    snap_t.push_back(t_final);
    snaps.push_back(T);

    // Steady-state Gauss-Seidel + SOR
    double omega = 1.85;
    int n = 0;
    double max_change;

    do {
        n++;
        max_change = 0.0;

        for (int i = 1; i < row - 1; i++) {
            for (int j = 1; j < col - 1; j++) {
                double S_old = S[i][j];
                double S_gs = (S[i - 1][j] + S[i + 1][j] + S[i][j - 1] + S[i][j + 1]) / 4.0;
                S[i][j] = (1.0 - omega) * S_old + omega * S_gs;

                double diff = fabs(S[i][j] - S_old);
                if (diff > max_change) max_change = diff;
            }
        }
    } while (max_change > tol && n < 1000000);

    double max_diff = 0.0;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            double d = fabs(T[i][j] - S[i][j]);
            if (d > max_diff) max_diff = d;
        }
    }

    cout << "\nTransient converged in " << steps << " steps (dt = " << dt << " s, t = " << t_final << " s)\n";
    for (int i = row - 1; i >= 0; i--) {
        for (int j = 0; j < col; j++) {
            printf("%.2f\t", T[i][j]);
        }
        cout << "\n";
    }

    cout << "\nSteady-state converged in " << n << " iterations\n";
    for (int i = row - 1; i >= 0; i--) {
        for (int j = 0; j < col; j++) {
            printf("%.2f\t", S[i][j]);
        }
        cout << "\n";
    }

    cout << "\nCenter T -> transient: " << T[row / 2][col / 2] << ", steady: " << S[row / 2][col / 2] << "\n";
    cout << "Max diff: " << max_diff << "\n";
    cout << "Tau: " << tau << " s | Snapshots at t = ";
    for (int k = 0; k < (int)snap_t.size(); k++) {
        cout << snap_t[k];
        if (k + 1 < (int)snap_t.size()) cout << ", ";
        else cout << " s\n";
    }

    cout << "\nOutput filename (without .csv): ";
    string fname;
    cin >> fname;

    string csv_file = fname + ".csv";
    ofstream fptr(csv_file);
    fptr << "x,y,T_transient,T_steady";
    for (int k = 0; k < (int)snap_t.size(); k++) {
        fptr << ",T@" << snap_t[k];
    }
    fptr << "\n";

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            fptr << j << "," << i << "," << T[i][j] << "," << S[i][j];
            for (int k = 0; k < (int)snap_t.size(); k++) {
                fptr << "," << snaps[k][i][j];
            }
            fptr << "\n";
        }
    }
    fptr.close();
    cout << "Saved " << csv_file << "\n";

    string cmd = "python visualizer.py " + csv_file;
    if (system(cmd.c_str()) == 0) {
        string open_cmd = "start " + fname + ".png";
        system(open_cmd.c_str());
    }

    return 0;
}