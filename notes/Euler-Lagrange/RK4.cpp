#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// Parámetros físicos consistentes
const double m = 0.1;
const double k = 1.0;
const double omega = sqrt(k / m); // omega = sqrt(10) ~ 3.1622
const double omega2 = k / m;

// Parámetros de simulación
const double h = 0.01;
const double T = 50.0;

struct State {
    double x;
    double v;
};

int menu() {
    int opc;
    cout << "=== EDO Oscilador Armónico Simple ===" << endl;
    cout << "1. Solución Analítica Pura" << endl;
    cout << "2. Comparar Euler Explícito vs Analítica" << endl;
    cout << "3. Comparar Euler Mejorado (Heun) vs Analítica" << endl;
    cout << "4. Comparar Runge-Kutta 4 (RK4) vs Analítica" << endl;
    cout << "Seleccione una opción: ";
    cin >> opc;
    return opc;
}

// Derivadas del sistema hamiltoniano
State deriv(const State& s) {
    return {s.v, -omega2 * s.x};
}

// Solución exacta bajo condiciones iniciales x(0)=1, v(0)=0
State analytical(double t) {
    return {cos(omega * t), -omega * sin(omega * t)};
}

// Integradores Numéricos
State euler_step(State s) {
    State ds = deriv(s);
    return {s.x + ds.x * h, s.v + ds.v * h};
}

State euler_mejorado_step(State s) {
    State ds1 = deriv(s);
    State s_pred = {s.x + ds1.x * h, s.v + ds1.v * h};
    State ds2 = deriv(s_pred);
    return {
        s.x + 0.5 * h * (ds1.x + ds2.x),
        s.v + 0.5 * h * (ds1.v + ds2.v)
    };
}

State rk4_step(State s) {
    State k1 = deriv(s);
    
    State s2 = {s.x + 0.5 * h * k1.x, s.v + 0.5 * h * k1.v};
    State k2 = deriv(s2);
    
    State s3 = {s.x + 0.5 * h * k2.x, s.v + 0.5 * h * k2.v};
    State k3 = deriv(s3);
    
    State s4 = {s.x + h * k3.x, s.v + h * k3.v};
    State k4 = deriv(s4);
    
    return {
        s.x + (h / 6.0) * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x),
        s.v + (h / 6.0) * (k1.v + 2.0 * k2.v + 2.0 * k3.v + k4.v)
    };
}

double compute_energy(const State& s) {
    return 0.5 * m * s.v * s.v + 0.5 * k * s.x * s.x;
}

// Generación del archivo analítico de referencia basal
void generate_analytical_reference() {
    ofstream file("datos_analiticos.dat");
    for (double t = 0.0; t <= T; t += h) {
        State exact = analytical(t);
        file << t << " " << exact.x << " " << exact.v << " " << compute_energy(exact) << "\n";
    }
    file.close();
}

// Script de Gnuplot multi-curva (Muestra x y v tanto numéricas como analíticas)
void generate_gnuplot_script(const string& numeric_file, const string& method_name) {
    ofstream gp("osciladorarmonico.gp");
    gp << "set term pngcairo size 1000,700 font 'Arial,10'\n";
    gp << "set output 'osciladorarmonico.png'\n";
    gp << "set title 'Harmonic oscilator' font 'Arial,12'\n";
    gp << "set xlabel 'Time (t)'\n";
    gp << "set ylabel 'X(t)'\n";
    gp << "set grid lc rgb '#d3d3d3' lt 1 lw 0.5\n";
    
    // Definición de estilos de línea (Líneas continuas para numérico, discontinuas para analítico)
    gp << "set style line 1 lc rgb '#0066cc' lt 1 lw 2\n";  // x Numérico
    gp << "set style line 2 lc rgb '#cc0000' lt 1 lw 2\n";  // v Numérico
    gp << "set style line 3 lc rgb '#404040' lt 0 lw 1.5\n"; // x Analítico (Punteado)
    gp << "set style line 4 lc rgb '#808080' lt 0 lw 1.5\n"; // v Analítico (Punteado)
    
    if (numeric_file.empty()) {
        // Caso analítico puro
        gp << "plot 'datos_analiticos.dat' using 1:2 with lines title 'x(t) Analitycal' lc rgb '#0066cc' lw 2, \\\n";
        gp << "     'datos_analiticos.dat' using 1:3 with lines title 'v(t) Analitycal' lc rgb '#cc0000' lw 2\n";
    } else {
        // Superposición de curvas
        gp << "plot 'datos_analiticos.dat' using 1:2 with lines linestyle 3 title 'x(t) Analytical', \\\n";
        gp << "     'datos_analiticos.dat' using 1:3 with lines linestyle 4 title 'v(t) Analytical', \\\n";
        gp << "     '" << numeric_file << "' using 1:2 with lines linestyle 1 title 'x(t) " << method_name << "', \\\n";
        gp << "     '" << numeric_file << "' using 1:3 with lines linestyle 2 title 'v(t) " << method_name << "'\n";
    }
    gp.close();
}

int main() {
    int choice = menu();
    State current_state = {1.0, 0.0}; 
    string filename = "";
    string method_name = "";

    // Generar siempre la referencia exacta independientemente de la opción elegida
    generate_analytical_reference();

    switch (choice) {
        case 1:
            method_name = "Analítica";
            break;

        case 2:
            filename = "datos_euler.dat";
            method_name = "Euler Explícito";
            {
                ofstream file(filename);
                for (double t = 0.0; t <= T; t += h) {
                    file << t << " " << current_state.x << " " << current_state.v << " " << compute_energy(current_state) << "\n";
                    current_state = euler_step(current_state);
                }
                file.close();
            }
            break;

        case 3:
            filename = "datos_euler_mejorado.dat";
            method_name = "Euler Mejorado (Heun)";
            {
                ofstream file(filename);
                for (double t = 0.0; t <= T; t += h) {
                    file << t << " " << current_state.x << " " << current_state.v << " " << compute_energy(current_state) << "\n";
                    current_state = euler_mejorado_step(current_state);
                }
                file.close();
            }
            break;

        case 4:
            filename = "datos_rk4.dat";
            method_name = "Runge-Kutta 4";
            {
                ofstream file(filename);
                for (double t = 0.0; t <= T; t += h) {
                    file << t << " " << current_state.x << " " << current_state.v << " " << compute_energy(current_state) << "\n";
                    current_state = rk4_step(current_state);
                }
                file.close();
            }
            break;

        default:
            cout << "Opción inválida." << endl;
            return 1;
    }

    cout << "Cálculo completado." << endl;
    generate_gnuplot_script(filename, method_name);
    
    cout << "Invocando Gnuplot para renderizar 'osciladorarmonico.png'..." << endl;
    system("gnuplot osciladorarmonico.gp");

    return 0;
}
