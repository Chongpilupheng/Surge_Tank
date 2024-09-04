#include <stdio.h>
#include <math.h>

#define G 9.81 //gravitational acceleration m2/s
#define At 23.25 //tunnel cross-sectional area m2
#define L 1964 //pipe length m
#define QTi 56 //intiial turbine flow m3/s
#define hf 1.22 //tunnel head loss m
#define QTf 112 //final turbine flow m3/s
#define TG 5 //time for linearly changing flow from QTi to QTf s
#define AS 148.8 //horizontal tank area at elevation m2
#define Z1 1.22

double derivative(double Qt) {
    double c = hf / pow(QTi, 2);
    return (G * At * (Z1 - c * Qt * fabs(Qt))) / L;
}

double runge_kutta(double Qt, double h) {
    double k1, k2, k3, k4;

    k1 = h * derivative(Qt);
    k2 = h * derivative(Qt + 0.5 * h);
    k3 = h * derivative(Qt + 0.5 * h);
    k4 = h * derivative(Qt + h);

    return Qt + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

int main() {
    double qtur; //turbine flow m3/s
    double t = 0; //intial value of time for qtur
    double Z2; // tank wateer level above upstream reservoir
    double T = 0; // time
    double Qt = QTi; // Initial value of Qt
    double h = 0.1; // Increment value
    double threshold = 100; // Stop when Qt exceeds this value
    
    printf("T\tZ2\tQt\tqtur\n"); // printing the variable
    printf("(s)\t(m)\t(m3/s)\t(m3/s)\n"); //units of the variable
    
    for (; Qt < threshold; Qt += 0.1) {
        if (qtur<112)
        {
        qtur = QTi + (QTf - QTi) * t / TG;
        t = t + 0.1;
        }
        else qtur = 112;
        
        Z2 = Z1 + 0.1 * (Qt - qtur)/ AS;
        
        printf("%.3lf\t%.3lf\t%.3lf\t%.3lf\n", T, Z2, Qt, qtur);
        Qt = runge_kutta(Qt, h);
        T = T + 1;
    }
    
    return 0;
}