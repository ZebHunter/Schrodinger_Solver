#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include "AbstractArray.h"

typedef std::complex<double> CMPLX;
typedef arrayt<CMPLX> arrayc;

//TODO:
// Заполнение массива приближений psi
// Заполнение массива координат
// Заполнения массива потенциальной энергии
// Проведение компресии массивов ( выбор элементов для построение
// Исправление main функции для проведения нормальных расчетов
// Графики!!!

// Приведенная постоянная Планка, Дж * с
const double h_bar = 1.05457266 * pow(10, -34);//1e-34

double L = 10.0;
double Lt = 5.0;
// cout << "Введите длину исследуемого промежутка: ";
//cin >> L;

// Количество точек аппроксимации
uint64_t N = 5000;
//cout << "Введите количество точек: ";
//cin >> N;

// Временной шаг
const uint64_t T = 20000;
//cout << "Введите количичество временных шагов: ";
//cin >> T;

// Масса, вводимая из консоли (если нужно, посторим окно ввода с форматом)
double m = 1;
//cout << "Введите массу: ";
//cin >> m;

// Коэффициент нелинейности в уравнении Шредингера
double g = 2.;
//cout << "Введите коэффициент нелинейности: ";
//cin >> g;

// Колебания
double w = 2.;
//cout << "Введите колебательную фазу: ";
//cin >> w;

// Начальное значение
double x0 = 0;

double alpha = 1.;
double delta = 0.5;
double sigma = 1.;

int64_t cutn = N / 100;
int64_t cutm = T / 100;

arrayc teta(N);
arrayc en(N);
arrayc V(N);

// Шаг по x
double dx = L / (N - 1);

// Шаг по t
double dt = L / (T - 1);

const CMPLX I(0.0, 1.0);

std::ofstream outfile_psi_re0("/psi-re0.dat");
std::ofstream outfile_psi_im0("/psi-im0.dat");
std::ofstream outfile_psi_sq0("/psi-sq0.dat");
std::ofstream outfile_psi_sq("/psi-sq25.dat");
std::ofstream outfile_v("/v.dat");
std::ofstream outfile_x("/x.dat");

arrayc psi(N);


double calc_V(double x) {
    if (x > x0) {
        return m * pow(w, 2) * pow(x, 2) / 2;
    }
    else return 0;
}


// Implements the implicit Crank-Nicolson method
void cn_method(arrayc& psi, arrayt<double>& V) {

    CMPLX b = 2.0 - I * pow(dx, 2) / (alpha * sigma * dt) - (1 - delta) * (g * pow(dx, 2) / (alpha * sigma) * pow(psi(2), 2) + calc_V(dx));
    CMPLX d = (sigma - 1.0) / sigma * (psi(1) + psi(3)) + psi(2) * (I * pow(dx, 2) / (dt * sigma * alpha) +
        2.0 * (1.0 - sigma) / sigma - delta * pow(dx, 2) / (alpha * sigma) * g * (psi(2)));

    teta(3) = 1.0 / b;
    en(3) = -d / b;

    for (size_t i = 3; i < N - 2; i++) {
        b = 2.0 - I * pow(dx, 2) / (alpha * sigma * dt) - (1 - delta) * (g * pow(dx, 2) / (alpha * sigma) * pow(psi(i), 2) + calc_V(dx));
        d = (sigma - 1.0) / sigma * (psi(i - 1) + psi(i + 1)) + psi(i) * (I * pow(dx, 2) / (dt * sigma * alpha) +
            2.0 * (1.0 - sigma) / sigma - delta * pow(dx, 2) / (alpha * sigma) * g * (psi(i)));

        teta(i + 1) = 1.0 / (b - teta(i));
        en(i + 1) = (en(i) - d) / (b - teta(i));
    }

    b = 2.0 - I * pow(dx, 2) / (alpha * sigma * dt) - (1 - delta) * (g * pow(dx, 2) / (alpha * sigma) * pow(psi(N), 2) + calc_V(dx));
    d = (sigma - 1.0) / sigma * (psi(N - 2) + psi(N)) + psi(2) * (I * pow(dx, 2) / (dt * sigma * alpha) +
        2.0 * (1.0 - sigma) / sigma - delta * pow(dx, 2) / (alpha * sigma) * g * (psi(N - 1)));

    en(N - 1) = -d / b;
    psi(N - 1) = en(N - 1);

    for (size_t i = 2; i < N - 2; i++) {
        psi(N - 1) = teta(N - i + 1) + en(N - i + 1);
    }
}

int main() {

    // define arrays for wavefunction, potential
    arrayc psi(N);
    arrayt<double> V(N);

    // initialize psi and V, populate x values
    for (int j = 0; j < N; j++) {

        V(j) = calc_V(j * dx);
        outfile_psi_re0 << psi(j).real() << std::endl;
        outfile_psi_im0 << psi(j).imag() << std::endl;
        outfile_psi_sq0 << abs((psi(j) * psi(j)).real()) << std::endl;
        outfile_v << V(j) << std::endl;
    }

    cn_method(psi, V);

    // store final values in outfiles
    for (int j = 0; j < N; j++) {
        outfile_psi_sq << abs((psi(j) * psi(j)).real()) << std::endl;
        outfile_x << j * dx << std::endl;
    }
    outfile_psi_re0 << std::endl; 
    outfile_psi_im0 << std::endl;
    outfile_psi_sq0 << std::endl;
    outfile_psi_sq << std::endl;
    outfile_v << std::endl;
    outfile_x << std::endl;

    outfile_psi_re0.close();
    outfile_psi_im0.close();
    outfile_psi_sq0.close();
    outfile_psi_sq.close();
    outfile_v.close();
    outfile_x.close();

    return 0;

} // end main