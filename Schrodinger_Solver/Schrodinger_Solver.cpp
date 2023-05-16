#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <GL/glut.h>
#include "AbstractArray.h"

typedef complex<double> CMPLX;
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


arrayc psi(N / cutn, T / cutm);

arrayt<double> X(N);

void rideX(arrayt<double>& x) {
    for (size_t i = 0; i < N; i++) {
        if (i == 0) x(i) = x0; 
        else x(i) = x0 + i * dx;
    }
}

void function(arrayt<double>& x, arrayc& psi) {
    for (size_t i = 0; i < N; i++) {
        psi(i) = sin(x(i));
    }
}

void compression(arrayc& begin, int k) {
    if (k % cutm == 0) {
        for (size_t i = 0; i < N; i++) {
            if (k % cutn == 0) {
                psi(i / cutn, k / cutm) = begin(i);
            }
        }
    }
}

double calc_V(double x) {
    if (x > x0) {
        return m * pow(w, 2) * pow(x, 2) / 2;
    }
    else return 0;
}

arrayc psi1(N);

// Implements the implicit Crank-Nicolson method
void cn_method(arrayc psi) {
    //for (size_t j = 0; j < T / cutm; j++) {
        psi(1) = 0;
        psi(N) = 0;
        CMPLX b = 2.0 - I * pow(dx, 2) / (alpha * sigma * dt) - (1 - delta) * (g * pow(dx, 2) / (alpha * sigma) * pow((psi(2)), 2));
        CMPLX d = (sigma - 1.0) / sigma * (psi(1) + psi(3)) + psi(2) * (I * pow(dx, 2) / (dt * sigma * alpha) +
            2.0 * (1.0 - sigma) / sigma - delta * pow(dx, 2) / (alpha * sigma) * g * pow((psi(2)), 2) - calc_V(dx) * pow(dx, 2) / (alpha * sigma));

        teta(3) = 1.0 / b;
        en(3) = -d / b;

        for (size_t i = 3; i < N - 2; i++) {
            b = 2.0 - I * pow(dx, 2) / (alpha * sigma * dt) - (1 - delta) * (g * pow(dx, 2) / (alpha * sigma) * pow((psi(i)), 2));
            d = (sigma - 1.0) / sigma * (psi(i - 1) + psi(i + 1)) + psi(i) * (I * pow(dx, 2) / (dt * sigma * alpha) +
                2.0 * (1.0 - sigma) / sigma - delta * pow(dx, 2) / (alpha * sigma) * g * pow((psi(i)), 2) - calc_V(dx) * pow(dx, 2) / (alpha * sigma));

            teta(i + 1) = 1.0 / (b - teta(i));
            en(i + 1) = (en(i) - d) / (b - teta(i));
        }

        b = 2.0 - I * pow(dx, 2) / (alpha * sigma * dt) - (1 - delta) * (g * pow(dx, 2) / (alpha * sigma) * pow(psi(N), 2));
        d = (sigma - 1.0) / sigma * (psi(N - 2) + psi(N)) + psi(N - 1) * (I * pow(dx, 2) / (dt * sigma * alpha) +
            2.0 * (1.0 - sigma) / sigma - delta * pow(dx, 2) / (alpha * sigma) * g * pow((psi(N - 1)), 2) - calc_V(dx) * pow(dx, 2) / (alpha * sigma));

        en(N - 1) = -d / b;
        psi(N - 1) = en(N - 1);

        for (size_t i = 2; i < N - 2; i++) {
            psi(N - i) = teta(N - i + 1) + en(N - i + 1);
        }
   // }
}

//const int X_COORD = 50;// X - размерность ] должны
//const int Y_COORD = 50;// Y - размерность ] быть равными
const float ITERATIONS = 0.00005;// прорисовка графика (чем меньше тем лучше)

int x_off = L;// начало
int y_off = Lt;// оси координат

//исходная функция


void drawgrid(float SERIF_OFFSET, float SERIF_DISTANCE) {
    glBegin(GL_LINES);
    //задаем цвета
    glColor3f(0.0, 0.0, 0.0);

    //рисуем координатные оси
    //горизонталь
    glVertex2f(0.0, Lt);
    glVertex2f(2*L, Lt);
    //засечки по горизонтали
    int p = Lt;
    for (int i = L; i < 2*L; i += SERIF_DISTANCE, p -= SERIF_DISTANCE) {
        glVertex2f(i, Lt);
        glVertex2f(i, Lt + SERIF_OFFSET);

        glVertex2f(p, Lt);
        glVertex2f(p, Lt + SERIF_OFFSET);
    }
    //вертикаль
    int t = Lt;
    glVertex2f(L, 2*Lt);
    glVertex2f(L, 0.0);
    //засечки по вертикали
    for (int i = Lt; i < 2*Lt; i += SERIF_DISTANCE, t -= SERIF_DISTANCE) {
        glVertex2f(L, i);
        glVertex2f(Lt + SERIF_OFFSET, i);

        glVertex2f(L, t);
        glVertex2f(Lt + SERIF_OFFSET, t);
    }
    glEnd();
}

void drawfunc() {
    //рисуем график
    glBegin(GL_POINTS);
    double j = 0;
    glColor3f(0.8, 0.0, 0.8);
    for (int i = 0; i < L * 2; i++) {
        //перерасчитываем координаты
        j = psi1(i).real();
        glVertex2f(x_off + X(i), y_off + j);//не убирать x и y!! это оффсет по осям!
    }
    glEnd();
}



void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    cout << "Osnovnie toshki po vashemu grafiku: \n";

    drawgrid(0.3, 5);
    drawfunc();

    glutSwapBuffers();

    glFlush();
}

int main(int argc,char **argv) {
    int k = 2;

    rideX(X);

    
    function(X, psi1);
    

    /*for (size_t i = 0; i < T; i++) {
        cn_method(psi1);
        compression(psi1, k);
        
    }
    */
    

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutInitWindowPosition(500, 200);
    glutCreateWindow("GLUT_TESTING_APP");

    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //пространство координат
    glOrtho(0.0, L*2, 0.0, Lt * 2, -1.0, 1.0);

    glutDisplayFunc(display);
    glutMainLoop();

    return 0;
} // end main