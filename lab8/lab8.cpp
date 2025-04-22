#include <iostream>
#include <cmath>
#include <numbers>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>

using namespace std;

double exact_derivative(double x)
{
    return -sin(x);
}

long double long_exact_derivative(long double x)
{
    return -sinl(x);
}

template <typename T>
T forward_diff(T x, T h, T (*f)(T))
{
    return (f(x + h) - f(x)) / h;
}

template <typename T>
T backward_diff(T x, T h, T (*f)(T))
{
    return (f(x) - f(x - h)) / h;
}

template <typename T>
T central_diff(T x, T h, T (*f)(T))
{
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

template <typename T>
T three_point_forward(T x, T h, T (*f)(T))
{
    return (-3.0 / 2.0 * f(x) + 2 * f(x + h) - 1.0 / 2.0 * f(x + 2.0 * h)) / h;
}

template <typename T>
T three_point_backward(T x, T h, T (*f)(T))
{
    return (1.0 / 2.0 * f(x - 2.0 * h) - 2.0 * f(x - h) + 3.0 / 2.0 * f(x)) / h;
}

double obliczSrednia(const std::vector<double>& vec)
{
    double suma = 0;
    for (double liczba : vec)
    {
        suma += liczba;
    }
    return suma / vec.size();
}


int main()
{
    ofstream file1("double_errors.txt");
    ofstream file2("long_double_errors.txt");
    if (!file1.is_open() || !file2.is_open())
    {
        cerr << "Nie mozna otworzyc pliku do zapisu!" << endl;
        return 1;
    }
    int num_steps = 10000;
    vector<double> h_values, errors[7];
    for (int i = 0; i < num_steps; i++)
    {
        double x_left = 0.0;
        double x_middle = numbers::pi / 4.0;
        double x_right = numbers::pi / 2.0;


        double h_min = 10e-16;
        double h_max = 10e-1;

        // Logarytmiczny rozkład h między h_min a h_max
        double h = h_min * pow(10.0, i * log10(h_max / h_min) / (num_steps - 1));

        // punkt początkowy
        double forward_left = forward_diff(x_left, h, cos);
        double forward_left_error = abs(exact_derivative(x_left) - forward_left);

        double three_point_forward_left = three_point_forward(x_left, h, cos);
        double three_point_forward_left_error = abs(exact_derivative(x_left) - three_point_forward_left);

        // punkt centralny
        double forward_middle = forward_diff(x_middle, h, cos);
        double forward_middle_error = abs(exact_derivative(x_middle) - forward_middle);

        double backward_middle = backward_diff(x_middle, h, cos);
        double backward_middle_error = abs(exact_derivative(x_middle) - backward_middle);

        double central_middle = central_diff(x_middle, h, cos);
        double central_middle_error = abs(exact_derivative(x_middle) - central_middle);

        // punkt końcowy
        double backward_right = backward_diff(x_right, h, cos);
        double backward_right_error = abs(exact_derivative(x_right) - backward_right);

        double three_point_backward_right = three_point_backward(x_right, h, cos);
        double three_point_backward_right_error = abs(exact_derivative(x_right) - three_point_backward_right);

        h_values.push_back(h);
        errors[0].push_back(forward_left_error);
        errors[1].push_back(three_point_forward_left_error);
        errors[2].push_back(forward_middle_error);
        errors[3].push_back(backward_middle_error);
        errors[4].push_back(central_middle_error);
        errors[5].push_back(backward_right_error);
        errors[6].push_back(three_point_backward_right_error);

        file1 << log10(h) << " " << log10(forward_left_error) << " " << log10(three_point_forward_left_error) << " " <<
            log10(forward_middle_error) << " " << log10(backward_middle_error) << " " << log10(central_middle_error) <<
            " " << log10(backward_right_error) << " " << log10(three_point_backward_right_error) << endl;
    }
    cout<< "Double" << endl;

    double forward_left_error_p = log10(errors[0][9999]/errors[0][8000]) / log10(h_values[9999]/h_values[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy dwupunktowa punktu poczatkowego: " << forward_left_error_p << endl;

    double three_point_forward_left_error_p = log10(errors[1][9999]/errors[1][9000]) / log10(h_values[9999]/h_values[9000]);
    cout << "Rzad dokladnosci przyblizenia roznicy trzypunktowej punktu poczatkowego: " << three_point_forward_left_error_p << endl;


    double forward_middle_error_p = log10(errors[2][9999]/errors[2][8000]) / log10(h_values[9999]/h_values[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy progresywnej punktu centralnego: " << forward_middle_error_p << endl;

    double backward_middle_error_p = log10(errors[3][9999]/errors[3][8000]) / log10(h_values[9999]/h_values[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy wstecznej punktu centralnego: " << backward_middle_error_p << endl;

    double central_middle_error_p = log10(errors[4][9999]/errors[4][8000]) / log10(h_values[9999]/h_values[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy centralnej punktu centralnego: " << central_middle_error_p << endl;

    double backward_right_error_p = log10(errors[5][9999]/errors[5][8000]) / log10(h_values[9999]/h_values[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy dwupunktowa punktu koncowego: " << backward_right_error_p << endl;

    double three_point_backward_right_error_p = log10(errors[6][9999]/errors[6][8000]) / log10(h_values[9999]/h_values[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy trzypunktowej punktu koncowego: " << three_point_backward_right_error_p << endl << endl;

    vector<long double> h_values_long, errors_long[7];
    for (int i = 0; i < num_steps; i++)
    {
        long double x_left = 0.0;
        long double x_middle = numbers::pi / 4.0;
        long double x_right = numbers::pi / 2.0;


        long double h_min = 10e-19;
        long double h_max = 10e-1;

        // Logarytmiczny rozkład h między h_min a h_max
        long double h = h_min * pow(10.0, i * log10l(h_max / h_min) / (num_steps - 1));

        // punkt początkowy
        long double forward_left = forward_diff(x_left, h, cosl);
        long double forward_left_error = abs(long_exact_derivative(x_left) - forward_left);

        long double three_point_forward_left = three_point_forward(x_left, h, cosl);
        long double three_point_forward_left_error = abs(long_exact_derivative(x_left) - three_point_forward_left);

        // punkt centralny
        long double forward_middle = forward_diff(x_middle, h, cosl);
        long double forward_middle_error = abs(long_exact_derivative(x_middle) - forward_middle);

        long double backward_middle = backward_diff(x_middle, h, cosl);
        long double backward_middle_error = abs(long_exact_derivative(x_middle) - backward_middle);

        long double central_middle = central_diff(x_middle, h, cosl);
        long double central_middle_error = abs(long_exact_derivative(x_middle) - central_middle);

        // punkt końcowy
        long double backward_right = backward_diff(x_right, h, cosl);
        long double backward_right_error = abs(long_exact_derivative(x_right) - backward_right);

        long double three_point_backward_right = three_point_backward(x_right, h, cosl);
        long double three_point_backward_right_error = abs(long_exact_derivative(x_right) - three_point_backward_right);

        h_values_long.push_back(h);
        errors_long[0].push_back(forward_left_error);
        errors_long[1].push_back(three_point_forward_left_error);
        errors_long[2].push_back(forward_middle_error);
        errors_long[3].push_back(backward_middle_error);
        errors_long[4].push_back(central_middle_error);
        errors_long[5].push_back(backward_right_error);
        errors_long[6].push_back(three_point_backward_right_error);

        // cout << scientific <<  "h = " << h << endl
        //         << scientific << "forward = " << forward_left << " error = " << forward_left_error << endl
        //         << scientific <<  "backward = " << backward_right << " error = " << backward_right_error << endl
        //         << scientific << "central = " << central_middle << " error = " << central_middle_error << endl
        //         << scientific <<  "three_point_forward = " << three_point_forward_left << " error = " << three_point_forward_left_error << endl
        //         << scientific <<  "three_point_backward = " << three_point_backward_right << " error = " << three_point_backward_right_error  << endl << endl;

        file2 << log10l(h) << " " << log10l(forward_left_error) << " " << log10l(three_point_forward_left_error) << " "
            << log10l(forward_middle_error) << " " << log10l(backward_middle_error) << " " <<
            log10l(central_middle_error) << " " << log10l(backward_right_error) << " " << log10l(
                three_point_backward_right_error) << endl;
    }
    cout << "Long double" << endl;
    long double forward_left_error_p_long = log10l(errors_long[0][9999]/errors_long[0][8000]) / log10l(h_values_long[9999]/h_values_long[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy dwupunktowa punktu poczatkowego: " << forward_left_error_p_long << endl;

    long double three_point_forward_left_error_p_long = log10l(errors_long[1][9999]/errors_long[1][9000]) / log10l(h_values_long[9999]/h_values_long[9000]);
    cout << "Rzad dokladnosci przyblizenia roznicy trzypunktowej punktu poczatkowego: " << three_point_forward_left_error_p_long << endl;


    long double forward_middle_error_p_long = log10l(errors_long[2][9999]/errors_long[2][8000]) / log10l(h_values_long[9999]/h_values_long[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy progresywnej punktu centralnego: " << forward_middle_error_p_long << endl;

    long double backward_middle_error_p_long = log10l(errors_long[3][9999]/errors_long[3][8000]) / log10l(h_values_long[9999]/h_values_long[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy wstecznej punktu centralnego: " << backward_middle_error_p_long << endl;

    long double central_middle_error_p_long = log10l(errors_long[4][9999]/errors_long[4][8000]) / log10l(h_values_long[9999]/h_values_long[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy centralnej punktu centralnego: " << central_middle_error_p_long << endl;

    long double backward_right_error_p_long = log10l(errors_long[5][9999]/errors_long[5][8000]) / log10l(h_values_long[9999]/h_values_long[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy dwupunktowa punktu koncowego: " << backward_right_error_p_long<< endl;

    long double three_point_backward_right_error_p_long = log10l(errors_long[6][9999]/errors_long[6][8000]) / log10l(h_values_long[9999]/h_values_long[8000]);
    cout << "Rzad dokladnosci przyblizenia roznicy trzypunktowej punktu koncowego: " << three_point_backward_right_error_p_long << endl;
    return 0;
}
