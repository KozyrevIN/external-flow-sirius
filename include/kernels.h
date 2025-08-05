#include <cmath>

double kernel_base(double r){
    return exp(-r * r)/M_PI;
}
double kernel_2(double r){
    return kernel_base(r) * (-2);
}
double kernel_4(double r){
    return kernel_base(r) * (-6+2*std::pow(r, 2));
}
double kernel_6(double r){
    return kernel_base(r) * (-12+8*std::pow(r, 2)-std::pow(r, 4));
}
double kernel_8(double r){
    return kernel_base(r) * (-20+16*std::pow(r, 2)-4*std::pow(r, 4)+std::pow(r, 6));
}

unordered_map<int, decltype(kernel_base)> kernel_map = {
    {2, kernel_2},
    {4, kernel_4},
    {6, kernel_6},
    {8, kernel_8},
};