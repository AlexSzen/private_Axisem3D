// main.cpp
// created by Kuangdai on 26-Mar-2016 
// main 

#include "axisem.h"
#include <iostream>
#include <cstdio>
#include <ctime>

extern "C" void set_ftz();

int main(int argc, char *argv[]) {
    std::clock_t start;
    double duration;
    start = std::clock();

    // denormal float handling 
    set_ftz();
    
    // axisem main
    return axisem_main(argc, argv);

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"RUNTIME FROM MAIN: "<< duration <<'\n';    
}
