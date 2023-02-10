#include <iostream>
#include <random>

using namespace std;


template<typename T>
T randunif(T range_from, T range_to) {
    std::random_device                  rand_dev;                       // stock random number libraries
    std::mt19937                        generator(rand_dev());          // use the mersenne twister as the underlying RN generator
    std::uniform_real_distribution<T>   distr(range_from, range_to);
    return distr(generator);
}

template<typename T>
T randnorm(T mean, T sd) {
    std::random_device                  rand_dev;                       // stock random number libraries
    std::mt19937                        generator(rand_dev());          // use the mersenne twister as the underlying RN generator
    std::normal_distribution<T>         distr(mean, sd);
    return distr(generator);
}
