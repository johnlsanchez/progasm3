#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_set>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random> 
#include <cmath>

//#define _USE_MATH_DEFINES
#define MAX_ITER 25000

std::mt19937 gen(time(nullptr));

std::vector<uint64_t> merge_sort(std::vector<uint64_t> nums) {
    if (nums.size() == 1) {
        return nums;
    }

    std::vector<uint64_t> front = std::vector<uint64_t>(nums.begin(), nums.begin() + (nums.size() / 2));
    std::vector<uint64_t> back = std::vector<uint64_t>(nums.begin() + (nums.size() / 2), nums.end());

    std::vector<uint64_t> front_sort = merge_sort(front);
    std::vector<uint64_t> back_sort = merge_sort(back);

    std::vector<uint64_t> merged;
    while (front_sort.size() != 0 || back_sort.size() != 0) {
        if (front_sort.size() == 0) {
            merged.push_back(back_sort.at(0));
            back_sort.erase(back_sort.begin());
        } else if (back_sort.size() == 0) {
            merged.push_back(front_sort.at(0));
            front_sort.erase(front_sort.begin());
        } else if (front_sort.at(0) < back_sort.at(0)) {
            merged.push_back(front_sort.at(0));
            front_sort.erase(front_sort.begin());
        } else {
            merged.push_back(back_sort.at(0));
            back_sort.erase(back_sort.begin());
        }
    }
    return merged;
}

uint32_t binary_insert(std::vector<uint64_t> nums, uint64_t val, uint32_t lh, uint32_t rh) {
    if (rh <= lh) {
        //printf("here1\n");
        return (val > nums.at(lh)) ? (lh + 1) : lh;
    }
    if (rh == lh + 1) {
        //printf("here2\n");
        return (nums.at(lh) > val) ? lh : rh;
    }
    uint32_t mid = (lh + rh) / 2;
    if (val == nums.at(mid)) {
        //printf("here3\n");
        return mid;
    }
    if (val > nums.at(mid)) {
        //printf("%i %i %i\n", val, mid, rh);
        return binary_insert(nums, val, mid, rh);
    }
    //printf("here\n");
    return binary_insert(nums, val, lh, mid);
}

uint64_t karmarkar_karp(std::vector<uint64_t> nums) {
    // sort nums -> O(nlog(n))
    std::vector<uint64_t> sorted = merge_sort(nums);

    // for each n, remove largest 2 and insert difference using binary heap
    // O(nlog(n)) if using binary heap
    while (sorted.size() > 2) {
        // pull largest two
        // TODO CHANGE
        uint64_t t1 = sorted.at(sorted.size() - 1);
        uint64_t t2 = sorted.at(sorted.size() - 2);
        sorted.pop_back();
        sorted.pop_back();

        // find difference
        t1 -= t2;

        // insert into binary heap
        // TODO CHANGE
        uint32_t n = binary_insert(sorted, t1, 0, sorted.size());
        if (n != sorted.size()) {
            auto it = sorted.begin() + n;
            sorted.insert(it, t1);
        } else {
            sorted.push_back(t1);
        }
    }
    // return difference of last two values
    return sorted.at(1) - sorted.at(0);
}

// print vector of booleans as 1s and 0s
void print(std::vector<bool> a) {
    for (uint32_t i = 0; i < a.size(); i++) {
        printf("%i ", a[i] ? 1 : 0);
    }
    printf("\n");
}

// Calculate residue for non-prepartition algorithms
// sum over elements of nums[i] * S[i]
uint64_t residue(std::vector<uint64_t> nums, std::vector<bool> S) {
    // declare running sum of value for both groups
    uint64_t sums = 0;
    uint64_t suml = 0;

    for (uint32_t i = 0; i < nums.size(); ++i) {
        if (S.at(i)) {
            sums += nums.at(i); // if positive, add to sums
        } else {
            suml += nums.at(i); // if negative, add to suml
        }
    }
    // return difference with respect to unsigned largest
    return (sums > suml) ? (sums - suml) : (suml - sums);
}

std::vector<bool> generate_random_sol(uint32_t size) {
    // randomly generate a list of size size
    std::vector<bool> S;
    for (uint32_t i = 0; i < size; ++i) {
        S.push_back((gen() % 2 == 0));
    }
    return S;
}

uint64_t repeated_random_signs(std::vector<uint64_t> nums) {
    const uint32_t size = nums.size();

    // generate random solution S
    std::vector<bool> S = generate_random_sol(size);
    uint64_t current_residue = residue(nums, S);
    
    for (uint32_t i = 0; i < MAX_ITER; i++)
    {
        // generate random solution S_prime
        std::vector<bool> S_prime = generate_random_sol(size);
        uint64_t prime_residue = residue(nums, S_prime);

        // compare S and S_prime, if S_prime better -> S = S_prime
        if (prime_residue < current_residue)
        {
            S = S_prime;
            current_residue = prime_residue;
        }
    }
    return current_residue;
}

uint64_t hill_climbing_signs(std::vector<uint64_t> nums) {
    const uint32_t size = nums.size();

    // generate random solution S
    std::vector<bool> S = generate_random_sol(size);
    uint64_t current_residue = residue(nums, S);

    for (uint32_t i = 0; i < MAX_ITER; i++) {
        // generate random a, b s.t. a != b
        int a = gen() % nums.size();
        int b = gen() % nums.size();
        while (a == b) {
            b = gen() % nums.size();
        }

        // flip S in one place, with 1/2 prob of another
        std::vector<bool> S_prime = S;
        S_prime.at(a) = !S_prime.at(a);
        S_prime.at(b) = (gen() % 2 == 0) ? S_prime.at(b) : !S_prime.at(b);
        uint64_t prime_residue = residue(nums, S_prime);

        // compare, if S_prime better than S -> S = S_prime
        if (prime_residue < current_residue)
        {
            S = S_prime;
            current_residue = prime_residue;
        }
    }
    return current_residue;
}

double cooling_schedule(uint32_t i) {
    return (pow(10, 10) * pow(0.8, i / 300));
}

uint64_t simulated_annealing_signs(std::vector<uint64_t> nums) {
    // random distribution
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    const uint32_t size = nums.size();

    // start with random solutions
    std::vector<bool> S = generate_random_sol(size), S_dp = S;
    uint64_t current_residue = residue(nums, S);
    uint64_t dp_residue = current_residue;

    for (uint32_t i = 0; i < MAX_ITER; i++) {
        // generate random a, b s.t. a != b
        int a = gen() % nums.size();
        int b = gen() % nums.size();
        while (a == b) {
            b = gen() % nums.size();
        }

        // random flip at a, with 1/2 of flipping at b
        std::vector<bool> S_prime = S;
        S_prime.at(a) = !S_prime.at(a);
        S_prime.at(b) = (gen() % 2 == 0) ? S_prime.at(b) : !S_prime.at(b);
        uint64_t prime_residue = residue(nums, S_prime);

        // if better, set S = S_prime
        if (prime_residue < current_residue) {
            S = S_prime;
            current_residue = prime_residue;
        } else {
            // else consider flipping at random according to cooling schedule and difference in residue
            double fn = std::exp(-(((double)(prime_residue - current_residue)) / (cooling_schedule(i))));
            if (dis(gen) < fn) {
                S = S_prime;
                current_residue = prime_residue;
            }
        }
        // if S better than S_dp -> S_dp = S
        if (current_residue < dp_residue) {
            S_dp = S;
            dp_residue = current_residue;
        }
    }
    // return best solution so far
    return dp_residue;
}

std::vector<uint8_t> generate_prepartition(uint32_t size) {
    // randomly generate a list of size size
    std::vector<uint8_t> P;
    for (uint32_t i = 0; i < size; ++i) {
        P.push_back((uint8_t)(gen() % size));
    }
    return P;
}

std::vector<uint64_t> calculate_aprime(std::vector<uint64_t> nums, std::vector<uint8_t> s_nums) {
    // initialize target vector
    std::vector<uint64_t> a_prime(nums.size(), 0);

    // populate with values such that p[i] = sum over nums where num[j] = i
    for (uint32_t i = 0; i < nums.size(); i++) {
        a_prime.at(s_nums.at(i)) += nums.at(i);
    }

    // remove zeros
    std::vector<uint64_t> real_a_prime;
    for (uint32_t i = 0; i < nums.size(); i++) {
        uint64_t val = a_prime.at(i);
        if (val) {
            real_a_prime.push_back(val);
        }
    }

    return a_prime;
}

uint64_t repeated_random_prepartition(std::vector<uint64_t> nums) {
    // start with random solution P
    std::vector<uint8_t> P = generate_prepartition(nums.size());
    uint64_t current_residue = karmarkar_karp(calculate_aprime(nums, P));

    // iterate from [0, MAX_ITER)
    for (uint32_t i = 0; i < MAX_ITER; i++)
    {
        // generate a P_prime, completely at random
        std::vector<uint8_t> P_prime = generate_prepartition(nums.size());
        uint64_t prime_residue = karmarkar_karp(calculate_aprime(nums, P_prime));

        // if P_prime residue better than P, swap
        if (prime_residue < current_residue)
        {
            P = P_prime;
            current_residue = prime_residue;
        }
    }
    return current_residue;
}

uint64_t hill_climbing_prepartition(std::vector<uint64_t> nums) {
    // start with random solution P
    std::vector<uint8_t> P = generate_prepartition(nums.size());
    const int size = nums.size();
    uint64_t current_residue = karmarkar_karp(calculate_aprime(nums, P));

    // iterate from [0, MAX_ITER)
    for (uint32_t i = 0; i < MAX_ITER; i++)
    {
        // randomly generate a, b in [0, size) such that P[a] != b
        int a = gen() % size;
        int b = gen() % size;
        while (P.at(a) == b) {
            b = gen() % size;
        }

        // generate a P_prime, random neighbor of P
        std::vector<uint8_t> P_prime = P;
        P_prime.at(a) = b;
        uint64_t prime_residue = karmarkar_karp(calculate_aprime(nums, P_prime));

        // if P_prime residue better than P, swap
        if (prime_residue < current_residue)
        {
            P = P_prime;
            current_residue = prime_residue;
        }
    }
    return current_residue;
}

uint64_t simulated_annealing_prepartition(std::vector<uint64_t> nums) {
    const int size = nums.size();

    // random distribution 
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // start with random solution P set to P_dp
    std::vector<uint8_t> P, P_dp;
    P = generate_prepartition(size);
    uint64_t current_residue = karmarkar_karp(calculate_aprime(nums, P));
    uint64_t dp_residue = current_residue;
    P_dp = P;

    // Iterate MAX_ITER times
    for (uint32_t i = 0; i < MAX_ITER; i++) {
        // find random a, b s.t. a != b
        int a = gen() % size;
        int b = gen() % size;
        while (a == b) {
            b = gen() % size;
        }

        // randomly change P_p[a] = b
        std::vector<uint8_t> P_prime = P;
        P_prime.at(a) = b;
        uint64_t prime_residue = karmarkar_karp(calculate_aprime(nums, P_prime));

        // if better solution change P
        if (prime_residue < current_residue) {
            P = P_prime;
            current_residue = prime_residue;
        } else {
            // if not, consider randomly pursuing it anyways
            double fn = std::exp(-(((double)(prime_residue - current_residue)) / (cooling_schedule(i))));
            if (dis(gen) < fn) {
                P = P_prime;
                current_residue = prime_residue;
            }
        }
        // if current solution is better, then set P_dp = P, P_dp:best so far
        if (current_residue < dp_residue) {
            P_dp = P;
            dp_residue = current_residue;
        }
    }
    return dp_residue;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: ./partition 0 algorithm inputfile\n");
    }

    // argv[] -> variables
    uint8_t debug = std::stoi(argv[1]);
    uint8_t algo = std::stoi(argv[2]);
    std::vector<uint64_t> input;

    // I/O
    std::string line;
    std::ifstream file(argv[3]);
    if (file.is_open()){
        while (getline(file, line)) {
            input.push_back(std::stoll(line));
        }
        file.close();
    }

    // debug - run all tests once
    if (debug) {
        printf("Karmarkar-Karp: %llu\n", karmarkar_karp(input));                                    // :)
        printf("Repeated Random: %llu\n", repeated_random_signs(input));                            // :)
        printf("Hill-Climbing: %llu\n", hill_climbing_signs(input));                                // :)
        printf("Annealing: %llu\n", simulated_annealing_signs(input));                              // :)
        printf("Prepartition Repeated Random: %llu\n", repeated_random_prepartition(input));        // :)
        printf("Prepartition Hill-Climbing: %llu\n", hill_climbing_prepartition(input));            // :)
        printf("Prepartition Simulated Annealing: %llu\n", simulated_annealing_prepartition(input));// :D
        return 0;
    }

    // else run specified algorithm according to:
    // 0 Karmarkar-Karp
    // 1 Repeated Random
    // 2 Hill Climbing
    // 3 Simulated Annealing
    // 11 Prepartitioned Repeated Random
    // 12 Prepartitioned Hill Climbing
    // 13 Prepartitioned Simulated Annealing
    if (algo == 0) {
        printf("%llu\n", karmarkar_karp(input));       
    } else if (algo == 1) {
        printf("%llu\n", repeated_random_signs(input));
    } else if (algo == 2) {
        printf("%llu\n", hill_climbing_signs(input));
    } else if (algo == 3) {
        printf("%llu\n", simulated_annealing_signs(input));
    } else if (algo == 11) {
        printf("%llu\n", repeated_random_prepartition(input));
    } else if (algo == 12) {
        printf("%llu\n", hill_climbing_prepartition(input));
    } else if (algo == 13) {
        printf("%llu\n", simulated_annealing_prepartition(input));
    }
    return 0;
}