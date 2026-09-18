# fastSnark

Cryptographic library for creating non-interactive proofs of computations in C++.

## Features

* **Fast:** Creating of proof has **O(n)** time complexity.
* **Light:** Size of proof has **O(1)** memory complexity.
* **Custom:** Make custom operations to proof using basis R1CS.

## Dependecies

For building and working a library you should get:

* **C++17** compiler or greater (`g++`, `clang`)
* **CMake** version 3.14 or greater
* [nlohmann/json](https://github.com/nlohmann/json)
* [herumi/mcl](https://github.com/herumi/mcl)
* [OpenSSL (libcrypto)](https://openssl.org)

## Build and install

For building and installing release version of library you need:

```bash
git clone git@github.com:SIA1985/snrk.git
cd snrk

cmake -B build -S . -DCMAKE_BUILD_TYPE=Release

cmake --build build

sudo cmake --install build
```

## Quick start

After installation you can use library by including only one header: snrk/snrk.h.
In function snrk::init() the following parameters are set, each of them can be personally set after call of snrk::init():
* **MultiThreading** = true/false - turns on/off multicore processing 
* **SplinePartition** = 2.. - number of dots in spline segment (6 is optimal)

### Code example (`main.cpp`):

```cpp
#include <snrk/snrk.h>
#include <iostream>
#include <csignal>


void sigFpeHandler(int signum) {
    std::cout << "Making of proof error!" << std::endl;
    exit(signum);
}

void sigTermHandler(int signum) {
    std::cout << "Internal error!" << std::endl;
    exit(signum);
}

int main() {
    if (std::signal(SIGFPE, sigFpeHandler) == SIG_ERR) {
        exit(2);
    }

    if (std::signal(SIGTERM, sigTermHandler) == SIG_ERR) {
        exit(3);
    }


    snrk::init();

    /*Create custom operations*/
    enum Ops : int { 
   	    pls,
        prd,
    };

    snrk::Gates = {
    {pls, []GateLambdaHeader {
        return {{snrk::BaseGateType::Sum, {a, b}, c}};
    }
    },
    {prd, []GateLambdaHeader {
        return {{snrk::BaseGateType::Product, {a, b}, c}};
    }
    }};


    /*Create circuit and add computations*/
    snrk::Circuit c({/*inputs*/});

    auto v = snrk::CircuitValue(5.1);

    c.addGate(pls, {1, 2}, 3);
    c.addGate(prd, {2, 2}, 4);
    c.addGate(pls, {v, 1}, 6.1);
    c.addGate(prd, {2, v}, 10.2);


    /*Cryptographic prepare proof*/
    snrk::GlobalParams gp(42);

    auto pp = gp.PP();

    snrk::CircutParams cp(c, pp);


    /*Create and check proof*/
    snrk::ProverProof proof(cp, pp, "abc");

    if (proof.check(gp.VP(), "abc")) {
        std::cout << "Correct proof!" << std::endl;
    } else {
        std::cout << "Incorrect proof!" << std::endl;
    }
}
```

### Manual compilation via g++:
Compiling via bash:

```bash
g++ main.cpp -o app -lsnrk -lcrypto -lmcl
```

## Versions
0.1 Unsoundness.
0.2 Soundness. Memory complexity improved from **O(logN)** to **O(1)**.

## 📄 License

This project is distributed under the **MIT** license. Detailed information, as well as the copyrights of the authors of the third-party libraries used (OpenSSL, MCL, nlohmann/json), can be found in the file [LICENSE](LICENSE).

