#include "../tools/getArgs.hpp"
#include "../BCE/BCE.h"
#include <string>
#include <chrono>
#include <iostream>

int main(int argc, char* argv[]) {
    argsController ac(argc, argv);

    if(!ac.exist("-f")) {
        printf("file path: -f\n");
        return 0;
    }

    if(!ac.exist("-d")) {
        printf("graph order: -d(core, two)\n");
        return 0;
    }

    uint32_t ls = 1, rs = 1;
    if(ac.exist("-l")) ls = std::stoi(ac["-l"]);
    if(ac.exist("-r")) rs = std::stoi(ac["-r"]);

    std::string fPath = ac["-f"];
    std::string order = ac["-d"];

    std::cout << fPath << ' ' << order << std::endl;
    std::cout << "l " << ls << std::endl;
    std::cout << "r " << rs << std::endl;

    if(ac.exist("noUVM")) {
        BCE bce(fPath, 0, order, ls, rs);

        auto t1 = std::chrono::steady_clock::now();
            
        bce.run();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    }
    else {
        BCE bce(fPath, 1, order, ls, rs);

        auto t1 = std::chrono::steady_clock::now();
            
        bce.run();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    }
    

    return 0;
}