#include "../tools/getArgs.hpp"
// #include "../plex/biplex.h"
#include "../plex/biplexv2.h"
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

    int graphMode = (ac.exist("noUVM") == false);

    int k = 2;
    if(ac.exist("-k")) k = atoi(ac["-k"].c_str());

    uint64_t outPutT = INT64_MAX;
    if(ac.exist("-outPutT")) outPutT = atoi(ac["-outPutT"].c_str());

    std::cout << "file path: " << fPath << " "
              << "k: " << k << " "
              << "order: " << order << " "
              << "l : " << ls << " "
              << "r : " << rs << " "
              << "mode: " << graphMode << std::endl;

    // if(ac.exist("-v2")) {
        std::cout << "-v2" << std::endl;

        biplexv2 bcev2(fPath, graphMode, k, order, ls, rs, outPutT);
 
        auto t1 = std::chrono::steady_clock::now();
            
        bcev2.run();

        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "time:" << duration.count() << "ms" << std::endl;
    // }
    // else {
    //     biplex bce(fPath, graphMode, k);
 
    //     auto t1 = std::chrono::steady_clock::now();
            
    //     bce.run();

    //     auto t2 = std::chrono::steady_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    //     std::cout << "time:" << duration.count() << "ms" << std::endl;
    // }

    


    return 0;
}