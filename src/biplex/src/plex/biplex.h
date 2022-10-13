#ifndef BIPLEX_H
#define BIPLEX_H

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"
#include <string>
#include <vector>

class biplex {
private:
    biGraph * g;
    int k;
    uint64_t totalMaximalBiPlex = 0;

    LinearSet C[2], X[2], R[2];

    std::vector<std::vector<uint32_t>> ws;//buffer for nonNeighbors in cand of pivot

    struct parameter{
        uint32_t value[2];
        // parameter(parameter & t) {
        //     value[0] = t[0];
        //     value[1] = t[1];
        // }

        uint32_t & operator [] (uint32_t i) { return value[i]; }
    };

    struct nonNeiMatainer{
        std::vector<uint32_t> buffer[2];
        std::vector<uint32_t> cntNonNei[2];
        // std::vector<std::vector<uint32_t>> cand[2];
        uint32_t n[2];
        int k;

        void init(uint32_t n1, uint32_t n2, int k) {
            this->n[0] = n1; this->n[1] = n2; 
            this->k = k;
            buffer[0].resize(n1 * k);
            buffer[1].resize(n2 * k);
            cntNonNei[0].resize(n1 + 1);
            cntNonNei[1].resize(n2 + 1);
        }
        
        //u is a t-side vertex
        void addNonNei(int t, uint32_t u, uint32_t nonNei) {
// if(t == 1 && u == 0)
// printf("iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii\n");
            buffer[t][u * k + cntNonNei[t][u]++] = nonNei;
        }

        uint32_t getCntNonNei(int t, uint32_t u) { 
            return cntNonNei[t][u]; 
        }
        void pop(uint32_t t, uint32_t u) {
// if(t == 0 && u == 0)
// printf("pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp\n");
// if(cntNonNei[t][u]==0) {
//     printf("%u %upppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp\n", t, u);
// }
            --cntNonNei[t][u];
        }

        uint32_t * getBuffer(int t, uint32_t u) {
            return buffer[t].data() + u * k;
        }
    } nonNei;


    void plexEnum(uint32_t deep, parameter r, parameter c, parameter x, parameter sx);
    void impPlexBranch(uint32_t deep, parameter r, parameter c, parameter x, parameter sx, uint32_t t);

    parameter plexUpdateC(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
                        parameter c);
    parameter plexUpdateX(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
                        parameter x, parameter sx);

public:
    biplex(const std::string & fPath, int mode = 1, int k = 2):k(k) {
        if(mode == 0) g = new biGraph(fPath, 0);
        else g = new biGraph(fPath);
        printf("load graph n1 %u n2 %u m %u\n", g->n[0], g->n[1], g->m);fflush(stdout);  
        ws.resize(g->n[0] + g->n[1]);

        nonNei.init(g->n[0], g->n[1], k);
        parameter tmpC = {g->n[0], g->n[1]};

        for(int t = 0; t <= 1; t++) {
            C[t].resize(g->n[t]);
            X[t].resize(g->n[t]);
            R[t].resize(g->n[t]);
        }
    }
    ~biplex() { delete g; }

    void run();
};

#endif

//     struct nonNeiMatainer{
//         std::vector<uint32_t> buffer[2];
//         //cntNonNei[deep][i] is the count of non-neighbors of C[i]
//         std::vector<std::vector<uint32_t>> cntNonNei[2];
//         // std::vector<std::vector<uint32_t>> cand[2];
//         uint32_t n[2];
//         int k;

//         void init(uint32_t n1, uint32_t n2, int k) {
//             this->n[0] = n1; this->n[1] = n2; 
//             this->k = k;
//             buffer[0].resize(n1 * k);
//             buffer[1].resize(n2 * k);
//             cntNonNei[0].resize(n1 + 1);
//             cntNonNei[1].resize(n2 + 1);
//             cntNonNei[0][0].resize(n1);
//             cntNonNei[1][0].resize(n2);
//             // cand[0][0].resize(n1);
//             // cand[0][1].resize(n2);

//             // for(uint32_t i = 0; i < n1; i++) cand[0][0][i] = i;
//             // for(uint32_t i = 0; i < n2; i++) cand[0][1][i] = i;
//         }

//         void resizeNextDeep(uint32_t deep, parameter & c) {
// // printf("deep %u, c[0] %u, c[1] %u\n", deep, c[0], c[1]);fflush(stdout);
// // printf("cntNonNei[0][deep + 1].size() %u, cntNonNei[0][deep].size %u\n", 
// //     cntNonNei[0][deep + 1].size(), cntNonNei[0][deep].size());fflush(stdout);
// // printf("cntNonNei[1][deep + 1].size() %u, cntNonNei[1][deep].size %u\n", 
// //     cntNonNei[0][deep + 1].size(), cntNonNei[0][deep].size());fflush(stdout);
//             if(cntNonNei[0][deep + 1].size() < c[0]) cntNonNei[0][deep + 1].resize(c[0]);
//             if(cntNonNei[1][deep + 1].size() < c[1]) cntNonNei[1][deep + 1].resize(c[1]);
//         }

//         void addNonNei(int t, uint32_t u, uint32_t nonNei, uint32_t deep, LinearSet * C) {
//             buffer[t][u * k + cntNonNei[t][deep][C[t].pos(u)]++] = nonNei;
//         }
//         void addNonNeiIdx(int t, uint32_t u, uint32_t nonNei, uint32_t deep, uint32_t i) {
//             buffer[t][u * k + cntNonNei[t][deep][i]++] = nonNei;
//         }

//         uint32_t getCntNonNei(int t, uint32_t deep, uint32_t u, LinearSet * C) { 
//             return cntNonNei[t][deep][C[t].pos(u)]; 
//         }
//         uint32_t getCntNonNeiIdx(int t, uint32_t deep, uint32_t i) {
//             return cntNonNei[t][deep][i]; 
//         }
//         // void pop(uint32_t u, uint32_t deep) {
//         //     --cntNonNei[deep][u];
//         // }

//         // uint32_t * operator [] (int t, uint32_t u) {
//         //     return buffer[t].data() + u * k;
//         // }
//         uint32_t * getBuffer(int t, uint32_t u) {
//             return buffer[t].data() + u * k;
//         }
//     } nonNeiC, nonNeiR, nonNeiX;
    