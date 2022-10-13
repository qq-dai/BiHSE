#include "biplex.h"
// #define DEBUG

void biplex::run() {
#ifdef DEBUG
g->print();printf("\n");
#endif
    plexEnum(0, {0, 0}, {g->n[0], g->n[1]}, {0, 0}, {0, 0});

    printf("total cnt : %llu\n", totalMaximalBiPlex);
}

void biplex::plexEnum(uint32_t deep, parameter r, parameter c, parameter x, parameter sx) {
#ifdef DEBUG
printf("\n");
for(unsigned i = 0; i < deep; i++) printf(" ");
printf("deep %u\n", deep);fflush(stdout);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
#endif
    
    for(int t = 0; t <= 1; t++) {
        uint32_t z = t ^ 1;
        if(c[t] == 0) {
            if(x[z] == sx[z]) {
                parameter newC = c, newR = r, newX = x, newSX = sx;
                for(uint32_t i = 0; i < newC[z]; ) {
                    uint32_t u = C[z][i];

                    if(nonNei.getCntNonNei(z, u) == 0) {
                        R[z].changeTo(u, newR[z]++);
                        C[z].changeTo(u, --newC[z]);

                        for(uint32_t j = 0; j < newC[t]; j++) {
                            uint32_t v = C[t][j];
                            if(!g->connect(u, v, z)) {
                                nonNei.addNonNei(t, v, u);
                            }
                        }

                        uint32_t tmp = x[t];
                        for(uint32_t j = newSX[t]; j < tmp; ) {
                            uint32_t v = X[t][j];
                            if(g->connect(u, v, z)) {
                                X[t].changeTo(v, --tmp);
                            }
                            else j++;
                        }
                        newSX[t] = tmp;
                    }
                    else i++;
                }

                impPlexBranch(deep, newR, newC, newX, newSX, t);

                for(uint32_t i = 0; i < c[z]; i++) {
                    uint32_t u = C[z][i];

                    if(nonNei.getCntNonNei(z, u) == 0) {
                        for(uint32_t j = 0; j < newC[t]; j++) {
                            uint32_t v = C[t][j];
                            if(!g->connect(u, v, z)) {
                                nonNei.pop(t, v);
                            }
                        }
                    }
                }
            }
            return;
        }
    }
#ifdef DEBUG
// printf("deep2 %u\n", deep);fflush(stdout);
#endif
    //find the pivot with maximal degree in cand
    uint32_t maxI[2] = {0}, maxE[2] = {0}, maxUV[2] = {C[0][0], C[1][0]};

    for(uint32_t t = 0; t < 2; t++) {
        uint32_t z = t ^ 1;//t + z = 1

        for(uint32_t i = 0; i < c[t]; i++) {//scan all vertices in C
            uint32_t u = C[t][i];
            uint32_t e = 0;

            if(g->deg(u, t) > c[z]) {
                for(uint32_t j = 0; j < c[z]; j++) {
                    uint32_t v = C[z][j];
                    if(g->connect(u, v, t)) {
                        e++;
                    }
                }
            }
            else {
                for(uint32_t j = g->p[t][u]; j < g->p[t][u + 1]; j++) {
                    uint32_t v = g->e[t][j];
                    if(C[z].pos(v) < c[z]) {//in C
                        e++;
                    }
                }
            }

            if(e > maxE[z]) {
                maxE[z] = e;
                maxI[t] = i;
                maxUV[t] = u;
            }
        }
    }
#ifdef DEBUG
// printf("deep3 %u\n", deep);fflush(stdout);
#endif
    //t is the side of pivot, and S[z].cSize() - maxE[z] is smaller than S[t].cSize() - maxE[t]
    uint32_t t = 1;
    if(c[0] - maxE[0] >= c[1] - maxE[1]) {
        t = 0;
    }
    uint32_t z = t ^ 1;

    if(maxE[0] == 0) {
#ifdef DEBUG
printf("no edge\n");
#endif
        plexEnum(deep + 1, r, {0, c[1]}, x, sx);
#ifdef DEBUG
printf("    return1 deep %u\n", deep);
#endif


        plexEnum(deep + 1, r, {c[0], 0}, x, sx);
#ifdef DEBUG
printf("    return2 deep %u\n", deep);
#endif
        return;
    }

    //ws[deep] is the nonNeighbors of pivot in cand
#ifdef DEBUG
printf("pivot %u, t %u, degree %u\n", maxUV[t], t, maxE[z]);fflush(stdout);
#endif
    //run pivot
    R[t].changeTo(maxUV[t], r[t]++);
    C[t].changeTo(maxUV[t], --c[t]);

    //make sure newC[0] points to U
    parameter newC = plexUpdateC(t, deep, r, maxUV[t], c);
    parameter newSX = plexUpdateX(t, deep, r, maxUV[t], x, sx);
#ifdef DEBUG
// printf("deep4.1 %u\n", deep);fflush(stdout);
#endif
    uint32_t wsSize = c[z] - newC[z];
#ifdef DEBUG
// printf("deep4.1 %u, wsSize %u\n", deep, wsSize);fflush(stdout);
#endif
    if(ws[deep].size() < wsSize) {
#ifdef DEBUG
// printf("deep4.1 %u, wsSize %u\n", deep, wsSize);fflush(stdout);
#endif
        ws[deep].resize(wsSize * 2);
    }
#ifdef DEBUG
// printf("deep4.2 %u\n", deep);fflush(stdout);
#endif
    memcpy(ws[deep].data(), C[z].begin() + newC[z], sizeof(uint32_t) * wsSize);
#ifdef DEBUG
// printf("deep5 %u\n", deep);fflush(stdout);
#endif

    //update non nei in C X R
    for(uint32_t i = 0; i < wsSize; i++) {
        uint32_t w = ws[deep][i];

        nonNei.addNonNei(z, w, maxUV[t]);
// printf("PivotaddNonNei %u %u, after degree %u iiiiiiiiiiiiiiiiiiiiiiiiiii\n", 
    // z, w, nonNei.getCntNonNei(z, w));
    }
    for(uint32_t i = sx[z]; i < x[z]; i++) {
        uint32_t v = X[z][i];
        if(!g->connect(maxUV[t], v, t))
            nonNei.addNonNei(z, v, maxUV[t]);
    }
    for(uint32_t i = 0; i < r[z]; i++) {
        uint32_t v = R[z][i];
        if(!g->connect(maxUV[t], v, t))
            nonNei.addNonNei(z, v, maxUV[t]);
    }
    
    plexEnum(deep + 1, r, newC, x, newSX);
#ifdef DEBUG
printf("    return3 deep %u\n", deep);

#endif

    //rm pivot
    X[t].changeTo(maxUV[t], x[t]++);
    // R[t].changeTo(maxUV[t], --r[t]);
    --r[t];
    
    for(uint32_t i = 0; i < wsSize; i++) {
        uint32_t w = ws[deep][i];
        nonNei.pop(z, w);
    }
    for(uint32_t i = sx[z]; i < x[z]; i++) {
        uint32_t v = X[z][i];
        if(!g->connect(maxUV[t], v, t))
            nonNei.pop(z, v);
    }
    for(uint32_t i = 0; i < r[z]; i++) {
        uint32_t v = R[z][i];
        if(!g->connect(maxUV[t], v, t)) {
            nonNei.pop(z, v);
        }
    }

#ifdef DEBUG

printf("pivot %u, t %u, degree %u\n", maxUV[t], t, maxE[z]);fflush(stdout);

printf("    before h node, deep %u  i %u\n", deep, 0);
printf("deep %u\n", deep);fflush(stdout);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("Ws:");
for(uint32_t i = 0; i < wsSize; i++) printf("%u ", ws[deep][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
#endif

    for(uint32_t i = 0; i < wsSize; i++) {
        uint32_t w = ws[deep][i];

        R[z].changeTo(w, r[z]++);
        C[z].changeTo(w, --c[z]);

        parameter newC = plexUpdateC(z, deep, r, w, c);
        parameter newSX = plexUpdateX(z, deep, r, w, x, sx);
// printf("     w %u %u\n", w, t);
// printf("oldC %u %u\n", c[0], c[1]);
// printf("newC %u %u\n", newC[0], newC[1]);
// for(uint32_t i = newSX[t]; i < x[t]; i++) printf("w%u-%u:%d ", w, X[t][i], g->connect(w, X[t][i], z));printf("\n");
#ifdef DEBUG
printf("    in h node, deep %u  i %u, t%u\n", deep, i, t);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("aC_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("aC_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < newC[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < newC[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = newSX[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = newSX[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
fflush(stdout);
#endif
        //scan c x r
        for(uint32_t i = 0; i < newC[t]; i++)
            if(!g->connect(w, C[t][i], z)) {
// printf("addCNei:t %u c %u nonNei %u, %u\n", t, C[t][i], w, nonNei.getCntNonNei(t, C[t][i]));fflush(stdout);
                nonNei.addNonNei(t, C[t][i], w);
            }
        for(uint32_t i = newSX[t]; i < x[t]; i++) {
            if(!g->connect(w, X[t][i], z)) {
// printf("addXNei:%u %u %u, %u\n", t, X[t][i], w, nonNei.getCntNonNei(t, X[t][i]));fflush(stdout);
                nonNei.addNonNei(t, X[t][i], w);
            }
        }
            
        for(uint32_t i = 0; i < r[t]; i++)
            if(!g->connect(w, R[t][i], z))
                nonNei.addNonNei(t, R[t][i], w);

        plexEnum(deep + 1, r, newC, x, newSX);

#ifdef DEBUG
printf("    return4 deep %u\n", deep);
#endif

        for(uint32_t i = 0; i < newC[t]; i++) {
            if(!g->connect(w, C[t][i], z)) {
// printf("NonNeiPop1 %u %u, deg %u\n", t, C[t][i], nonNei.getCntNonNei(t, C[t][i]));
                nonNei.pop(t, C[t][i]);
            }
        }
        for(uint32_t i = newSX[t]; i < x[t]; i++) {
            if(!g->connect(w, X[t][i], z)) {
// printf("NonNeiPop2 t %u X %u, deg %u\n", t, X[t][i], nonNei.getCntNonNei(t, X[t][i]));
                nonNei.pop(t, X[t][i]);
            }
        }
        for(uint32_t i = 0; i < r[t]; i++) {
            if(!g->connect(w, R[t][i], z)) {
// printf("NonNeiPop3 %u %u, deg %u\n", t, R[t][i], nonNei.getCntNonNei(t, R[t][i]));
                nonNei.pop(t, R[t][i]);
            }
        }

        X[z].changeTo(w, x[z]++);
        // R[z].changeTo(w, --r[z]);
        --r[z];

#ifdef DEBUG
printf("    out h node, deep %u  i %u, t%u, w%u\n", deep, i, t, w);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("aC_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("aC_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < newC[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < newC[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = newSX[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = newSX[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
fflush(stdout);
#endif
    }
    
    for(uint32_t i = 0; i < wsSize; i++) X[z].changeTo(ws[deep][i], --x[z]);
    X[t].changeTo(maxUV[t], --x[t]);
}

void biplex::impPlexBranch(uint32_t deep, parameter r, parameter c, parameter x, parameter sx, uint32_t t) {
#ifdef DEBUG
for(unsigned i = 0; i < deep; i++) printf(" ");
printf("imp_deep %u\n", deep);fflush(stdout);
printf("R_l:");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
printf("R_r:");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n");
printf("C_l:");
for(uint32_t i = 0; i < c[0]; i++) printf("%u ", C[0][i]);printf("\n");
printf("C_r:");
for(uint32_t i = 0; i < c[1]; i++) printf("%u ", C[1][i]);printf("\n");
printf("X_l:");
for(uint32_t i = sx[0]; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("X_r:");
for(uint32_t i = sx[1]; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
printf("aX_l:");
for(uint32_t i = 0; i < x[0]; i++) printf("%u ", X[0][i]);printf("\n");
printf("aX_r:");
for(uint32_t i = 0; i < x[1]; i++) printf("%u ", X[1][i]);printf("\n");
fflush(stdout);
#endif
    
    uint32_t z = t ^ 1;
    //move C[z] into R[z], let X[t]=0

    if(c[t] == 0) {
        if(x[t] - sx[t] + x[z] - sx[z] == 0 && r[t] > 0 && r[z] > 0) {
            totalMaximalBiPlex++;
#ifdef DEBUG
printf("    get a answer\n");
for(uint32_t i = 0; i < r[0]; i++) printf("%u ", R[0][i]);printf("\n");
for(uint32_t i = 0; i < r[1]; i++) printf("%u ", R[1][i]);printf("\n\n");
#endif
        }

        return ;
    }
    
    uint32_t pivot = C[z][0], pivotIdx = 0;
    uint32_t nd = nonNei.getCntNonNei(z, pivot);
    for(uint32_t i = 1; i < c[z]; i++) {
        uint32_t v = C[z][i];
        if(nd > nonNei.getCntNonNei(z, v)) {
            nd = nonNei.getCntNonNei(z, v);
            pivot = v;
            pivotIdx = i;
        }
    }
    
    //build Cv/Pv
    uint32_t wsSize = 0;

    uint32_t costkd = 0, cntNonNeiPivot = nd;
    for(uint32_t i = 0; i < cntNonNeiPivot; i++) {
        uint32_t v = nonNei.getBuffer(z, pivot)[i];
        costkd += g->deg(v, t);
    }
    if(costkd < c[z] * cntNonNeiPivot) {
        wsSize = 0;
        for(uint32_t i = 0; i < cntNonNeiPivot; i++) {
            uint32_t v = nonNei.getBuffer(z, pivot)[i];
            uint32_t tmp = c[z];
            
            for(uint32_t j = g->p[t][v]; j < g->p[t][v + 1]; j++) {
                uint32_t w = g->e[t][j];

                if(wsSize <= C[z].pos(w) && C[z].pos(w) < c[z]) {
                    C[z].changeTo(w, --tmp);
                }
            }

            wsSize = tmp;
        }
    }
    else {
        wsSize = c[z];
        for(uint32_t j = 0; j < c[z]; j++) {
            uint32_t u = C[z][j];
            bool connectAll = true;

            for(uint32_t i = 0; i < cntNonNeiPivot; i++) {
                uint32_t v = nonNei.getBuffer(z, pivot)[i];
                if(!g->connect(u, v, z)) {
                    connectAll = false;
                    break;
                }
            }

            if(connectAll) {
                C[z].changeToByPos(j, --wsSize);
            }
        }
    }

    if(ws[deep].size() < wsSize) {
        ws[deep].resize(wsSize * 2);
    }
    memcpy(ws[deep].data(), C[z].begin(), sizeof(uint32_t) * wsSize);

    for(uint32_t i = 0; i < wsSize; i++) {
        uint32_t v = ws[deep][i];

        R[z].changeTo(v, r[z]++);
        C[z].changeTo(v, --c[z]);

        parameter newC = plexUpdateC(z, deep, r, v, c);
        parameter newSX = plexUpdateX(z, deep, r, v, x, sx);

        // for(uint32_t i = 0; i < newC[t]; i++)
        //     if(!g->connect(v, C[t][i], z)) nonNei.addNonNei(t, C[t][i], v);
        // for(uint32_t i = 0; i < newX[t]; i++) 
        //     if(!g->connect(v, X[t][i], z)) nonNei.addNonNei(t, X[t][i], v);
        for(uint32_t i = 0; i < r[t]; i++) 
            if(!g->connect(v, R[t][i], z)) 
                nonNei.addNonNei(t, R[t][i], v);

        // plexEnum(deep + 1, r, newC, newX);
        impPlexBranch(deep + 1, r, newC, x, newSX, t);

        // for(uint32_t i = 0; i < newC[t]; i++) if(!g->connect(w, C[t][i], z)) nonNei.pop(t, C[t][i]);
        // for(uint32_t i = 0; i < newX[t]; i++) if(!g->connect(w, X[t][i], z)) nonNei.pop(t, X[t][i]);
        for(uint32_t i = 0; i < r[t]; i++) 
            if(!g->connect(v, R[t][i], z)) 
                nonNei.pop(t, R[t][i]);

        X[z].changeTo(v, x[z]++);
        // R[z].changeTo(v, --r[z]);
        --r[z];
    }

    for(uint32_t i = 0; i < wsSize; i++) X[z].changeTo(ws[deep][i], --x[z]);
}

biplex::parameter biplex::plexUpdateC(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
    parameter c) 
{    
    uint32_t z = t ^ 1;
    parameter newC;
    newC[t] = c[t];
    newC[z] = 0;

    for(uint32_t i = 0; i < c[z]; i++) {
        uint32_t v = C[z][i];
        if(g->connect(u, v, t)) {
            C[z].changeTo(v, newC[z]++);
// printf("inUpdateC:%u C %u, newC %u\n", u, v, newC[z]);
        }
        else {
            if(nonNei.getCntNonNei(z, v) < k - 1) {
                // nonNei.addNonNeiIdx(z, v, u, deep + 1, i);
                C[z].changeTo(v, newC[z]++);
// printf("inUpdateCnonNei:%u C %u, newC %u\n", u, v, newC[z]);
// printf("nonNei.getCntNonNei(z, v) %u\n", nonNei.getCntNonNei(z, v));
            }
        }
    }

    if(newC[t]> 0)
    for(uint32_t i = 0; i < nonNei.getCntNonNei(t, u); i++) {
        uint32_t v = nonNei.getBuffer(t, u)[i];
        
        if(nonNei.getCntNonNei(z, v) == k) {
            if(g->deg(v, z) < c[t]) {
                uint32_t tmp = 0;
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    if(C[t].pos(w) < newC[t]) {
                        C[t].changeTo(w, tmp++);
                    }
                }
                newC[t] = tmp;
            }
            else {
                for(uint32_t j = 0; j < newC[t]; ) {
                    uint32_t w = C[t][j];
                    
                    if(!g->connect(v, w, z)) {
                        C[t].changeTo(w, --newC[t]);
                    }
                    else j++;
                }
            }
        }
    }

    return newC;
}

biplex::parameter biplex::plexUpdateX(uint32_t t, uint32_t deep, parameter r, uint32_t u, 
                        parameter x, parameter sx)
{
    uint32_t z = t ^ 1;
    parameter newSX;
    newSX[t] = sx[t];
    newSX[z] = x[z];

    for(uint32_t i = sx[z]; i < newSX[z]; ) {
        uint32_t v = X[z][i];
        if(g->connect(u, v, t)) {
            X[z].changeTo(v, --newSX[z]);
        }
        else {
            if(nonNei.getCntNonNei(z, v) + 1 < k) {
                // nonNei.addNonNeiIdx(z, v, u, deep + 1, i);
                X[z].changeTo(v, --newSX[z]);
            }
            else i++;
        }
    }

    if(sx[t] < x[t])
    for(uint32_t i = 0; i < nonNei.getCntNonNei(t, u); i++) {
        uint32_t v = nonNei.getBuffer(t, u)[i];
        
        if(nonNei.getCntNonNei(z, v) == k) {
            if(g->deg(v, z) < x[t] - newSX[t]) {
                uint32_t tmp = x[t];
                newSX[t] = sx[t];
                for(uint32_t j = g->p[z][v]; j < g->p[z][v + 1]; j++) {
                    uint32_t w = g->e[z][j];
                    uint32_t pw = X[t].pos(w);
                    if(newSX[t] <= pw && pw < tmp) {
                        X[t].changeTo(w, --tmp);
                    }
                }
                newSX[t] = tmp;
            }
            else {
                newSX[t] = x[t];
                for(uint32_t j = sx[t]; j < newSX[t]; ) {
                    uint32_t w = X[t][j];
                    
                    if(g->connect(v, w, z)) {
                        X[t].changeTo(w, --newSX[t]);
                    }
                    else j++;
                }
            }
        }
    }

    return newSX;
}
