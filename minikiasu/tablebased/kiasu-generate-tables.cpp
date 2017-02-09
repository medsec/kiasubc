// ===================================================================
// @last-modified 2016-06-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include "kiasu.h"

uint8_t SBOX[KIASU_SBOXLEN] = { 
    14,  4, 11,  2,  3,  8,  0,  9,  1, 10,  7, 15,  6, 12,  5, 13 
};
uint8_t INVERSE_SBOX[KIASU_SBOXLEN] = { 
     6,  8,  3,  4,  1, 14, 12, 10,  5,  7,  9,  2, 13, 15,  0, 11
};

static const uint8_t times2[KIASU_SBOXLEN] = { 
    0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x3, 0x1, 0x7, 0x5, 0xb, 0x9, 0xf, 0xd
};
static const uint8_t times3[KIASU_SBOXLEN] = { 
    0x0, 0x3, 0x6, 0x5, 0xc, 0xf, 0xa, 0x9, 0xb, 0x8, 0xd, 0xe, 0x7, 0x4, 0x1, 0x2
};
static const uint8_t times9[KIASU_SBOXLEN] = {
    0x0, 0x9, 0x1, 0x8, 0x2, 0xb, 0x3, 0xa, 0x4, 0xd, 0x5, 0xc, 0x6, 0xf, 0x7, 0xe
};
static const uint8_t timesB[KIASU_SBOXLEN] = {
    0x0, 0xb, 0x5, 0xe, 0xa, 0x1, 0xf, 0x4, 0x7, 0xc, 0x2, 0x9, 0xd, 0x6, 0x8, 0x3
};
static const uint8_t timesD[KIASU_SBOXLEN] = {
    0x0, 0xd, 0x9, 0x4, 0x1, 0xc, 0x8, 0x5, 0x2, 0xf, 0xb, 0x6, 0x3, 0xe, 0xa, 0x7
};
static const uint8_t timesE[KIASU_SBOXLEN] = {
    0x0, 0xe, 0xf, 0x1, 0xd, 0x3, 0x2, 0xc, 0x9, 0x7, 0x6, 0x8, 0x4, 0xa, 0xb, 0x5
};

// ---------------------------------------------------------------------

static uint8_t gf_16_mult(uint8_t a, uint8_t b) {
    uint8_t p = 0, i = 0, hbs = 0;

    for (i = 0; i < 4; i++) {
        if (b & 1) {
            p ^= a;
        }
        hbs = a & 0x08;
        a <<= 1;

        if (hbs) {
            a ^= 0x13; // = x^4 + x + 1 = 0001 0011
        }
        
        b >>= 1;
    }
    return (uint8_t)p;
}

#define times2_macro(x) gf_16_mult(x, 0x2)
#define times3_macro(x) gf_16_mult(x, 0x3)
#define times9_macro(x) gf_16_mult(x, 0x9)
#define timesB_macro(x) gf_16_mult(x, 0xB)
#define timesD_macro(x) gf_16_mult(x, 0xD)
#define timesE_macro(x) gf_16_mult(x, 0xE)

// ---------------------------------------------------------------------

static void compute_byte_inverse_sbox() {
    puts("static const uint8_t BYTE_INVERSE_SBOX[256] = {");

    for (int i = 0; i < 16; ++i) {
        uint8_t y1 = INVERSE_SBOX[i];

        for (int j = 0; j < 16; ++j) {
            uint8_t y2 = INVERSE_SBOX[j];
            uint8_t y = (y1 << 4) | y2;

            printf("0x%02x, ", y);
        }
        puts("");
    }

    puts("};");
}

// ---------------------------------------------------------------------

static void compute_mult_tables() {
    uint8_t i, j;

    puts("static const uint8_t times2[KIASU_SBOXLEN] = {");
    for (i = 0; i < KIASU_SBOXLEN; ++i) {
        j = times2_macro(i);
        printf("0x%x, ", j);
    }
    puts("\n};");

    puts("static const uint8_t times3[KIASU_SBOXLEN] = {");
    for (i = 0; i < KIASU_SBOXLEN; ++i) {
        j = times3_macro(i);
        printf("0x%x, ", j);
    }
    puts("\n};");

    puts("static const uint8_t times9[KIASU_SBOXLEN] = {");
    for (i = 0; i < KIASU_SBOXLEN; ++i) {
        j = times9_macro(i);
        printf("0x%x, ", j);
    }
    puts("\n};");
    
    puts("static const uint8_t timesB[KIASU_SBOXLEN] = {");
    for (i = 0; i < KIASU_SBOXLEN; ++i) {
        j = timesB_macro(i);
        printf("0x%x, ", j);
    }
    puts("\n};");

    puts("static const uint8_t timesD[KIASU_SBOXLEN] = {");
    for (i = 0; i < KIASU_SBOXLEN; ++i) {
        j = timesD_macro(i);
        printf("0x%x, ", j);
    }
    puts("\n};");
    
    puts("static const uint8_t timesE[KIASU_SBOXLEN] = {");
    for (i = 0; i < KIASU_SBOXLEN; ++i) {
        j = timesE_macro(i);
        printf("0x%x, ", j);
    }
    puts("\n};");
}

// ---------------------------------------------------------------------

static void compute_tables() {
    uint32_t z_1, z_2, z_3, z_4;

    puts("uint64_t Te1[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_1 = 0L;
        z_1 = ((uint32_t)times2[y_1] << 12) 
            | ((uint32_t)y_1 <<  8)
            | ((uint32_t)y_1 <<  4)
            | ((uint32_t)times3[y_1] << 0); // Row 1
        z_1 <<= 16;

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_2 = ((uint32_t)times3[y_2] << 12) 
                | ((uint32_t)times2[y_2] << 8)
                | ((uint32_t)y_2 << 4)
                |  (uint32_t)y_2; // Row 2

            printf("0x%08x%08xL, ", z_1, z_2);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te2[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_3 = ((uint32_t)y_1 << 12) 
            | ((uint32_t)times3[y_1] <<  8)
            | ((uint32_t)times2[y_1] <<  4)
            | ((uint32_t)y_1); // Row 3
        z_3 <<= 16;

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_4 = ((uint32_t)y_2 << 12) 
                | ((uint32_t)y_2 << 8)
                | ((uint32_t)times3[y_2] << 4)
                |  (uint32_t)times2[y_2]; // Row 4

            printf("0x%08x%08xL, ", z_4, z_3);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te3[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_1 = ((uint32_t)times2[y_1] << 12) 
            | ((uint32_t)y_1 <<  8)
            | ((uint32_t)y_1 <<  4)
            |  (uint32_t)times3[y_1]; // Row 1

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_2 = ((uint32_t)times3[y_2] << 12) 
                | ((uint32_t)times2[y_2] <<  8)
                | ((uint32_t)y_2 <<  4)
                |  (uint32_t)y_2; // Row 2
            z_3 = (z_2 << 16) | z_1;

            printf("0x%08x%08xL, ", z_3, 0);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te4[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_3 = ((uint32_t)y_1 << 12) 
            | ((uint32_t)times3[y_1] <<  8)
            | ((uint32_t)times2[y_1] <<  4)
            | ((uint32_t)y_1); // Row 3

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_4 = ((uint32_t)y_2 << 12) 
                | ((uint32_t)y_2 << 8)
                | ((uint32_t)times3[y_2] << 4)
                |  (uint32_t)times2[y_2]; // Row 4
            z_1 = (z_4 << 16) | z_3;

            printf("0x%08x%08xL, ", 0, z_1);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te5[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_1 = ((uint32_t)times2[y_1] << 12) 
            | ((uint32_t)y_1 <<  8)
            | ((uint32_t)y_1 <<  4)
            |  (uint32_t)times3[y_1]; // Row 1
        z_1 <<= 16;

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_2 = ((uint32_t)times3[y_2] << 12) 
                | ((uint32_t)times2[y_2] <<  8)
                | ((uint32_t)y_2 <<  4)
                |  (uint32_t)y_2; // Row 2

            printf("0x%08x%08xL, ", z_2, z_1);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te6[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_3 = ((uint32_t)y_1 << 12) 
            | ((uint32_t)times3[y_1] <<  8)
            | ((uint32_t)times2[y_1] <<  4)
            | ((uint32_t)y_1); // Row 3
        z_3 <<= 16;

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_4 = ((uint32_t)y_2 << 12) 
                | ((uint32_t)y_2 << 8)
                | ((uint32_t)times3[y_2] << 4)
                |  (uint32_t)times2[y_2]; // Row 4

            printf("0x%08x%08xL, ", z_3, z_4);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te7[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_1 = ((uint32_t)times2[y_1] << 12) 
            | ((uint32_t)y_1 <<  8)
            | ((uint32_t)y_1 <<  4)
            |  (uint32_t)times3[y_1]; // Row 1

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_2 = ((uint32_t)times3[y_2] << 12) 
                | ((uint32_t)times2[y_2] <<  8)
                | ((uint32_t)y_2 <<  4)
                |  (uint32_t)y_2; // Row 2
            z_3 = (z_2 << 16) | z_1;

            printf("0x%08x%08xL, ", 0, z_3);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Te8[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        uint8_t y_1 = SBOX[x_1];

        z_3 = ((uint32_t)y_1 << 12) 
            | ((uint32_t)times3[y_1] <<  8)
            | ((uint32_t)times2[y_1] <<  4)
            | ((uint32_t)y_1); // Row 3

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            uint8_t y_2 = SBOX[x_2];

            z_4 = ((uint32_t)y_2 << 12) 
                | ((uint32_t)y_2 << 8)
                | ((uint32_t)times3[y_2] << 4)
                |  (uint32_t)times2[y_2]; // Row 4
            z_1 = (z_4 << 16) | z_3;

            printf("0x%08x%08xL, ", z_1, 0);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td1[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_1 = ((uint32_t)timesE[x_1] << 12) 
            | ((uint32_t)times9[x_1] <<  8)
            | ((uint32_t)timesD[x_1] <<  4)
            |  (uint32_t)timesB[x_1]; // Row 1

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_2 = ((uint32_t)timesB[x_2] << 12) 
                | ((uint32_t)timesE[x_2] << 8)
                | ((uint32_t)times9[x_2] << 4)
                |  (uint32_t)timesD[x_2]; // Row 2
            // z_3 = (z_1 << 16) | z_2;

            // printf("0x%08x%08xL, ", z_3, 0);
            z_3 = (z_1 ^ z_2) << 16;
            printf("0x%08x%08xL, ", z_3, 0);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td2[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_3 = ((uint32_t)timesD[x_1] << 12) 
            | ((uint32_t)timesB[x_1] <<  8)
            | ((uint32_t)timesE[x_1] <<  4)
            |  (uint32_t)times9[x_1]; // Row 3

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_4 = ((uint32_t)times9[x_2] << 12) 
                | ((uint32_t)timesD[x_2] << 8)
                | ((uint32_t)timesB[x_2] << 4)
                |  (uint32_t)timesE[x_2]; // Row 4
            // z_1 = (z_3 << 16) | z_4;

            // printf("0x%08x%08xL, ", 0, z_1);
            z_1 = (z_3 ^ z_4) << 16;
            printf("0x%08x%08xL, ", z_1, 0);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td3[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_1 = ((uint32_t)timesE[x_1] << 12) 
            | ((uint32_t)times9[x_1] <<  8)
            | ((uint32_t)timesD[x_1] <<  4)
            |  (uint32_t)timesB[x_1]; // Row 1

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_2 = ((uint32_t)timesB[x_2] << 12) 
                | ((uint32_t)timesE[x_2] << 8)
                | ((uint32_t)times9[x_2] << 4)
                |  (uint32_t)timesD[x_2]; // Row 2
            // z_2 <<= 16;

            // printf("0x%08x%08xL, ", z_1, z_2);
            z_3 = z_1 ^ z_2;
            printf("0x%08x%08xL, ", z_3, 0);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td4[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_3 = ((uint32_t)timesD[x_1] << 12) 
            | ((uint32_t)timesB[x_1] <<  8)
            | ((uint32_t)timesE[x_1] <<  4)
            |  (uint32_t)times9[x_1]; // Row 3

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_4 = ((uint32_t)times9[x_2] << 12) 
                | ((uint32_t)timesD[x_2] << 8)
                | ((uint32_t)timesB[x_2] << 4)
                |  (uint32_t)timesE[x_2]; // Row 4
            // z_4 <<= 16;

            // printf("0x%08x%08xL, ", z_4, z_3);
            z_1 = z_3 ^ z_4;
            printf("0x%08x%08xL, ", z_1, 0);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td5[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_1 = ((uint32_t)timesE[x_1] << 12) 
            | ((uint32_t)times9[x_1] <<  8)
            | ((uint32_t)timesD[x_1] <<  4)
            |  (uint32_t)timesB[x_1]; // Row 1

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_2 = ((uint32_t)timesB[x_2] << 12) 
                | ((uint32_t)timesE[x_2] << 8)
                | ((uint32_t)times9[x_2] << 4)
                |  (uint32_t)timesD[x_2]; // Row 2
            // z_3 = (z_1 << 16) | z_2;

            // printf("0x%08x%08xL, ", 0, z_3);
            z_3 = (z_1 ^ z_2) << 16;
            printf("0x%08x%08xL, ", 0, z_3);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td6[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_3 = ((uint32_t)timesD[x_1] << 12) 
            | ((uint32_t)timesB[x_1] <<  8)
            | ((uint32_t)timesE[x_1] <<  4)
            |  (uint32_t)times9[x_1]; // Row 3

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_4 = ((uint32_t)times9[x_2] << 12) 
                | ((uint32_t)timesD[x_2] << 8)
                | ((uint32_t)timesB[x_2] << 4)
                |  (uint32_t)timesE[x_2]; // Row 4
            // z_1 = (z_3 << 16) | z_4;

            // printf("0x%08x%08xL, ", z_1, 0);
            z_1 = (z_3 ^ z_4) << 16;
            printf("0x%08x%08xL, ", 0, z_1);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td7[256] = {");

    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_1 = ((uint32_t)timesE[x_1] << 12) 
            | ((uint32_t)times9[x_1] <<  8)
            | ((uint32_t)timesD[x_1] <<  4)
            |  (uint32_t)timesB[x_1]; // Row 1

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_2 = ((uint32_t)timesB[x_2] << 12) 
                | ((uint32_t)timesE[x_2] << 8)
                | ((uint32_t)times9[x_2] << 4)
                |  (uint32_t)timesD[x_2]; // Row 2
            // z_2 <<= 16;

            // printf("0x%08x%08xL, ", z_2, z_1);
            z_3 = z_1 ^ z_2;
            printf("0x%08x%08xL, ", 0, z_3);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");

    puts("uint64_t Td8[256] = {");

    // Bottom table
    for (uint8_t x_1 = 0; x_1 < 16; ++x_1) {
        z_3 = ((uint32_t)timesD[x_1] << 12) 
            | ((uint32_t)timesB[x_1] <<  8)
            | ((uint32_t)timesE[x_1] <<  4)
            |  (uint32_t)times9[x_1]; // Row 3

        for (uint8_t x_2 = 0; x_2 < 16; ++x_2) {
            z_4 = ((uint32_t)times9[x_2] << 12) 
                | ((uint32_t)timesD[x_2] << 8)
                | ((uint32_t)timesB[x_2] << 4)
                |  (uint32_t)timesE[x_2]; // Row 4
            // z_4 <<= 16;

            // printf("0x%08x%08xL, ", z_3, z_4);
            z_1 = z_3 ^ z_4;
            printf("0x%08x%08xL, ", 0, z_1);

            if (x_2 % 4 == 3) {
                puts("");
            }
        }
    }
    puts("};");
    
}

// ---------------------------------------------------------------------

int main() {
    compute_byte_inverse_sbox();
    compute_mult_tables();
    compute_tables();
}

