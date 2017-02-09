// ===================================================================
// @last-modified 2016-09-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstring>

#include "kiasu.h"
#include "kiasu-tables.h"

static const uint8_t SBOX[KIASU_SBOXLEN] = { 
    14,  4, 11,  2,  3,  8,  0,  9,  1, 10,  7, 15,  6, 12,  5, 13 
};
// static uint8_t INVERSE_SBOX[KIASU_SBOXLEN] = { 
//      6,  8,  3,  4,  1, 14, 12, 10,  5,  7,  9,  2, 13, 15,  0, 11
// };
static const uint8_t TIMES2[KIASU_SBOXLEN] = { 
    0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x3, 0x1, 0x7, 0x5, 0xb, 0x9, 0xf, 0xd
};
static const uint8_t TIMES3[KIASU_SBOXLEN] = { 
    0x0, 0x3, 0x6, 0x5, 0xc, 0xf, 0xa, 0x9, 0xb, 0x8, 0xd, 0xe, 0x7, 0x4, 0x1, 0x2
};

// ---------------------------------------------------------------------

static void precompute_trail_table(
    uint8_t trails[KIASU_SBOXLEN][KIASU_SBOXLEN][KIASU_SBOXLEN]) {

    for (size_t i = 0; i < KIASU_SBOXLEN; ++i) {
        for (size_t j = 0; j < KIASU_SBOXLEN; ++j) {
            for (size_t k = 0; k < KIASU_SBOXLEN; ++k) {
                trails[i][j][k] = 0;
            }
        }
    }

    for (uint8_t x = 0; x < KIASU_SBOXLEN; ++x) {
        uint8_t y = SBOX[x];

        for (uint8_t x_ = 0; x_ < KIASU_SBOXLEN; ++x_) {
            uint8_t y_ = SBOX[x_];
            uint8_t delta_x = (x ^ x_) & 0xF;
            uint8_t delta_y = (y ^ y_) & 0xF;
        
            trails[delta_x][delta_y][x] = 1;
        }
    }
}

static uint8_t sbox_trails[KIASU_SBOXLEN][KIASU_SBOXLEN][KIASU_SBOXLEN];

// ---------------------------------------------------------------------

#ifdef DEBUG
static void print_hex(const char label[], 
                      const uint8_t* x, 
                      const size_t n) {
    printf("%s ", label);

    for (size_t i = 0; i < n; ++i) {
        printf("%02x", x[i]);
    }

    puts("");
}

// ---------------------------------------------------------------------

static void print_block(const char label[], const uint8_t* x) {
    print_hex(label, x, KIASU_BLOCKLEN);
}
#endif

// ---------------------------------------------------------------------

void invert_mixcolumns(block out, const block in) {
    const uint64_t state = Td1[in[0]] ^ Td2[in[1]] ^ Td3[in[2]] ^ Td4[in[3]]
        ^ Td5[in[4]] ^ Td6[in[5]] ^ Td7[in[6]] ^ Td8[in[7]];
    out[0] = (state >> 56) & 0xFF;
    out[1] = (state >> 48) & 0xFF;
    out[2] = (state >> 40) & 0xFF;
    out[3] = (state >> 32) & 0xFF;
    out[4] = (state >> 24) & 0xFF;
    out[5] = (state >> 16) & 0xFF;
    out[6] = (state >>  8) & 0xFF;
    out[7] =  state        & 0xFF;
}

// ---------------------------------------------------------------------

void invert_shiftrows(block out, const block in) {
    out[0] = (in[0] & 0xF0) | (in[6] & 0x0F);
    out[1] = (in[5] & 0xF0) | (in[3] & 0x0F);
    out[2] = (in[2] & 0xF0) | (in[0] & 0x0F);
    out[3] = (in[7] & 0xF0) | (in[5] & 0x0F);

    out[4] = (in[4] & 0xF0) | (in[2] & 0x0F);
    out[5] = (in[1] & 0xF0) | (in[7] & 0x0F);
    out[6] = (in[6] & 0xF0) | (in[4] & 0x0F);
    out[7] = (in[3] & 0xF0) | (in[1] & 0x0F);
}

// ---------------------------------------------------------------------

static inline void invert_subbytes(block inout) {
    for (uint32_t i = 0; i < KIASU_BLOCKLEN; ++i) {
        inout[i] = BYTE_INVERSE_SBOX[inout[i]];
    }
}

// ---------------------------------------------------------------------

static inline void mixcolumns(block out, const block in) {
    uint8_t x1, x2, x3, x4;
    uint8_t y1, y2, y3, y4;

    for (uint32_t i = 0; i < KIASU_BLOCKLEN; i += 2) {
        x1 = (in[i  ] >> 4) & 0xF;
        x2 =  in[i  ]       & 0xF;
        x3 = (in[i+1] >> 4) & 0xF;
        x4 =  in[i+1]       & 0xF;

        y1 = TIMES2[x1] ^ TIMES3[x2] ^        x3  ^        x4;
        y2 =        x1  ^ TIMES2[x2] ^ TIMES3[x3] ^        x4;
        y3 =        x1  ^        x2  ^ TIMES2[x3] ^ TIMES3[x4];
        y4 = TIMES3[x1] ^        x2  ^        x3  ^ TIMES2[x4];

        out[i  ] = (y1 << 4) | y2;
        out[i+1] = (y3 << 4) | y4;
    }
}

// ---------------------------------------------------------------------

static const uint8_t TWEAK_SHIFTS_AFTER_SBOX[KIASU_SBOXLEN] = {
    0x0, 0x8, 0x1, 0x8, 0x2, 0x2, 0x1, 0x1, 0x4, 0x4, 0x6, 0x1, 0x2, 0x2, 0x1, 0x1 
};

void prepare_delta(block out, const uint8_t tweak_shift) {
    const uint8_t tweak_shift_after_sbox = 
        TWEAK_SHIFTS_AFTER_SBOX[(tweak_shift >> 4) & 0xF] << 4;

    block DELTA_IN = { 
        tweak_shift_after_sbox, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 
    };
    mixcolumns(out, DELTA_IN);
    out[0] ^= tweak_shift;
}

// ---------------------------------------------------------------------

static inline void shiftrows(block out, const block in) {
    out[0] = (in[0] & 0xF0) | (in[2] & 0x0F);
    out[1] = (in[5] & 0xF0) | (in[7] & 0x0F);
    out[2] = (in[2] & 0xF0) | (in[4] & 0x0F);
    out[3] = (in[7] & 0xF0) | (in[1] & 0x0F);

    out[4] = (in[4] & 0xF0) | (in[6] & 0x0F);
    out[5] = (in[1] & 0xF0) | (in[3] & 0x0F);
    out[6] = (in[6] & 0xF0) | (in[0] & 0x0F);
    out[7] = (in[3] & 0xF0) | (in[5] & 0x0F);
}

// ---------------------------------------------------------------------

static inline void subbytes(block inout) {
    for (uint32_t i = 0; i < KIASU_BLOCKLEN; ++i) {
        inout[i] = (SBOX[(inout[i] >> 4) & 0xF] << 4)
                 | (SBOX[inout[i] & 0xF]);
    }
}

// ---------------------------------------------------------------------

static inline void xor_block(block out, 
                             const block a, 
                             const block b) {
    for (uint32_t i = 0; i < KIASU_BLOCKLEN; ++i) {
        out[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------

void kiasu_encrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key) {
    block tmp;

    memcpy(out, in, KIASU_BLOCKLEN);
    subbytes(out);
    shiftrows(tmp, out);
    mixcolumns(out, tmp);
    xor_block(out, out, round_key);
    xor_block(out, out, tweak);
}

// ---------------------------------------------------------------------

void kiasu_encrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds) {
    // If num_rounds != KIASU_ROUNDS, we start directly from the colliding round
    // where we do not want to break the collision from XORing the tweak to the
    // state before the first round.
    uint64_t state;
    size_t i;

    if (num_rounds == KIASU_ROUNDS) {
        xor_block(out, in, ctx->key[0]);
        xor_block(out, out, tweak);
        i = 1;
    } else {
        xor_block(out, in, ctx->key[1]);
        i = 2;
    }

    for (; i <= KIASU_ROUNDS; ++i) {
        state = Te1[out[0]] ^ Te2[out[1]] ^ Te3[out[2]] ^ Te4[out[3]]
            ^ Te5[out[4]] ^ Te6[out[5]] ^ Te7[out[6]] ^ Te8[out[7]];

        out[0] = ((state >> 56) & 0xFF) ^ ctx->key[i][0] ^ tweak[0];
        out[1] = ((state >> 48) & 0xFF) ^ ctx->key[i][1] ^ tweak[1];
        out[2] = ((state >> 40) & 0xFF) ^ ctx->key[i][2] ^ tweak[2];
        out[3] = ((state >> 32) & 0xFF) ^ ctx->key[i][3] ^ tweak[3];

        out[4] = ((state >> 24) & 0xFF) ^ ctx->key[i][4] ^ tweak[4];
        out[5] = ((state >> 16) & 0xFF) ^ ctx->key[i][5] ^ tweak[5];
        out[6] = ((state >>  8) & 0xFF) ^ ctx->key[i][6] ^ tweak[6];
        out[7] =  (state        & 0xFF) ^ ctx->key[i][7] ^ tweak[7];
    }
}

// ---------------------------------------------------------------------

void kiasu_decrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key) {
    block tmp;

    xor_block(out, in, round_key);
    xor_block(out, out, tweak);
    invert_mixcolumns(tmp, out);
    invert_shiftrows(out, tmp);
    invert_subbytes(out);
}

// ---------------------------------------------------------------------

void kiasu_decrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds) {
    const uint8_t first_round = KIASU_ROUNDS - num_rounds;
    block tmp;

    xor_block(out, in, ctx->key[KIASU_ROUNDS]);
    xor_block(out, out, tweak);

    invert_mixcolumns(tmp, out);
    invert_shiftrows(out, tmp);
    invert_subbytes(out);

    for (uint32_t i = KIASU_ROUNDS-1; i > first_round; --i) {
        xor_block(out, out, ctx->key[i]);
        xor_block(out, out, tweak);
        invert_mixcolumns(tmp, out);
        invert_shiftrows(out, tmp);
        invert_subbytes(out);
    }

    xor_block(out, out, ctx->key[first_round]);
    xor_block(out, out, tweak);
}

// ---------------------------------------------------------------------

static void expand_subkey(uint8_t* out, const uint8_t* in) {
    const uint8_t x12 = SBOX[ in[6]       & 0x0F] << 4;
    const uint8_t x13 = SBOX[(in[7] >> 4) & 0x0F];
    const uint8_t x14 = SBOX[ in[7]       & 0x0F] << 4;
    const uint8_t x15 = SBOX[(in[6] >> 4) & 0x0F];

    const uint8_t x6 = (x12 & 0xF0) | (x13 & 0x0F);
    const uint8_t x7 = (x14 & 0xF0) | (x15 & 0x0F);

    out[0] = in[0] ^ x6;
    out[1] = in[1] ^ x7;

    out[2] = in[2] ^ out[0];
    out[3] = in[3] ^ out[1];
    out[4] = in[4] ^ out[2];
    out[5] = in[5] ^ out[3];
    out[6] = in[6] ^ out[4];
    out[7] = in[7] ^ out[5];
}

// ---------------------------------------------------------------------

void kiasu_key_setup(kiasu_ctx* ctx, const block key) {
    uint8_t* current_key = (uint8_t*)ctx->key;
    uint8_t* previous_key;
    memcpy(current_key, key, KIASU_BLOCKLEN);

    for (int i = 1; i < KIASU_KEYS; ++i) {
        previous_key = current_key;
        current_key += KIASU_BLOCKLEN;
        expand_subkey(current_key, previous_key);
    }

    precompute_trail_table(sbox_trails);
}

// ---------------------------------------------------------------------

void get_sbox_trails(uint8_t keys[KIASU_SBOXLEN], 
                     const uint8_t delta_x, 
                     const uint8_t delta_y) {
    for (size_t i = 0; i < KIASU_SBOXLEN; ++i) {
        keys[i] = sbox_trails[delta_x & 0xF][delta_y & 0xF][i];
    }
}
