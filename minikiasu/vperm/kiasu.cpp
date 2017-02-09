// ===================================================================
// @last-modified 2016-06-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstring>

#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>

#include "kiasu.h"
#include "kiasu-tables.h"

uint8_t SBOX[KIASU_SBOXLEN] = { 
    14,  4, 11,  2,  3,  8,  0,  9,  1, 10,  7, 15,  6, 12,  5, 13 
};
uint8_t INVERSE_SBOX[KIASU_SBOXLEN] = { 
     6,  8,  3,  4,  1, 14, 12, 10,  5,  7,  9,  2, 13, 15,  0, 11
};

static uint8_t times2[KIASU_SBOXLEN] = { 
    0x0, 0x2, 0x4, 0x6, 0x8, 0xa, 0xc, 0xe, 0x3, 0x1, 0x7, 0x5, 0xb, 0x9, 0xf, 0xd
};
static uint8_t times3[KIASU_SBOXLEN] = { 
    0x0, 0x3, 0x6, 0x5, 0xc, 0xf, 0xa, 0x9, 0xb, 0x8, 0xd, 0xe, 0x7, 0x4, 0x1, 0x2
};

// Input transform shuffles
static const __m128i PiccoloInShuffleOdd  = _mm_setr_epi8(0, 0xff, 1, 0xff, 2, 0xff, 3, 0xff, 4, 0xff, 5, 0xff, 6, 0xff, 7, 0xff);
static const __m128i PiccoloInShuffleEven = _mm_setr_epi8(0xff, 0, 0xff, 1, 0xff, 2, 0xff, 3, 0xff, 4, 0xff, 5, 0xff, 6, 0xff, 7);

// Output transform shuffles
static const __m128i PiccoloOutShuffleOdd  = _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff);
static const __m128i PiccoloOutShuffleEven = _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff);

// ShiftRows shuffles
static const __m128i SHIFT_ROWS = _mm_setr_epi8(0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11);
static const __m128i INVERSE_SHIFT_ROWS = _mm_setr_epi8(0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3);

// Piccolo SBox, SBOX*2, and SBOX*3 (high and low nibbles)
static const __m128i PiccoloSBoxL = _mm_setr_epi8(0x0e, 0x04, 0x0b, 0x02, 0x03, 0x08, 0x00, 0x09, 0x01, 0x0a, 0x07, 0x0f, 0x06, 0x0c, 0x05, 0x0d);
static const __m128i TwoMulPiccoloSBoxL = _mm_setr_epi8(0x0f, 0x08, 0x05, 0x04, 0x06, 0x03, 0x00, 0x01, 0x02, 0x07, 0x0e, 0x0d, 0x0c, 0x0b, 0x0a, 0x09);
static const __m128i ThreeMulPiccoloSBoxL = _mm_setr_epi8(0x01, 0x0c, 0x0e, 0x06, 0x05, 0x0b, 0x00, 0x08, 0x03, 0x0d, 0x09, 0x02, 0x0a, 0x07, 0x0f, 0x04);

// Piccolo inverse SBox (high and low nibbles)
static const __m128i PiccoloInvSBoxL = _mm_setr_epi8(0x6, 0x8, 0x3, 0x4, 0x1, 0xe, 0xc, 0xa, 0x5, 0x7, 0x9, 0x2, 0xd, 0xf, 0x0, 0xb);
// Piccolo Multiplication by 9, B, D, and E
// Piccolo 9 multiplication (high and low nibbles)
static const __m128i NineMulPiccoloL = _mm_setr_epi8(0x0, 0x9, 0x1, 0x8, 0x2, 0xb, 0x3, 0xa, 0x4, 0xd, 0x5, 0xc, 0x6, 0xf, 0x7, 0xe);
static const __m128i BMulPiccoloL = _mm_setr_epi8(0x0, 0xb, 0x5, 0xe, 0xa, 0x1, 0xf, 0x4, 0x7, 0xc, 0x2, 0x9, 0xd, 0x6, 0x8, 0x3);
static const __m128i DMulPiccoloL = _mm_setr_epi8(0x0, 0xd, 0x9, 0x4, 0x1, 0xc, 0x8, 0x5, 0x2, 0xf, 0xb, 0x6, 0x3, 0xe, 0xa, 0x7);
static const __m128i EMulPiccoloL = _mm_setr_epi8(0x0, 0xe, 0xf, 0x1, 0xd, 0x3, 0x2, 0xc, 0x9, 0x7, 0x6, 0x8, 0x4, 0xa, 0xb, 0x5);

// MixCoumns shuffles
static const __m128i PiccoloOneShufa  = _mm_setr_epi8(2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13);
static const __m128i PiccoloOneShufb  = _mm_setr_epi8(3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14);
static const __m128i PiccoloThreeShuf = _mm_setr_epi8(1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12);

static const __m128i PiccoloAndMaskL  = _mm_setr_epi8(0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f);
static const __m128i PiccoloAndMaskH  = _mm_setr_epi8(0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0, 0xf0);

#define load(a)       _mm_loadu_si128(a)
#define pshufb(a, b)  _mm_shuffle_epi8(a, b)
#define store(a, b)   _mm_storeu_si128((__m128i*)a, b)
#define vand(a, b)    _mm_and_si128(a, b)
#define vxor(a, b)    _mm_xor_si128(a, b)

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
// ---------------------------------------------------------------------

static void print128(const char label[], const __m128i x) {
    uint8_t y[16];
    store((__m128i*)y, x);
    print_hex(label, y, 16);
}
#endif

// ---------------------------------------------------------------------

void invert_mixcolumns(block out, const block in) {
    const uint64_t state = Td1[in[0]] 
        ^ Td2[in[1]] 
        ^ Td3[in[2]] 
        ^ Td4[in[3]]
        ^ Td5[in[4]] 
        ^ Td6[in[5]] 
        ^ Td7[in[6]] 
        ^ Td8[in[7]];
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
    for (uint8_t i = 0; i < KIASU_BLOCKLEN; ++i) {
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

        y1 = times2[x1] ^ times3[x2] ^        x3  ^        x4;
        y2 =        x1  ^ times2[x2] ^ times3[x3] ^        x4;
        y3 =        x1  ^        x2  ^ times2[x3] ^ times3[x4];
        y4 = times3[x1] ^        x2  ^        x3  ^ times2[x4];

        out[i  ] = (y1 << 4) | y2;
        out[i+1] = (y3 << 4) | y4;
    }
}

// ---------------------------------------------------------------------

// TWEAK_SHIFTS_AFTER_SBOX[alpha] = beta, where
// beta is our own chosen output differences from the Piccolo S-box such that
// Pr[alpha -> beta] = 2^{-2} (i.e. is maximal) through the S-box.
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

static __m128i format_input(const uint64_t* in) {
    // Given an eight-byte Mini-KIASU-BC state:
    // Input:           [0_h, 0_l, 1_h, 1_l, ..., 7_h, 7_l], 
    //   of  8 bytes, (0_h, 0_l) being Byte 0 and so on.
    // Expected output: [0_h, 0_l, 1_h, 1_l, ..., 7_h, 7_l], 
    //   of 16 bytes, each taking a byte (and only the low nibbles used)

    __m128i tmp1;
    __m128i tmp4;

    tmp1 = _mm_set_epi64((__m64)0L, (__m64)*in);
    tmp4 = tmp1;

    // Keep the low parts
    tmp1 = vand(tmp1, PiccoloAndMaskL);

    // Keep the high parts
    tmp4 = vand(tmp4, PiccoloAndMaskH);

    // Shift the elements. Word1 high is put low, Word2 low is put high.
    tmp4 = _mm_srli_epi64(tmp4, 4);

    // Shuffle
    tmp1 = pshufb(tmp1, PiccoloInShuffleEven);
    tmp4 = pshufb(tmp4, PiccoloInShuffleOdd);

    return vxor(tmp4, tmp1);
}

// ---------------------------------------------------------------------

static uint64_t format_output(const __m128i in) {
    // Given an eight-byte Mini-KIASU-BC state:
    // Input:           [0_h, 0_l, 1_h, 1_l, ..., 7_h, 7_l], 
    //   of 16 bytes, each taking a byte (and only the low nibbles used)
    // Expected output: [0_h, 0_l, 1_h, 1_l, ..., 7_h, 7_l], 
    //   of  8 bytes, (0_h, 0_l) being Byte 0 and so on.

    // The low nibbles
    __m128i tmpL = in; 

    // The high nibbles, shifted to the high parts of the bytes
    __m128i tmpH = _mm_slli_epi64(in, 4); 
    tmpH = vand(tmpH, PiccoloAndMaskH);

    // Shuffle
    tmpL = pshufb(tmpL, PiccoloOutShuffleEven);
    tmpH = pshufb(tmpH, PiccoloOutShuffleOdd);

    tmpL  = vxor(tmpL, tmpH);
    return _mm_extract_epi64(tmpL, 0);
}

// ---------------------------------------------------------------------

void kiasu_encrypt_final_round(block out, 
                               const block in, 
                               const block tweak, 
                               const block round_key) {
    block tmp;

    memcpy(out, in, KIASU_BLOCKLEN);
    subbytes(out);
    shiftrows(tmp, out);
    xor_block(out, tmp, round_key);
    xor_block(out, out, tweak);
}

// ---------------------------------------------------------------------

static void internal_encrypt_round(__m128i* state, 
                                  const __m128i k,
                                  const __m128i t) {
    // TODO: Integrate SHIFT_ROWS into MixColumns Shuffling
    *state = pshufb(*state, SHIFT_ROWS);

    __m128i tmp1 = pshufb(PiccoloSBoxL, *state);
    __m128i tmp2 = pshufb(TwoMulPiccoloSBoxL, *state);
    __m128i tmp3 = pshufb(ThreeMulPiccoloSBoxL, *state);
    __m128i tmp4;

    // Shuffle for MixColumns 
    // No need to shuffle the 2*S[x_i] since these terms are on the diagonal.
    tmp4 = pshufb(tmp1, PiccoloOneShufa);
    tmp1 = pshufb(tmp1, PiccoloOneShufb);
    tmp3 = pshufb(tmp3, PiccoloThreeShuf);
    
    tmp1 = vxor(tmp1, tmp2);
    tmp3 = vxor(tmp3, tmp4);
    *state = vxor(tmp1, tmp3);

    // XOR all together
    *state = vxor(*state, k);
    *state = vxor(*state, t);
}

// ---------------------------------------------------------------------

static void internal_encrypt_final_round(__m128i* state, 
                                        const __m128i k,
                                        const __m128i t) {
    *state = pshufb(*state, SHIFT_ROWS);
    *state = pshufb(PiccoloSBoxL, *state);
    *state = vxor(*state, k);
    *state = vxor(*state, t);
}

// ---------------------------------------------------------------------

void kiasu_encrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds) {
    const uint64_t* m64 = (uint64_t*)in;
    const uint64_t* t64 = (uint64_t*)tweak;

    const __m128i t = format_input(t64);
    __m128i state = format_input(m64);

    state = vxor(state, ctx->k[0]);
    state = vxor(state, t);

    for (int i = 1; i < num_rounds; ++i) {
        internal_encrypt_round(&state, ctx->k[i], t);
    }

    if (num_rounds == KIASU_ROUNDS) {
        internal_encrypt_final_round(&state, ctx->k[num_rounds], t);
    } else {
        internal_encrypt_round(&state, ctx->k[num_rounds], t);
    }

    const uint64_t result = format_output(state);
    memcpy(out, &result, KIASU_BLOCKLEN);
}

// ---------------------------------------------------------------------

void kiasu_decrypt_final_round(block out, 
                               const block in, 
                               const block tweak, 
                               const block round_key) {
    block tmp;

    xor_block(out, in, round_key);
    xor_block(tmp, out, tweak);
    invert_shiftrows(out, tmp);
    invert_subbytes(out);
}

// ---------------------------------------------------------------------

static void internal_decrypt_round(__m128i* state, 
                                   const __m128i k,
                                   const __m128i t) {
    // XOR all together
    *state = vxor(*state, k);
    *state = vxor(*state, t);

    // MixColumns multiples
    __m128i tmp1 = pshufb(NineMulPiccoloL, *state);
    __m128i tmp2 = pshufb(BMulPiccoloL, *state);
    __m128i tmp3 = pshufb(DMulPiccoloL, *state);
    __m128i tmp4 = pshufb(EMulPiccoloL, *state);

    // Shuffle for MixColumns 
    // No need to shuffle the 0xe * x_i since these terms are on the diagonal.
    tmp1 = pshufb(tmp1, PiccoloOneShufb);
    tmp2 = pshufb(tmp2, PiccoloThreeShuf);
    tmp3 = pshufb(tmp3, PiccoloOneShufa);

    tmp1 = vxor(tmp1, tmp2);
    tmp3 = vxor(tmp3, tmp4);
    tmp1 = vxor(tmp1, tmp3);

    // TODO: Integrate SHIFT_ROWS into MixColumns Shuffling
    *state = pshufb(tmp1, INVERSE_SHIFT_ROWS);
    *state = pshufb(PiccoloInvSBoxL, *state);
}

// ---------------------------------------------------------------------

static void internal_decrypt_final_round(__m128i* state, 
                                         const __m128i k,
                                         const __m128i t) {
    *state = vxor(*state, k);
    *state = vxor(*state, t);
    *state = pshufb(*state, INVERSE_SHIFT_ROWS);
    *state = pshufb(PiccoloInvSBoxL, *state);
}

// ---------------------------------------------------------------------

void kiasu_decrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds) {
    const uint8_t first_round = KIASU_ROUNDS - num_rounds;
    const uint64_t* m64 = (uint64_t*)in;
    const uint64_t* t64 = (uint64_t*)tweak;

    const __m128i t = format_input(t64);
    __m128i state = format_input(m64);

    internal_decrypt_final_round(&state, ctx->k[KIASU_ROUNDS], t);

    for (uint32_t i = KIASU_ROUNDS-1; i > first_round; --i) {
        internal_decrypt_round(&state, ctx->k[i], t);
    }

    state = vxor(state, ctx->k[first_round]);
    state = vxor(state, t);

    const uint64_t result = format_output(state);
    memcpy(out, &result, KIASU_BLOCKLEN);
}

// ---------------------------------------------------------------------

// void derive_d(const block klast,
//               block d, 
//               const block c, 
//               const block tweak_c, 
//               const block tweak_d, 
//               const block DELTA) {
//     // ---------------------------------------------------------------
//     // Note: DELTA contains already 
//     // MC(tweak_shift) xor (tweak xor shifted_tweak).
//     // Thus, we do NOT xor the tweak-difference (tweak xor shifted_tweak)
//     // from Round 6 to the state one round backwards.
//     // ---------------------------------------------------------------

//     const __m128i tc = format_input((uint64_t*)tweak_c);
//     const __m128i td = format_input((uint64_t*)tweak_d);
//     const __m128i k = format_input((uint64_t*)klast);
//     const __m128i delta = format_input((uint64_t*)DELTA);
//     __m128i sc = format_input((uint64_t*)c);
//     __m128i sd = sc;

//     internal_decrypt_final_round(&sd, k, tc);
//     sd = vxor(sd, delta);
//     internal_encrypt_final_round(&sd, k, td);

//     sc = vand(sc, _mm_setr_epi8(0x0,0xf,0xf,0xf,0xf,0xf,0xf,0x0,0xf,0xf,0x0,0xf,0xf,0x0,0xf,0xf));
//     sd = vand(sd, _mm_setr_epi8(0xf,0x0,0x0,0x0,0x0,0x0,0x0,0xf,0x0,0x0,0xf,0x0,0x0,0xf,0x0,0x0));
//     sd = vxor(sc, sd);

//     const uint64_t result = format_output(sd);
//     memcpy(d, &result, KIASU_BLOCKLEN);
// }

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
    
    ctx->k[0] = format_input((uint64_t*)current_key);

    for (int i = 1; i < KIASU_KEYS; ++i) {
        previous_key = current_key;
        current_key += KIASU_BLOCKLEN;
        expand_subkey(current_key, previous_key);

        ctx->k[i] = format_input((uint64_t*)current_key);
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
