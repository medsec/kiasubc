#pragma once

#include <stdint.h>
#include <emmintrin.h>

// ---------------------------------------------------------------------

#define KIASU_SBOXLEN  16
#define KIASU_BLOCKLEN  8
#define KIASU_ROUNDS    7
#define KIASU_KEYS      8
#define KIASU_KEYLEN    8

typedef uint8_t block[KIASU_BLOCKLEN];

typedef struct {
    block key[KIASU_KEYS];
    __m128i k[KIASU_KEYS];
} kiasu_ctx;

// ---------------------------------------------------------------------

void prepare_delta(block out, const uint8_t tweak_shift);
void invert_mixcolumns(block out, const block in);
void invert_shiftrows(block out, const block in);

// ---------------------------------------------------------------------

void kiasu_encrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key);

// ---------------------------------------------------------------------

void kiasu_encrypt_final_round(block out, 
                               const block in, 
                               const block tweak, 
                               const block round_key);

// ---------------------------------------------------------------------

void kiasu_encrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds);

// ---------------------------------------------------------------------

void kiasu_decrypt_round(block out, 
                         const block in, 
                         const block tweak, 
                         const block round_key);

// ---------------------------------------------------------------------

void kiasu_decrypt_final_round(block out, 
                               const block in, 
                               const block tweak, 
                               const block round_key);

// ---------------------------------------------------------------------

void kiasu_decrypt(const kiasu_ctx* ctx, 
                   block out, 
                   const block in, 
                   const block tweak, 
                   const uint8_t num_rounds);

// ---------------------------------------------------------------------

void kiasu_key_setup(kiasu_ctx* ctx, const block key);

// ---------------------------------------------------------------------

void get_sbox_trails(uint8_t keys[KIASU_SBOXLEN], 
                     const uint8_t delta_x, 
                     const uint8_t delta_y);

// ---------------------------------------------------------------------

void compute_table();

// void derive_d(const block klast,
//               block d, 
//               const block c, 
//               const block tweak_c, 
//               const block tweak_d, 
//               const block DELTA);
