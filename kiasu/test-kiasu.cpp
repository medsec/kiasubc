// ===================================================================
// @last-modified 2016-06-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "kiasu.h"

// ---------------------------------------------------------------------
// Printing
// ---------------------------------------------------------------------

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
// Tests
// ---------------------------------------------------------------------

static void run_tests(const size_t num_texts, 
                      const block k, 
                      const block* expected_p, 
                      const block* tweaks,
                      const block* expected_c) {
    kiasu_ctx ctx;
    kiasu_key_setup(&ctx, k);
    int num_errors = 0;

    for (size_t i = 0; i < num_texts; ++i) {
        block c;
        kiasu_encrypt(&ctx, c, expected_p[i], tweaks[i], KIASU_ROUNDS);

        if (memcmp(expected_c[i], c, KIASU_BLOCKLEN)) {
            printf("Kiasu encryption test %zu failed\n", i);
            print_block("Expected", expected_c[i]);
            print_block("But was ", c);
            num_errors++;
        }
    }

    printf("%zu/%zu encryption tests correct\n", num_texts-num_errors, num_texts);
    num_errors = 0;

    for (size_t i = 0; i < num_texts; ++i) {
        block p;
        kiasu_decrypt(&ctx, p, expected_c[i], tweaks[i], KIASU_ROUNDS);
        
        if (memcmp(expected_p[i], p, KIASU_BLOCKLEN)) {
            printf("Kiasu decryption test %zu failed\n", i);
            print_block("Expected", expected_p[i]);
            print_block("But was ", p);
            num_errors++;
        }
    }

    printf("%zu/%zu decryption tests correct\n", num_texts-num_errors, num_texts);
}

// ---------------------------------------------------------------------

static void test_kiasu() {
    const block k = { 
        0x7d, 0x89, 0x55, 0x51, 0xa4, 0x5f, 0xa7, 0xce, 0x80, 0x25, 0x56, 0xe5, 0xff, 0x76, 0xf1, 0xcf
    };

#define NUM_KIASU_TEST_TEXTS 4
    const block expected_p[NUM_KIASU_TEST_TEXTS] = { 
        { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
        { 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF }, 
        { 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE, 0xFE }, 
        { 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0, 0xE0 }
    };

    const block expected_c[NUM_KIASU_TEST_TEXTS] = { 
        { 0xec, 0x7a, 0xfc, 0xce, 0x2b, 0x21, 0x39, 0xee, 0x82, 0x10, 0xda, 0x8f, 0x51, 0x98, 0xcb, 0x60 }, 
        { 0x21, 0xf4, 0x95, 0x6a, 0x44, 0x88, 0x54, 0xe8, 0x4c, 0x64, 0xf4, 0x42, 0xbb, 0x98, 0x41, 0xd5 }, 
        { 0x41, 0x12, 0xa0, 0xe7, 0x1d, 0x22, 0x4e, 0xa0, 0x9e, 0x6e, 0x60, 0xa6, 0x16, 0xec, 0xfe, 0x50 }, 
        { 0xe2, 0x9e, 0xbd, 0x96, 0x38, 0x5b, 0x32, 0xf9, 0x8e, 0xdb, 0xeb, 0xc8, 0x24, 0xc5, 0x63, 0x83 }
    };

    const block tweaks[NUM_KIASU_TEST_TEXTS] = { 
        { 0xe3, 0x23, 0x00, 0x00, 0xe0, 0x11, 0x00, 0x00, 0x84, 0x6d, 0x00, 0x00, 0x6f, 0x87, 0x00, 0x00 }, 
        { 0xe3, 0x23, 0x00, 0x00, 0xe0, 0x11, 0x00, 0x00, 0x84, 0x6d, 0x00, 0x00, 0x6f, 0x87, 0x00, 0x00 }, 
        { 0xe3, 0x23, 0x00, 0x00, 0xe0, 0x11, 0x00, 0x00, 0x84, 0x6d, 0x00, 0x00, 0x6f, 0x87, 0x00, 0x00 }, 
        { 0xe3, 0x23, 0x00, 0x00, 0xe0, 0x11, 0x00, 0x00, 0x84, 0x6d, 0x00, 0x00, 0x6f, 0x87, 0x00, 0x00 }
    };

    run_tests(NUM_KIASU_TEST_TEXTS, k, expected_p, tweaks, expected_c);
#undef NUM_KIASU_TEST_TEXTS
}

// ---------------------------------------------------------------------

int main() {
    test_kiasu();
    return 0;
}
