// ===================================================================
// @last-modified 2016-09-01
// This is free and unencumbered software released into the public domain.
// ===================================================================

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>

#include "kiasu.h"

// ---------------------------------------------------------------------
// Static Vars and Typedefs
// ---------------------------------------------------------------------

#ifdef DEBUG
static const block BASE_KEY = { 
    0x7d, 0x89, 0x55, 0x51, 0xa4, 0x5f, 0xa7, 0xce
};
#endif
static const block BASE_PLAINTEXT = { 
    0x0b, 0xfc, 0x8f, 0x5a, 0x35, 0xe7, 0x8f, 0x2d
};
static const block BASE_TWEAK = { 
    0x3a, 0x00, 0x4a, 0x00, 0xf1, 0x00, 0xd8, 0x00 
};
// The delta difference at the end of our boomerang trails.
static const uint64_t MASKS[7] = {
    0L, 
    0L, 
    0L, 
    0x0000FFFFFFFFFFFFL, 
    0x0FFFFFFFFFFFFFFFL,
    0x0FFFFFFFFFFFFFFFL,  
    0x0FFFF0FFFF0FFFF0L
};
static const uint8_t TWEAK_SHIFT = 0xa0;
static const uint8_t NUM_QUARTETS = 16;
static block DELTA;

typedef struct {
    block p;
    block q;
    uint64_t q64;
    block t;
    block c;
    block d;
} pair;

// ---------------------------------------------------------------------

typedef struct {
    block p;
    block p_;
    block q;
    block q_;
    block t;
    block t_;
    block c;
    block c_;
    block d;
    block d_;
} quartet;

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

static void xor_block(block out, 
                      const block a, 
                      const block b) {
    for (uint32_t i = 0; i < KIASU_BLOCKLEN; ++i) {
        out[i] = a[i] ^ b[i];
    }
}

// ---------------------------------------------------------------------

static void get_random(void* buffer, size_t num_bytes) {
    FILE* fp;
    fp = fopen("/dev/urandom", "r");
    fread(buffer, 1, num_bytes, fp);
    fclose(fp);
}

// ---------------------------------------------------------------------

static inline uint64_t to_uint64(const uint8_t* x) {
    return ((uint64_t)x[0] << 56) 
         | ((uint64_t)x[1] << 48)
         | ((uint64_t)x[2] << 40)
         | ((uint64_t)x[3] << 32)
         | ((uint64_t)x[4] << 24)
         | ((uint64_t)x[5] << 16)
         | ((uint64_t)x[6] <<  8)
         |  (uint64_t)x[7];
}

// ---------------------------------------------------------------------

static inline int compare_plaintexts(const pair& lhs, 
                                     const pair& rhs) {
    // Return -1 if B < A
    // Return  0 if B = A
    // Return +1 if B > A
    // Compares texts for the 48 bits that are NOT on the main diagonal.
    // It actually does not matter how we sort since our interest limits to 
    // the identification of pairs which are equal in all 48 bits NOT on the
    // main diagonal for finding quartets.
    return lhs.q64 < rhs.q64;
}

// ---------------------------------------------------------------------

static inline int are_equal(const pair* pair1, 
                            const pair* pair2) {
    return pair1->q64 == pair2->q64;
}

// ---------------------------------------------------------------------

static void create_plaintext(block p, 
                             const uint32_t i, 
                             const uint8_t num_rounds) {
    if (num_rounds == KIASU_ROUNDS) {
        // ---------------------------------------------------------------
        // For each plaintext, encode the counter into Nibbles 0,5,10,15.
        // ---------------------------------------------------------------

        p[0] = ((i >> 8) & 0xF0) | (p[0] & 0x0F);
        p[2] = ((i >> 8) & 0x0F) | (p[2] & 0xF0);
        p[5] = (i & 0xF0) | (p[5] & 0x0F);
        p[7] = (i & 0x0F) | (p[7] & 0xF0);
    } else {
        // ---------------------------------------------------------------
        // For each plaintext, encode the counter into the first column.
        // ---------------------------------------------------------------

        p[0] = (i >> 8) & 0xFF;
        p[1] = i & 0xFF;
    }
}

// ---------------------------------------------------------------------

static void create_pair(pair* target, 
                        const block p, 
                        const block c, 
                        const block q, 
                        const block d, 
                        const block t, 
                        const uint64_t mask) {
    memcpy(target->p, p, KIASU_BLOCKLEN);
    memcpy(target->c, c, KIASU_BLOCKLEN);
    memcpy(target->q, q, KIASU_BLOCKLEN);
    memcpy(target->d, d, KIASU_BLOCKLEN);
    memcpy(target->t, t, KIASU_BLOCKLEN);
    target->q64 = to_uint64(q) & mask;
}

// ---------------------------------------------------------------------

static void create_quartet(quartet* target, 
                           const pair* pair1, 
                           const pair* pair2) {
    memcpy(target->p, pair1->p, KIASU_BLOCKLEN);
    memcpy(target->c, pair1->c, KIASU_BLOCKLEN);
    memcpy(target->q, pair1->q, KIASU_BLOCKLEN);
    memcpy(target->d, pair1->d, KIASU_BLOCKLEN);
    memcpy(target->t, pair1->t, KIASU_BLOCKLEN);

    memcpy(target->p_, pair2->p, KIASU_BLOCKLEN);
    memcpy(target->c_, pair2->c, KIASU_BLOCKLEN);
    memcpy(target->q_, pair2->q, KIASU_BLOCKLEN);
    memcpy(target->d_, pair2->d, KIASU_BLOCKLEN);
    memcpy(target->t_, pair2->t, KIASU_BLOCKLEN);
}

// ---------------------------------------------------------------------

static void find_num_quartets(std::vector<pair>& pairs, 
                              const uint32_t num_pairs, 
                              uint32_t* num_quartets, 
                              quartet* quartets) {
    uint32_t num_equal = 0;
    pair* current;
    pair* previous;

    for (uint32_t i = 1; i < num_pairs; ++i) {
        previous = &(pairs[i-1]);
        current = &(pairs[i]);

        if (are_equal(current, previous)) {
            num_equal++;

            if (num_equal + *num_quartets < NUM_QUARTETS) {
                quartet quart;
                create_quartet(&quart, current, previous);
                quartets[num_equal + *num_quartets] = quart;
#ifdef DEBUG
                // print_quartet(current, previous);
#endif
            }
        } else if (num_equal > 0) {
            *num_quartets += ((num_equal + 1) * num_equal) / 2;
            num_equal = 0;
        }
    }

    if (num_equal > 0) {
        *num_quartets += ((num_equal + 1) * num_equal) / 2;
        num_equal = 0;
    }
}

// ---------------------------------------------------------------------

static void find_quartets(std::vector<pair>& pairs, 
                          const uint32_t num_pairs, 
                          uint32_t* num_quartets, 
                          quartet* quartets) {
    std::sort(pairs.begin(), pairs.end(), compare_plaintexts);
    find_num_quartets(pairs, num_pairs, num_quartets, quartets);
}

// ---------------------------------------------------------------------

static void find_key_byte(const size_t index, 
                          const int is_msb_nibble, 
                          const block delta_x, 
                          const block delta_y, 
                          const size_t key_index,
                          const block p,
                          const block t,
                          uint8_t key_candidates[4][KIASU_SBOXLEN]) {
    uint8_t keys[KIASU_SBOXLEN];
    uint8_t key_candidate;

    const uint8_t dx = (delta_x[index] >> (is_msb_nibble ? 4 : 0)) & 0xF;
    const uint8_t dy = (delta_y[index] >> (is_msb_nibble ? 4 : 0)) & 0xF;

    get_sbox_trails(keys, dx, dy);

    for (size_t i = 0; i < KIASU_SBOXLEN; ++i) {
        if (keys[i] != 0) {
            key_candidate = ((p[index] ^ t[index]) 
                >> (is_msb_nibble ? 4 : 0));
            key_candidate = (key_candidate ^ i) & 0xF;
            key_candidates[key_index][key_candidate]++;
        }
    }
}

// ---------------------------------------------------------------------

static void find_keys_for_quartet(const quartet* quart, 
                                  uint8_t key_candidates[4][KIASU_SBOXLEN]) {
    // Given a quartet (P, P', Q, Q') and their tweaks (T_P, T_{P'}, T_Q,
    // T_{Q'}), we assume that the pairs (P, T_P) and (P', T_{P'}) collide 
    // after the first round: Round(P, T) = Round(P', T').

    // Similarly, we assume that the pairs (Q, T_Q) and (Q', T_{Q'}) collide 
    // after the first round: Round(Q, T) = Round(Q', T').

    // We compute the difference before (Delta X) and after (Delta Y) SubBytes:
    // AddKey -> AddTweak := Delta X = X xor X' 
    //                               = (P xor P') xor (T_P xor T_{P'})
    // Y xor Y' = Delta Y <- ShiftRows^{-1} <- MixColumns^{-1} <- AddTweak^{-1}

    // Next, we show which possible values X, X' can be mapped to Y, Y' 
    // We derive the keys that map P -> X and P' -> X'
    // We repeat this procedure for Q -> X and Q' -> X'

    // For each possible value of a key byte, we increment a counter
    // Finally, we output the key-byte counters

    block tmp;
    block delta_x;
    xor_block(delta_x, quart->p, quart->p_);

    block delta_y;
    xor_block(delta_y, quart->t, quart->t_);
    xor_block(delta_x, delta_x, delta_y);
    
    invert_mixcolumns(tmp, delta_y);
    invert_shiftrows(delta_y, tmp);

#ifdef DEBUG
    print_block("P,P': delta_x", delta_x);
    print_block("P,P': delta_y", delta_y);
#endif

    find_key_byte(0, 1, delta_x, delta_y, 0, quart->p, quart->t, key_candidates); 
    find_key_byte(2, 0, delta_x, delta_y, 1, quart->p, quart->t, key_candidates); 
    find_key_byte(5, 1, delta_x, delta_y, 2, quart->p, quart->t, key_candidates); 
    find_key_byte(7, 0, delta_x, delta_y, 3, quart->p, quart->t, key_candidates); 

    // We repeat this for Q -> X and Q' -> X'
    // Note: We construct pairs so that T_Q = T_P xor TWEAK_SHIFT. 
    // Hence, it also holds that T_{Q'} = T_{P'} xor TWEAK_SHIFT. 
    // It follows that T_Q xor T_{Q'} = T_P xor T_{P'}.
    xor_block(delta_x, quart->q, quart->q_);
    xor_block(delta_y, quart->t, quart->t_);
    xor_block(delta_x, delta_x, delta_y);
    
    invert_mixcolumns(tmp, delta_y);
    invert_shiftrows(delta_y, tmp);

#ifdef DEBUG
    print_block("Q,Q': delta_x", delta_x);
    print_block("Q,Q': delta_y", delta_y);
#endif

    block t_q;
    memcpy(t_q, quart->t, KIASU_BLOCKLEN);
    t_q[0] ^= TWEAK_SHIFT;

    find_key_byte(0, 1, delta_x, delta_y, 0, quart->q, t_q, key_candidates); 
    find_key_byte(2, 0, delta_x, delta_y, 1, quart->q, t_q, key_candidates); 
    find_key_byte(5, 1, delta_x, delta_y, 2, quart->q, t_q, key_candidates); 
    find_key_byte(7, 0, delta_x, delta_y, 3, quart->q, t_q, key_candidates); 
}

// ---------------------------------------------------------------------

static uint8_t find_max(const uint8_t* array, const size_t n) {
    uint8_t max = 0;
    size_t index = 0;

    for (size_t i = 0; i < n; ++i) {
        if (max < array[i]) {
            max = array[i];
            index = i;
        }
    }

    return index;
}

// ---------------------------------------------------------------------

static void find_keys(const quartet* quartets, 
                      const uint32_t num_quartets, 
                      block recovered_key) {
    uint8_t key_candidates[4][KIASU_SBOXLEN];

    for (uint32_t i = 0; i < 4; ++i) {
        for (uint32_t j = 0; j < KIASU_SBOXLEN; ++j) {
            key_candidates[i][j] = 0;
        }
    }

    for (uint32_t i = 0; i < num_quartets; ++i) {
        find_keys_for_quartet(&(quartets[i]), key_candidates);
    }

    uint8_t max[4];

    for (int i = 0; i < 4; ++i) {
        max[i] = find_max(&(key_candidates[i][0]), KIASU_SBOXLEN);
    }

    memset(recovered_key, 0, KIASU_KEYLEN);
    recovered_key[0] = (max[0] << 4) & 0xF0;
    recovered_key[2] =  max[1]       & 0x0F;
    recovered_key[5] = (max[2] << 4) & 0xF0;
    recovered_key[7] =  max[3]       & 0x0F;

    for (uint32_t i = 0; i < 4; ++i) {
        printf("Key nibble %2d \n", i);

        for (uint32_t j = 0; j < KIASU_SBOXLEN; ++j) {
            printf("%4d", j);
        }

        puts("");

        for (uint32_t j = 0; j < KIASU_SBOXLEN; ++j) {
            printf("%4d", key_candidates[i][j]);
        }

        puts("");
    }
}

// ---------------------------------------------------------------------

static void get_key(block key) {
#ifdef DEBUG
    memcpy(key, BASE_KEY, KIASU_KEYLEN);
#else
    // Initialize Kiasu-BC with a random key.
    get_random(key, KIASU_KEYLEN);
#endif
}

// ---------------------------------------------------------------------

static void perform_experiment(block key, 
                               const uint32_t num_sets, 
                               const uint32_t num_texts_per_set, 
                               const uint32_t num_structures,
                               const uint8_t num_rounds, 
                               uint32_t* num_quartets, 
                               quartet* quartets) {
    kiasu_ctx ctx;
    kiasu_key_setup(&ctx, key);

    block t;
    block shifted_tweak;
    memcpy(t, BASE_TWEAK, KIASU_BLOCKLEN);
    memcpy(shifted_tweak, BASE_TWEAK, KIASU_BLOCKLEN);

    puts("------------------");
    print_block("Base K ", key);
    print_block("Base P ", BASE_PLAINTEXT);
    print_block("Base T ", t);
    print_block("DELTA  ", DELTA);

    const uint32_t num_pairs = num_sets * num_texts_per_set;
    std::vector<pair> pairs;
    const uint64_t mask = MASKS[num_rounds];

    // ---------------------------------------------------------------
    // Build texts
    // ---------------------------------------------------------------

    for (uint32_t k = 0; k < num_structures; ++k) {
        // ---------------------------------------------------------------
        // Setup base plaintext: random for each structure
        // ---------------------------------------------------------------

        block p;
        block q;
        block c;
        block d;
        get_random(p, KIASU_BLOCKLEN);
        
        for (uint32_t i = 0; i < num_texts_per_set; ++i) {
            create_plaintext(p, i, num_rounds);

            for (uint32_t j = 0; j < num_sets; ++j) { // For each tweak
                t[0] = ((j << 4) & 0xF0) | (t[0] & 0x0F);
                shifted_tweak[0] = t[0] ^ TWEAK_SHIFT;

                kiasu_encrypt(&ctx, c, p, t, num_rounds);
                xor_block(d, c, DELTA);
                kiasu_decrypt(&ctx, q, d, shifted_tweak, num_rounds);

                pair target;
                create_pair(&target, p, c, q, d, t, mask);
                pairs.push_back(target);
#ifdef DEBUG
                print_block("p ", p);
                print_block("q ", target.q);
                print_block("q_", (uint8_t*)&(target.q64));
                print_block("c ", c);
                print_block("d ", d);
                print_block("t ", t);
                print_block("t'", shifted_tweak);
                puts("------------------");
#endif  
            }
        }

        // ---------------------------------------------------------------
        // Correct quartets must be in the same structure. 
        // For num_rounds < KIASU_ROUNDS, we could also perform quartet 
        // search for each plaintext separately (since for those, we start 
        // from the end of Round 1 where P = P'). 
        // For num_rounds = KIASU_ROUNDS, we start from begin of Round 1.
        // ---------------------------------------------------------------

        find_quartets(pairs, num_pairs, num_quartets, quartets);

        printf("Structure %2u/%2u: %u quartets \n", 
            k+1, num_structures, *num_quartets);
    }
}

// ---------------------------------------------------------------------

static void perform_experiments(const uint32_t num_keys, 
                                const uint32_t num_texts_per_set, 
                                const uint32_t num_sets, 
                                const uint32_t num_structures,
                                const uint8_t num_rounds) {
    block key;
    
    for (uint32_t i = 0; i < num_keys; ++i) {
        get_key(key);
        
        uint32_t num_quartets = 0;
        quartet quartets[NUM_QUARTETS];

        perform_experiment(
            key, num_sets, num_texts_per_set, 
            num_structures, num_rounds, &num_quartets, quartets
        );

        puts("Experiment finished");

        // -----------------------------------------------------------------
        // Append key recovery for full #rounds
        // -----------------------------------------------------------------
        
        if (num_quartets > NUM_QUARTETS) {
            num_quartets = NUM_QUARTETS;
        }

        if ((num_rounds == KIASU_ROUNDS) && (num_quartets > 0)) {
            block recovered_key;
            find_keys(quartets, num_quartets, recovered_key);
            print_block("Found K", recovered_key);
        }

        print_block("Base K ", key);
    }
}

// ---------------------------------------------------------------------
// User interface, main, and options.
// ---------------------------------------------------------------------

static void print_usage(const char* name) {
    puts("Tests the 6-round boomerang distinguisher for Kiasu-BC"\
         "with 64-bit Mini-Kiasu-BC. Mini-Kiasu-BC uses the AES-128 "\
         "key schedule on nibbles, the Piccolo S-box and the AES MDS matrix "\
         "with x^4 + x + 1 as irreducible polynomial.");
    printf("Usage: %s <#keys> <#texts per set> <#sets> [<#rounds>]\n", name);
    puts(" [<#structures>");
    puts("- <#keys>: #Random keys to test. E.g. 1 or 20");
    puts("- <#texts per set>: #Texts per set. E.g. 2^{16}");
    puts("- <#sets>: #Sets. Should be <= 2^4");
    puts("- <#rounds>: #Rounds of the distinguisher. Should be 1-6. Default: 5");
    puts("- <#structures>: #Structures. Default: 1");
}

// ---------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        exit(1);
    }

    prepare_delta(DELTA, TWEAK_SHIFT);

    // ---------------------------------------------------------------------
    // Interface
    // ---------------------------------------------------------------------

    const uint32_t num_keys = atoi(argv[1]);
    const uint32_t num_texts_per_set = atoi(argv[2]);
    const uint32_t num_sets = atoi(argv[3]);
    uint32_t num_structures = 1;
    uint8_t num_rounds = 5;

    if (argc >= 5) {
        num_rounds = (uint8_t)(atoi(argv[4]) & 0xFF);
    }

    if (argc >= 6) {
        num_structures = atoi(argv[5]);
    }

    perform_experiments(
        num_keys, num_texts_per_set, num_sets, num_structures, num_rounds
    );
    
    return 0;
}
