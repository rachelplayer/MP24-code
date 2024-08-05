// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "examples.h"

using namespace std;
using namespace seal;

void example_bgv_basics()
{
    print_example_banner("Example: BGV Basics");

    /* Set number of trials. */
    int trials = 1;

    /* Set verbose to true for debugging. */
    bool verbose = false;

    /* Select parameters appropriate for our experiment */
    EncryptionParameters parms(scheme_type::bgv);

    size_t poly_modulus_degree = 4096;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    
    //size_t poly_modulus_degree = 8192;
    //parms.set_poly_modulus_degree(poly_modulus_degree);

    //size_t poly_modulus_degree = 16384;
    //parms.set_poly_modulus_degree(poly_modulus_degree);

    //size_t poly_modulus_degree = 32768;
    //parms.set_poly_modulus_degree(poly_modulus_degree);


    /*
    Use BFVDefault coeff_modulus. 
    */
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));

    /* Use same plain_modulus as used in the BGV Basics example. */
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));
    SEALContext context(parms);

    /*
    Print the parameters that we have chosen.
    */
    print_parameters(context);

    /* Also print exact plaintext modulus chosen. */
    auto &context_data = *context.key_context_data();
    std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;

    /* Generate keys */
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    /*
    Using a (quite redundant!) batch encoding
    */
    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();
    size_t row_size = slot_count / 2;

    /* Construct plaintext and ciphertext objects */
    Plaintext plain1;
    Plaintext plain2;
    Ciphertext encrypted1;
    Ciphertext encrypted2;
    Ciphertext encrypted3;
    Ciphertext encrypted4;

    /* Holders for the running total of the observed noises in ciphertexts */
    double total_fresh_observed(0);
    double total_add_observed(0);
    double total_mult_observed(0);
    double total_modswitch_observed(0);

    /* Gather data */
    for (int i = 0; i < trials; i++)
    {

         /*
         Here we create the following input plaintext matrices:
            [ i,  0,  0,  0,  0,  0, ...,  0 ]
            [ 0,  0,  0,  0,  0,  0, ...,  0 ]

            [ i+1,  0,  0,  0,  0,  0, ...,  0 ]
            [ 0,  0,  0,  0,  0,  0, ...,  0 ]
         */
         vector<uint64_t> pod_matrix1(slot_count, 0ULL);
         pod_matrix1[0] = i;
         vector<uint64_t> pod_matrix2(slot_count, 0ULL);
         pod_matrix2[0] = i+1;

         /* Encode the matrices into plaintexts. */
         batch_encoder.encode(pod_matrix1, plain1);
         batch_encoder.encode(pod_matrix2, plain2);

         /* Encrypt the plaintexts into ciphertexts */
         encryptor.encrypt(plain1, encrypted1);
         encryptor.encrypt(plain2, encrypted2);

         /* What is the noise growth after fresh encryption? */
         auto fresh_noise = decryptor.invariant_noise_budget(encrypted1);
         total_fresh_observed += fresh_noise;

         /* Add encrypted1 and encrypted2 together and store in encrypted3. */
         evaluator.add(encrypted1, encrypted2, encrypted3);

         /* What is the noise growth after addition? */
         auto add_noise = decryptor.invariant_noise_budget(encrypted3);
         total_add_observed += add_noise;

         /* Multiply encrypted3 by encrypted2 and store in encrypted4. */
         evaluator.multiply(encrypted3, encrypted2, encrypted4);

         /* What is the noise growth after multiplication? */
         auto mult_noise = decryptor.invariant_noise_budget(encrypted4);
         total_mult_observed += mult_noise;

         /* Modulus switch encrypted4 to next prime in the chain. */
        evaluator.mod_switch_to_next_inplace(encrypted4);

         /* What is the noise growth after mod switch? */
         auto modswitch_noise = decryptor.invariant_noise_budget(encrypted4);
         total_modswitch_observed += modswitch_noise;

    }

    /* Debugging: check that decryption is correct. */
    if(verbose)
    {
        cout << "Check correctness:" << endl;
        Plaintext decrypted_result;
        decryptor.decrypt(encrypted4, decrypted_result);
        vector<uint64_t> pod_result;
        batch_encoder.decode(decrypted_result, pod_result);
        print_matrix(pod_result, row_size);
    }

    /* Compute the mean of the observed noises */
    auto mean_fresh_observed = total_fresh_observed / trials;
    auto mean_add_observed = total_add_observed / trials;
    auto mean_mult_observed = total_mult_observed / trials;
    auto mean_modswitch_observed = total_modswitch_observed / trials;

    /* Print out the results */
    cout << "After fresh encryption:" << endl;
    cout << "Mean noise budget observed: " << mean_fresh_observed  << endl;    
    cout << endl;

    cout << "After addition:" << endl;
    cout << "Mean noise budget observed: " << mean_add_observed  << endl;    
    cout << endl;

    cout << "After multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult_observed  << endl;    
    cout << endl;

    cout << "After modulus switching:" << endl;
    cout << "Mean noise budget observed: " << mean_modswitch_observed  << endl;    
    cout << endl;

}
