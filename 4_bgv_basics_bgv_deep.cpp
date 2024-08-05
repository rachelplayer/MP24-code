// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "examples.h"

using namespace std;
using namespace seal;

void example_bgv_basics()
{
    print_example_banner("Example: BGV Basics");

    /* Set number of trials. */
    int trials = 10000;

    /* Set verbose to true for debugging. */
    bool verbose = false;

    /* Select parameters appropriate for our experiment:
       n < 16384 too small to support computation. */
    EncryptionParameters parms(scheme_type::bgv);

    size_t poly_modulus_degree = 16384;
    parms.set_poly_modulus_degree(poly_modulus_degree);

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
    Plaintext plain3;
    Plaintext plain4;
    Plaintext plain5;
    Plaintext plain6;
    Plaintext plain7;
    Plaintext plain8;
    Ciphertext encrypted1;
    Ciphertext encrypted2;
    Ciphertext encrypted3;
    Ciphertext encrypted4;
    Ciphertext encrypted5;
    Ciphertext encrypted6;
    Ciphertext encrypted7;
    Ciphertext encrypted8;
    Ciphertext encrypted9;
    Ciphertext encrypted10;
    Ciphertext encrypted11;
    Ciphertext encrypted12;
    Ciphertext encrypted13;
    Ciphertext encrypted14;
    Ciphertext encrypted15;

    /* Holders for the running total of the observed noises in ciphertexts */
    double total_fresh_observed(0);
    double total_mult1_observed(0);
    double total_mult2_observed(0);
    double total_mult3_observed(0);

    /* Gather data */
    for (int i = 0; i < trials; i++)
    {

         /*
         Here we create the input plaintext matrices 
         encrypting i+1, ...., i+8 respectively in the first slot.
         */

         vector<uint64_t> pod_matrix1(slot_count, 0ULL);
         pod_matrix1[0] = i+1;
         vector<uint64_t> pod_matrix2(slot_count, 0ULL);
         pod_matrix2[0] = i+2;
         vector<uint64_t> pod_matrix3(slot_count, 0ULL);
         pod_matrix3[0] = i+3;
         vector<uint64_t> pod_matrix4(slot_count, 0ULL);
         pod_matrix4[0] = i+4;
         vector<uint64_t> pod_matrix5(slot_count, 0ULL);
         pod_matrix5[0] = i+5;
         vector<uint64_t> pod_matrix6(slot_count, 0ULL);
         pod_matrix6[0] = i+6;
         vector<uint64_t> pod_matrix7(slot_count, 0ULL);
         pod_matrix7[0] = i+7;
         vector<uint64_t> pod_matrix8(slot_count, 0ULL);
         pod_matrix8[0] = i+8;

         /* Encode the matrices into plaintexts. */
         batch_encoder.encode(pod_matrix1, plain1);
         batch_encoder.encode(pod_matrix2, plain2);
         batch_encoder.encode(pod_matrix3, plain3);
         batch_encoder.encode(pod_matrix4, plain4);
         batch_encoder.encode(pod_matrix5, plain5);
         batch_encoder.encode(pod_matrix6, plain6);
         batch_encoder.encode(pod_matrix7, plain7);
         batch_encoder.encode(pod_matrix8, plain8);

         /* Encrypt the plaintexts into ciphertexts */
         encryptor.encrypt(plain1, encrypted1);
         encryptor.encrypt(plain2, encrypted2);
         encryptor.encrypt(plain3, encrypted3);
         encryptor.encrypt(plain4, encrypted4);
         encryptor.encrypt(plain5, encrypted5);
         encryptor.encrypt(plain6, encrypted6);
         encryptor.encrypt(plain7, encrypted7);
         encryptor.encrypt(plain8, encrypted8);

         /* What is the noise growth after fresh encryption? */
         auto fresh_noise = decryptor.invariant_noise_budget(encrypted1);
         total_fresh_observed += fresh_noise;

        /*  Multiply the ciphertexts pairwise and store the output in encrypted9, ... , encrypted12 */
         evaluator.multiply(encrypted1, encrypted2, encrypted9);
         evaluator.multiply(encrypted3, encrypted4, encrypted10);
         evaluator.multiply(encrypted5, encrypted6, encrypted11);
         evaluator.multiply(encrypted7, encrypted8, encrypted12);

         /* What is the noise growth after first multiplication? */
         auto mult1_noise = decryptor.invariant_noise_budget(encrypted9);
         total_mult1_observed += mult1_noise;

        /* Relinearize */
        evaluator.relinearize_inplace(encrypted9, relin_keys); 
        evaluator.relinearize_inplace(encrypted10, relin_keys); 
        evaluator.relinearize_inplace(encrypted11, relin_keys); 
        evaluator.relinearize_inplace(encrypted12, relin_keys); 

        /*  Multiply the ciphertexts pairwise and store the output in encrypted13, encrypted14 */
         evaluator.multiply(encrypted9, encrypted10, encrypted13);
         evaluator.multiply(encrypted11, encrypted12, encrypted14); 

         /* What is the noise growth after second multiplication? */
         auto mult2_noise = decryptor.invariant_noise_budget(encrypted13);
         total_mult2_observed += mult2_noise;

        /* Relinearize */
        evaluator.relinearize_inplace(encrypted13, relin_keys); 
        evaluator.relinearize_inplace(encrypted14, relin_keys); 

        /*  Multiply the ciphertexts encrypted13 and encrypted14 and stored in encrypted15 */
         evaluator.multiply(encrypted13, encrypted14, encrypted15);

         /* What is the noise growth after third multiplication? */
         auto mult3_noise = decryptor.invariant_noise_budget(encrypted15);
         total_mult3_observed += mult3_noise;
    }

    /* Debugging: check that decryption is correct. */
    if(verbose)
    {
        cout << "Check correctness:" << endl;
        Plaintext decrypted_result;
        decryptor.decrypt(encrypted15, decrypted_result);
        vector<uint64_t> pod_result;
        batch_encoder.decode(decrypted_result, pod_result);
        print_matrix(pod_result, row_size);
    }

    /* Compute the mean of the observed noises */
    auto mean_fresh_observed = total_fresh_observed / trials;
    auto mean_mult1_observed = total_mult1_observed / trials;
    auto mean_mult2_observed = total_mult2_observed / trials;
    auto mean_mult3_observed = total_mult3_observed / trials;

    /* Print out the results */
    cout << "After fresh encryption:" << endl;
    cout << "Mean noise budget observed: " << mean_fresh_observed  << endl;    
    cout << endl;

    cout << "After first multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult1_observed  << endl;    
    cout << endl;

    cout << "After second multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult2_observed  << endl;    
    cout << endl;

    cout << "After third multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult3_observed  << endl;    
    cout << endl;

}
