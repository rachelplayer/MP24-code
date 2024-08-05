/* 
    HElib experiments
    Some of this code is adapted from Section 5 of https://github.com/shaih/HElib/blob/master/doc/designDocument/he-library.pdf 
    Some of this code is adapted SEAL/examples/examples.cpp at https://github.com/microsoft/SEAL commit ba2d578
    This code requires the following changes to be made to HElib:
        - make Ctxt::tensorProduct public so we can do homomorphic multiplication without automatically mod switching or relinearizing
*/

#include <iostream>
#include <iomanip>

#include <helib/helib.h>
#include <helib/binaryArith.h>
#include <helib/intraSlot.h>

//#include "EncryptedArray.h"
//#include "FHE.h"
//#include "norms.h"

using namespace std;

/*
This function computes, for a given chain of operations, over a user-specified number of trials,
an average observed noise growth in ciphertexts.
*/
void test_noise(int trials);

/* Helper functions */
NTL::xdouble get_sum_of_squared_differences(NTL::xdouble mean, vector<NTL::xdouble> array, int size_of_array);
NTL::xdouble get_standard_dev(NTL::xdouble mean, vector<NTL::xdouble> array, int trials);
NTL::xdouble get_noise();
NTL::xdouble get_noise_budget(helib::Ctxt encrypted, helib::SecKey secret_key);
NTL::xdouble get_helib_estimated_noise_budget(helib::Ctxt encrypted);

NTL::xdouble get_sum_of_squared_differences(NTL::xdouble mean, vector<NTL::xdouble> array, int size_of_array)
{
    NTL::xdouble sum_of_squared_differences(0);
    NTL::xdouble temp(0);

    for (int i=0; i<size_of_array; i++)
    {
        temp = array[i] - mean;
        temp = temp * temp;
        sum_of_squared_differences += temp;
    }
    return sum_of_squared_differences;
}

NTL::xdouble get_standard_dev(NTL::xdouble mean, vector<NTL::xdouble> array, int trials)
{
    NTL::xdouble trials_copy(trials);
    NTL::xdouble denominator = trials_copy - 1;
    NTL::xdouble sum_sq_diff = get_sum_of_squared_differences(mean, array, trials);
    NTL::xdouble variance = sum_sq_diff / denominator;
    NTL::xdouble std_dev = sqrt(variance);
    return std_dev;
}

/* Inspired by the HElib debugging function decryptAndPrint */
NTL::xdouble get_noise(helib::Ctxt encrypted, helib::SecKey secret_key)
{
    NTL::ZZX plaintext, noise_poly;
    NTL::ZZ noise_zz;
    NTL::xdouble noise;
    secret_key.Decrypt(plaintext, encrypted, noise_poly);
    noise_zz = helib::largestCoeff(noise_poly);
    conv(noise, noise_zz);
    return noise;
}

NTL::xdouble get_noise_budget(helib::Ctxt encrypted, helib::SecKey secret_key)
{
    NTL::xdouble noise, log_noise, log_q, noise_budget;
    noise = get_noise(encrypted, secret_key);
    log_noise = log(noise)/log(NTL::xdouble(2));
    log_q = NTL::xdouble(encrypted.getContext().logOfProduct(encrypted.getPrimeSet())/log(2));
    noise_budget = log_q - log_noise - 1;
    return noise_budget;
}

NTL::xdouble get_helib_estimated_noise_budget(helib::Ctxt encrypted)
{
    NTL::xdouble helib_est_noise, log_helib_est_noise, log_q, helib_est_noise_budget;
    helib_est_noise = encrypted.getNoiseBound();
    log_helib_est_noise = log(helib_est_noise)/log(NTL::xdouble(2));
    log_q = NTL::xdouble(encrypted.getContext().logOfProduct(encrypted.getPrimeSet())/log(2));
    helib_est_noise_budget = log_q - log_helib_est_noise - 1;
    return helib_est_noise_budget;
}

int main()
{

    while (true)
    {
        cout << "\n HElib noise budget experiments:" << endl << endl;
        cout << "  1. Observed Noise Test" << endl;
        cout << "  0. Exit" << endl;

        int selection = 0;
        cout << endl << "Run example: ";
        if (!(cin >> selection))
        {
            cout << "Invalid option." << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }
        
        switch (selection)
        {
        case 1: {
            int trials;
            cout << "Trials: ";
            if (!(cin >> trials) || (trials < 1))
            {
                cout << "Invalid option." << endl;
                break;
            }
            test_noise(trials);
            break;
        }

        case 0: 
            return 0;

        default:
            cout << "Invalid option."<< endl;
            break;
        }
    }

    return 0;
}

void test_noise(int trials)
{
    NTL::xdouble trials_copy(trials);

    /* Set verbose to true for debugging. */
    bool verbose = false;

    /* Select parameters appropriate for our experiment */
    unsigned long m = 8192; // polynomial modulus n = 4096
    //unsigned long m = 16384; // polynomial modulus n = 8192
    //unsigned long m = 32768; // polynomial modulus n = 16384
    unsigned long p = 3;    // set plaintext modulus t = 3
    unsigned long s = 1;    // lower bound for number of plaintext slots
    
    /* So, we set the number of bits in the modulus chain according to HE Standard*/
    unsigned long bits;
    if (m == 4096)
    {
        bits = 54; 
    }
    else if (m == 8192)
    {
        bits = 109; 
    }
    else if (m == 16384)
    {
        bits = 218;
    }
    else
    {
        bits = 438;
    }

    /* Set other parameters to HElib defaults */
    unsigned long r = 1;    // Hensel lifting, default is 1
    unsigned long c = 2;    // columns in key switching matrix, default is 2 or 3
    unsigned long k = 80;   // security parameter, default is 80 (may not correspond to true bit security)
    
    /* Check that choice of m is ok */
    long check_m = helib::FindM(k, bits, c, p, r, s, m);
    if (check_m != m)
    {
        cout << "Could not select m = " << m << ". Using m = " << check_m << " instead." << endl;
        m = check_m;
    }

    /* Store parameters in context and construct chain of moduli */
    helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .bits(bits)
                               .c(c)
                               .build();


    // Print the context.
    context.printout();
    std::cout << std::endl;

    /* Generate keys */
    helib::SecKey secret_key(context);
    secret_key.GenSecKey();
    const helib::PubKey& public_key = secret_key;
     
    /* Construct plaintext and ciphertext objects */
    const helib::EncryptedArray& ea = context.getEA();
    helib::Ptxt<helib::BGV> plain1(context);
    helib::Ptxt<helib::BGV> plain2(context);
    helib::Ptxt<helib::BGV> plain3(context);
    helib::Ptxt<helib::BGV> plain4(context);
    helib::Ptxt<helib::BGV> plain5(context);
    helib::Ptxt<helib::BGV> plain6(context);
    helib::Ptxt<helib::BGV> plain7(context);
    helib::Ptxt<helib::BGV> plain8(context);           
    helib::Ctxt encrypted1(public_key);
    helib::Ctxt encrypted2(public_key);
    helib::Ctxt encrypted3(public_key);
    helib::Ctxt encrypted4(public_key);
    helib::Ctxt encrypted5(public_key);
    helib::Ctxt encrypted6(public_key);
    helib::Ctxt encrypted7(public_key);
    helib::Ctxt encrypted8(public_key);
    helib::Ctxt encrypted9(public_key);
    helib::Ctxt encrypted10(public_key);
    helib::Ctxt encrypted11(public_key);
    helib::Ctxt encrypted12(public_key);
    helib::Ctxt encrypted13(public_key);
    helib::Ctxt encrypted14(public_key);
    helib::Ctxt encrypted15(public_key);

    /* Holders for the running total of the observed noises in ciphertexts */
    NTL::xdouble total_fresh_observed(0);
    NTL::xdouble total_mult1_observed(0);
    NTL::xdouble total_mult2_observed(0);
    NTL::xdouble total_mult3_observed(0);

    /* Holders for the running total of the HElib estimated noises in ciphertexts */
    NTL::xdouble total_fresh_helib_est(0);
    NTL::xdouble total_mult1_helib_est(0);
    NTL::xdouble total_mult2_helib_est(0);
    NTL::xdouble total_mult3_helib_est(0);

    /* Holders for all the observed noises */
    vector<NTL::xdouble> array_fresh_observed;
    array_fresh_observed.reserve(trials);
    vector<NTL::xdouble> array_mult1_observed;
    array_mult1_observed.reserve(trials);
    vector<NTL::xdouble> array_mult2_observed;
    array_mult2_observed.reserve(trials);
    vector<NTL::xdouble> array_mult3_observed;
    array_mult3_observed.reserve(trials);

    /* Holders for all the HElib estimated noises */
    vector<NTL::xdouble> array_fresh_helib_est;
    array_fresh_helib_est.reserve(trials);
    vector<NTL::xdouble> array_mult1_helib_est;
    array_mult1_helib_est.reserve(trials);
    vector<NTL::xdouble> array_mult2_helib_est;
    array_mult2_helib_est.reserve(trials);
    vector<NTL::xdouble> array_mult3_helib_est;
    array_mult3_helib_est.reserve(trials);

    /* Gather noise data over user-specified number of trials */
    for (int i = 0; i < trials; i++)
    {
        /* Encode the values i+1, ..., i+8 into plaintexts */
        long value1 = i+1;
        long value2 = i+2;
        long value3 = i+3;
        long value4 = i+4;
        long value5 = i+5;
        long value6 = i+6;
        long value7 = i+7;
        long value8 = i+8;        
        plain1[0] = value1;
        plain2[0] = value2;
        plain3[0] = value3;
        plain4[0] = value4;
        plain5[0] = value5;
        plain6[0] = value6;
        plain7[0] = value7;
        plain8[0] = value8;

        /* Encrypt the plaintexts into ciphertexts */
        public_key.Encrypt(encrypted1, plain1);
        public_key.Encrypt(encrypted2, plain2);
        public_key.Encrypt(encrypted3, plain3);
        public_key.Encrypt(encrypted4, plain4);
        public_key.Encrypt(encrypted5, plain5);
        public_key.Encrypt(encrypted6, plain6);
        public_key.Encrypt(encrypted7, plain7);
        public_key.Encrypt(encrypted8, plain8);

        /* What is the observed noise growth at the fresh encryption of ciphertexts? */
        auto fresh_noise = get_noise_budget(encrypted1, secret_key);
        total_fresh_observed += fresh_noise;
        array_fresh_observed.push_back(fresh_noise);

        /* What is the HElib estimated noise growth at the fresh encryption of ciphertexts? */
        auto fresh_helib_est = get_helib_estimated_noise_budget(encrypted1);
        total_fresh_helib_est += fresh_helib_est;
        array_fresh_helib_est.push_back(total_fresh_helib_est);

        /*  Multiply the ciphertexts pairwise and store the output in encrypted9, ... , encrypted12 */
        encrypted9.tensorProduct(encrypted1, encrypted2);
        encrypted10.tensorProduct(encrypted3, encrypted4);
        encrypted11.tensorProduct(encrypted5, encrypted6);
        encrypted12.tensorProduct(encrypted7, encrypted8);

        /* What is the observed noise growth at the first multiplication of ciphertexts? */
        auto mult1_noise = get_noise_budget(encrypted9, secret_key);
        total_mult1_observed += mult1_noise;
        array_mult1_observed.push_back(mult1_noise);

        /* What is the HElib estimated noise growth at the first multiplication of ciphertexts? */
        auto mult1_helib_est = get_helib_estimated_noise_budget(encrypted9);
        total_mult1_helib_est += mult1_helib_est;
        array_mult1_helib_est.push_back(total_mult1_helib_est);

        /*  Multiply the ciphertexts pairwise and store the output in encrypted13, encrypted14 */
        encrypted13.tensorProduct(encrypted9, encrypted10);
        encrypted14.tensorProduct(encrypted11, encrypted12);

        /* What is the observed noise growth at the second multiplication of ciphertexts? */
        auto mult2_noise = get_noise_budget(encrypted13, secret_key);
        total_mult2_observed += mult2_noise;
        array_mult2_observed.push_back(mult2_noise);  

        /* What is the HElib estimated noise growth at the second multiplication of ciphertexts? */
        auto mult2_helib_est = get_helib_estimated_noise_budget(encrypted13);
        total_mult2_helib_est += mult2_helib_est;
        array_mult2_helib_est.push_back(total_mult2_helib_est);

        /*  Multiply the ciphertexts encrypted13 and encrypted14 and stored in encrypted15 */
        encrypted15.tensorProduct(encrypted13, encrypted14);

        /* What is the observed noise growth at the third multiplication of ciphertexts? */
        auto mult3_noise = get_noise_budget(encrypted15, secret_key);
        total_mult3_observed += mult3_noise;
        array_mult3_observed.push_back(mult3_noise);

        /* What is the HElib estimated noise growth at the third multiplication of ciphertexts? */
        auto mult3_helib_est = get_helib_estimated_noise_budget(encrypted15);
        total_mult3_helib_est += mult3_helib_est;
        array_mult3_helib_est.push_back(total_mult3_helib_est);

        if(verbose)
        {
            if(i == 2)
            {
                // Decrypt: expected result for i=2 is the product
                // ((3*4)*(5*6))*((7*8)*(9*10)) = 1814400 = 0 mod 3
                helib::Ptxt<helib::BGV> plaintext15(context);
                secret_key.Decrypt(plaintext15, encrypted15);
                std::cout << "Decrypted Result: " << plaintext15 << std::endl;
            }
        }

    }

    /* Compute the mean of the observed noises */
    auto mean_fresh_observed = total_fresh_observed / trials_copy;
    auto mean_mult1_observed = total_mult1_observed / trials_copy;
    auto mean_mult2_observed = total_mult2_observed / trials_copy;
    auto mean_mult3_observed = total_mult3_observed / trials_copy;

    /* Compute the mean of the HElib estimated noises */
    auto mean_fresh_helib_est = total_fresh_helib_est/ trials_copy;
    auto mean_mult1_helib_est = total_mult1_helib_est / trials_copy;
    auto mean_mult2_helib_est = total_mult2_helib_est / trials_copy;
    auto mean_mult3_helib_est = total_mult3_helib_est / trials_copy;

    /* Print out the results */
    cout << "After fresh encryption:" << endl;
    cout << "Mean noise budget observed: " << mean_fresh_observed  << endl;
    cout << "Mean HElib estimated noise budget: " << mean_fresh_helib_est << endl;        
    cout << endl;

    cout << "After first multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult1_observed  << endl;
    cout << "Mean HElib estimated noise budget: " << mean_mult1_helib_est << endl;            
    cout << endl;

    cout << "After second multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult2_observed  << endl;    
    cout << "Mean HElib estimated noise budget: " << mean_mult2_helib_est << endl;        
    cout << endl;

    cout << "After third multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult3_observed  << endl;    
    cout << "Mean HElib estimated noise budget: " << mean_mult3_helib_est << endl;        
    cout << endl;

}