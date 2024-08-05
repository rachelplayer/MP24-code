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
    /* Set verbose to true for debugging. */
    bool verbose = false;

    NTL::xdouble trials_copy(trials);

    /* Select parameters appropriate for our experiment */
    unsigned long m = 4096; // polynomial modulus n = 2048
    //unsigned long m = 8192; // polynomial modulus n = 4096
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

    /* Parameter set corresponding to n = 2048 does not support modulus switching */
    bool is_not_2048 = true;
    if (m == 4096)
    {
        is_not_2048 = false;
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
    helib::Ctxt encrypted1(public_key);
    helib::Ctxt encrypted2(public_key);
    helib::Ctxt encrypted3(public_key);

    /* Holders for the running total of the observed noises in ciphertexts */
    NTL::xdouble total_fresh_observed(0);
    NTL::xdouble total_add_observed(0);
    NTL::xdouble total_mult_observed(0);
    NTL::xdouble total_modswitch_observed(0);

    /* Holders for the running total of the HElib estimated noises in ciphertexts */
    NTL::xdouble total_fresh_helib_est(0);
    NTL::xdouble total_add_helib_est(0);
    NTL::xdouble total_mult_helib_est(0);
    NTL::xdouble total_modswitch_helib_est(0);

    /* Holders for all the observed noises */
    vector<NTL::xdouble> array_fresh_observed;
    array_fresh_observed.reserve(trials);
    vector<NTL::xdouble> array_add_observed;
    array_add_observed.reserve(trials);
    vector<NTL::xdouble> array_mult_observed;
    array_mult_observed.reserve(trials);
    vector<NTL::xdouble> array_modswitch_observed;
    array_modswitch_observed.reserve(trials);

    /* Holders for all the HElib estimated noises */
    vector<NTL::xdouble> array_fresh_helib_est;
    array_fresh_helib_est.reserve(trials);
    vector<NTL::xdouble> array_add_helib_est;
    array_add_helib_est.reserve(trials);
    vector<NTL::xdouble> array_mult_helib_est;
    array_mult_helib_est.reserve(trials);
    vector<NTL::xdouble> array_modswitch_helib_est;
    array_modswitch_helib_est.reserve(trials);

    /* Gather noise data over user-specified number of trials */
    for (int i = 0; i < trials; i++)
    {
        /* Encode the values i+1, i into plaintexts */
        long value1 = i+1;
        long value2 = i;      
        plain1[0] = value1;
        plain2[0] = value2;

        /* Encrypt the plaintexts into ciphertexts */
        public_key.Encrypt(encrypted1, plain1);
        public_key.Encrypt(encrypted2, plain2);


        if(verbose)
        {
            if(i == 2)
            {

                // What actually is log q?
                NTL::xdouble check_log_q;
                check_log_q = NTL::xdouble(encrypted1.getContext().logOfProduct(encrypted1.getPrimeSet())/log(2));
                std::cout << "Log q of fresh ciphertext:" << check_log_q << std::endl;

                // Decrypt the modified ciphertext into a new plaintext
                helib::Ptxt<helib::BGV> plaintext1(context);
                secret_key.Decrypt(plaintext1, encrypted1);
                helib::Ptxt<helib::BGV> plaintext2(context);
                secret_key.Decrypt(plaintext2, encrypted2);
                std::cout << "Operation: fresh encryption" << std::endl;
                std::cout << "Decrypted Result 1: " << plaintext1 << std::endl;
                std::cout << "Decrypted Result 2: " << plaintext2 << std::endl;
            }
        }

        /* What is the observed noise growth at the fresh encryption of ciphertexts? */
        auto fresh_noise = get_noise_budget(encrypted1, secret_key);
        total_fresh_observed += fresh_noise;
        array_fresh_observed.push_back(fresh_noise);

        /* What is the HElib estimated noise growth at the fresh encryption of ciphertexts? */
        auto fresh_helib_est = get_helib_estimated_noise_budget(encrypted1);
        total_fresh_helib_est += fresh_helib_est;
        array_fresh_helib_est.push_back(total_fresh_helib_est);

        /* Compute the homomorphic addition of encrypted1 and encrypted2. Done in place, adding encrypted2 into encrypted1 */
        encrypted1 += encrypted2;

        /* What is the observed noise growth after addition? */
        auto add_noise = get_noise_budget(encrypted1, secret_key);
        total_add_observed += add_noise;
        array_add_observed.push_back(add_noise);

        /* What is the HElib estimated noise growth after addition? */
        auto add_helib_est = get_helib_estimated_noise_budget(encrypted1);
        total_add_helib_est += add_helib_est;
        array_add_helib_est.push_back(total_add_helib_est);

        /* Compute the homomorphic multiplication of encrypted1 and encrypted2 and store the output in encrypted3 */
        encrypted3.tensorProduct(encrypted1, encrypted2);

        /* What is the observed noise growth after multiplication? */
        auto mult_noise = get_noise_budget(encrypted3, secret_key);
        total_mult_observed += mult_noise;
        array_mult_observed.push_back(mult_noise);

        /* What is the HElib estimated noise growth after multiplication? */
        auto mult_helib_est = get_helib_estimated_noise_budget(encrypted3);
        total_mult_helib_est += mult_helib_est;
        array_mult_helib_est.push_back(total_mult_helib_est);

        /* Modulus switch encrypted3 down to next modulus in chain */
        if(is_not_2048)
        {

            if(i == 0)
            {
                cout << "before mod switch: bit size of q is " << encrypted3.getContext().logOfProduct(encrypted3.getPrimeSet())/log(2) << endl;
                cout << endl;
            }

            helib::IndexSet natural_primes = encrypted3.naturalPrimeSet();
            encrypted3.modDownToSet(natural_primes);

            if(i == 0)
            {
                cout << "after mod switch: bit size of q is " << encrypted3.getContext().logOfProduct(encrypted3.getPrimeSet())/log(2) << endl;
                cout << endl;
            }

            if(verbose)
            {
                if(i == 2)
                {
                    // Decrypt
                    helib::Ptxt<helib::BGV> plaintext3(context);
                    secret_key.Decrypt(plaintext3, encrypted3);
                    std::cout << "Operation: modswitch" << std::endl;
                    std::cout << "Decrypted Result: " << plaintext3 << std::endl;
                }
            }

        }

        /* What is the observed noise growth after modulus switching? */
        auto modswitch_noise = get_noise_budget(encrypted3, secret_key);
        total_modswitch_observed += modswitch_noise;
        array_modswitch_observed.push_back(modswitch_noise);
        
        /* What is the HElib estimated noise growth after modulus switching? */
        auto modswitch_helib_est = get_helib_estimated_noise_budget(encrypted3);
        total_modswitch_helib_est += modswitch_helib_est;
        array_modswitch_helib_est.push_back(total_modswitch_helib_est);

    }

    /* Compute the mean of the observed noises */
    auto mean_fresh_observed = total_fresh_observed / trials_copy;
    auto mean_add_observed = total_add_observed / trials_copy;
    auto mean_mult_observed = total_mult_observed / trials_copy;
    auto mean_modswitch_observed = total_modswitch_observed / trials_copy;

    /* Compute the mean of the HElib estimated noises */
    auto mean_fresh_helib_est = total_fresh_helib_est/ trials_copy;
    auto mean_add_helib_est = total_add_helib_est / trials_copy;
    auto mean_mult_helib_est = total_mult_helib_est / trials_copy;
    auto mean_modswitch_helib_est = total_modswitch_helib_est / trials_copy;

    /* Print out the results */
    cout << "After fresh encryption:" << endl;
    cout << "Mean noise budget observed: " << mean_fresh_observed  << endl;    
    cout << "Mean HElib estimated noise budget: " << mean_fresh_helib_est << endl;        
    cout << endl;

    cout << "After addition:" << endl;
    cout << "Mean noise budget observed: " << mean_add_observed  << endl;
    cout << "Mean HElib estimated noise budget: " << mean_add_helib_est << endl;        
    cout << endl;

    cout << "After multiplication:" << endl;
    cout << "Mean noise budget observed: " << mean_mult_observed  << endl;  
    cout << "Mean HElib estimated noise budget: " << mean_mult_helib_est << endl;        
    cout << endl;

    if(is_not_2048)
    {
        cout << "After mod switch:" << endl;
        cout << "Mean noise budget observed: " << mean_modswitch_observed  << endl;    
        cout << "Mean HElib estimated noise budget: " << mean_modswitch_helib_est << endl;        
        cout << endl;
    }


}