#include <cmath>
#include "mat.hh"
#include "vec.hh"

#include <cstdio>
#include <cassert>

int main (int argc, char *argv[]) {
    // Define test matrices and test vectors
    double eles1[] = {1, 6, 3, 
                      8, 4, 2, 
                      1, 0, 3};

    double eles2[] = {5, 3, 2,
                      7, 3, 1,
                      7, 4, 3};

    mat test_mat_1 = mat(eles1);

    mat test_mat_2 = mat(eles2);
    vec test_vec_1 = vec(1, 4, 3);
    vec test_vec_2 = vec(7, 1, 4);

    // Define the answers
    // Inverse of mat_1
    double inv_eles1[] = {-1./11., 3./22., 0,
                         1./6., 0, -1./6.,
                         1./33, -1./22., 1./3.};

    mat inv_mat_1 = mat(inv_eles1);

    // mat_1 + mat_2
    double add_eles[] = {6, 9, 5,
                         15, 7, 3,
                         8, 4, 6};
    mat mat_add = mat(add_eles);

    // mat_1 - mat_2
    double sub_eles[] = {-4, 3, 1,
                         1, 1, 1,
                         -6, -4, 0};
    mat sub_mat = mat(sub_eles);
    
    // -mat_1
    double min_eles[] = {-1, -6, -3,
                         -8, -4, -2,
                         -1, 0, -3};
    mat min_mat_1 = mat(min_eles);

    // 6*mat_1
    double scale_eles[] = {6, 36, 18,
                           48, 24, 12,
                           6, 0, 18};
    mat six_mat_1 = mat(scale_eles);

    // mat_1*mat_2
    double mat_mul_eles[] = {68, 33, 17,
                             82, 44, 26,
                             26, 15, 11};
    mat mat_mul = mat(mat_mul_eles);

    // mat_1 / 3
    double div_eles[] = {1./3., 2., 1.,
                         8./3., 4./3, 2./3,
                         1./3, 0, 1.};
    mat mat_1_div_3 = mat(div_eles);

    // test_vec_1 + test_vec_2
    vec add_vec = vec(8, 5, 7);

    // test_vec_1 - test_vec_2
    vec sub_vec = vec(-6, 3, -1);

    // test_mat_1*test_vec_1
    vec mat_vec_mul = vec(34, 30, 10);

    // Better be true!
    assert((test_mat_1.inverse() - inv_mat_1).mod_sq() <= .01);
    assert(((test_mat_1 + test_mat_2) - mat_add).mod_sq() <= .01);
    assert(((test_mat_1 - test_mat_2) - sub_mat).mod_sq() <= .01);
    assert(((-test_mat_1 - min_mat_1).mod_sq() <= .01));
    assert(((6*test_mat_1 - six_mat_1).mod_sq() <= .01));
    assert((test_mat_1*test_mat_2 - mat_mul).mod_sq() <= .01);
    assert((test_mat_1/3. - mat_1_div_3).mod_sq() <= .01);
    assert(((test_vec_1 + test_vec_2 - add_vec).magnitude() <= .01));
    assert((test_vec_1 - test_vec_2 - sub_vec).magnitude() <= .01);
    assert((test_mat_1*test_vec_1 - mat_vec_mul).magnitude() <= .01);

    return 0;

}
