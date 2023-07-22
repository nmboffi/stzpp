out_file = open("stencil_test.cc", "w")
for ii in [-1, 0, 1]:
    for jj in [-1, 0, 1]:
        for kk in [-1, 0, 1]:
            ijk = 13 + ii + 3*jj + 9*kk
            if (abs(ii) + abs(jj) + abs(kk) < 3):
                stenc_str = "case " + str(ijk) + ":\n"
                for qq in [1, 2, 3]:
                    for ll in [1, 2, 3]:
                        stenc_str += "\tdouble st" + str(qq) + str(ll) + " = "
                        stenc_str += "U" if ll == 1 else "V" if ll == 2 else "W"
                        stenc_str += "xp" if ii == 1 else "xm" if ii == -1 else ""
                        stenc_str += "yp" if jj == 1 else "ym" if jj == -1 else ""
                        stenc_str += "zp" if kk == 1 else "zm" if kk == -1 else ""
                        stenc_str += "_" + str(qq) + "("
                        stenc_str += "T.a11, T.a12, T.a13, T.a21, T.a22, T.a23, T.a31, T.a32, T.a33, mu, lamb, D, dx, dy, dz);\n"

                stenc_str += "\ts_ent = mat3(st11, st12, st13, st21, st22, st23, st31, st32, st33);\n\tbreak;\n\n"

                out_file.write(stenc_str)
out_file.close()
