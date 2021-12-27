from random import randint
from hashlib import blake2b
from curves import Pallas
from fields import Fp, Fq
from poly import Poly

R_INV = Fp(0x21f1c4ff1e2278d570cb2996efc89a65ac9fba6a4077fc57cf3f8e8753a769a9)
R = Fp(0x3fffffffffffffffffffffffffffffff992c350be41914ad34786d38fffffffd)
R2 = R * R

#     "iso-pallas",
A = Fp(0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b) #* R
B = Fp(1265) # * R  #Fp(Fp.modulus)
# Z = Fp(0x0f7bdb65814179b44647aef782d5cdc851f64fc4dc888857ca330bcc09ac318e)
Z = Fp(-13)
theta = Fp(0x0f7bdb65814179b44647aef782d5cdc851f64fc4dc888857ca330bcc09ac318e)


def hash_to_field(curve_id, domain_prefix, message):
    CHUNKLEN = 64
    R_IN_BYTES = 128
    tail = bytes([22 + len(curve_id) + len(domain_prefix)])
    h0 = blake2b(digest_size=CHUNKLEN, person=b'\x00' * 16)
    h0.update(b''.join([b'\x00' * R_IN_BYTES, message, bytes([0, CHUNKLEN * 2, 0]), domain_prefix,
                        b'-', curve_id, b'_XMD:BLAKE2b_SSWU_RO_', tail]))
    b_0 = h0.digest()
    h1 = blake2b(digest_size=CHUNKLEN, person=b'\x00' * 16)
    h1.update(b''.join([h0.digest(), b'\x01', domain_prefix, b'-', curve_id,
                        b'_XMD:BLAKE2b_SSWU_RO_', tail]))
    b_1 = h1.digest()
    h2 = blake2b(digest_size=CHUNKLEN, person=b'\x00' * 16)
    h2.update(bytes(a ^ b for (a, b) in zip(b_0, b_1)))
    h2.update(b''.join([b'\x02', domain_prefix, b'-', curve_id, b'_XMD:BLAKE2b_SSWU_RO_', tail]))
    b_2 = h2.digest()
    buf0 = Fp(int.from_bytes(b_1, byteorder='big'))
    buf1 = Fp(int.from_bytes(b_2, byteorder='big'))
    return buf0, buf1


def map_to_curve_simple_swu(u):
    #print("OK theta * R: ", theta * R)
    #print("OK z * R: ", Z * R)
    z_u2 = (Z * u**2)
    #print("OK z_u2 * R: ", z_u2 * R)
    # 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
    tv1 = Fp(1) / (Z**2 * u**4 + Z * u**2)
    # 2.  x1 = (-B / A) * (1 + tv1)
    x1 = ((Fp(0)-B)/A) * (Fp(1) + tv1)
    #print("OK x1 is: ", x1)
    num_x1 = (Z**2 * u**4 + Z * u**2 + Fp(1)) * B
    #print("OK num_x1: ", num_x1)  # GOOD!
    # 3.  If tv1 == 0, set x1 = B / (Z * A)
    if tv1 == Fp(0): x1 = B / (Z * A)
    # 4. gx1 = x1^3 + A * x1 + B
    gx1 = x1**3 + A * x1 + B
    #print("OK gx1: ", gx1)
    # 5.  x2 = Z * u^2 * x1
    x2 = Z * u**2 * x1
    # 6. gx2 = x2^3 + A * x2 + B
    gx2 = x2**3 + A * x2 + B
    # 7.  If is_square(gx1), set x = x1 and y = sqrt(gx1)
    if gx1.is_square():
        x = x1
        #print("if gx1.square x: ", x1)
        y = gx1.sqrt()
    # 8.  Else set x = x2 and y = sqrt(gx2)
    else:
        x = x2
        y = gx2.sqrt()
    # 9.  If sgn0(u) != sgn0(y), set y = -y
    if u.sgn0 != y.sgn0: y = Fp(0) - y
    # 10. return (x, y)
    return x, y

def is_on_iso_pallas_curve(pt):
    a = Fp(0x18354a2eb0ea8c9c49be2d7258370742b74134581a27a59f92bb4b0b657a014b)
    b = Fp(1265)
    return pt.z * pt.y**2 == pt.x**3 + a*pt.x*pt.z**2 + b*pt.z**3


def iso_map(pt):
    iso = [Fp(0x0e38e38e38e38e38e38e38e38e38e38e4081775473d8375b775f6034aaaaaaab),
           Fp(0x3509afd51872d88e267c7ffa51cf412a0f93b82ee4b994958cf863b02814fb76),
           Fp(0x17329b9ec525375398c7d7ac3d98fd13380af066cfeb6d690eb64faef37ea4f7),
           Fp(0x1c71c71c71c71c71c71c71c71c71c71c8102eea8e7b06eb6eebec06955555580),
           Fp(0x1d572e7ddc099cff5a607fcce0494a799c434ac1c96b6980c47f2ab668bcd71f),
           Fp(0x325669becaecd5d11d13bf2a7f22b105b4abf9fb9a1fc81c2aa3af1eae5b6604),
           Fp(0x1a12f684bda12f684bda12f684bda12f7642b01ad461bad25ad985b5e38e38e4),
           Fp(0x1a84d7ea8c396c47133e3ffd28e7a09507c9dc17725cca4ac67c31d8140a7dbb),
           Fp(0x3fb98ff0d2ddcadd303216cce1db9ff11765e924f745937802e2be87d225b234),
           Fp(0x025ed097b425ed097b425ed097b425ed0ac03e8e134eb3e493e53ab371c71c4f),
           Fp(0x0c02c5bcca0e6b7f0790bfb3506defb65941a3a4a97aa1b35a28279b1d1b42ae),
           Fp(0x17033d3c60c68173573b3d7f7d681310d976bbfabbc5661d4d90ab820b12320a),
           Fp(0x40000000000000000000000000000000224698fc094cf91b992d30ecfffffde5)]


if __name__ == '__main__':

    tmp1, tmp2 = hash_to_field(b'pallas', b'z.cash:test', b'Trans rights now!')
    (act_x1, act_y1) = map_to_curve_simple_swu(tmp1)
    assert act_x1 == Fp(0x05c3482fe40155e152fdc0be06c4766b67a2b3d8d9bb64ee6137382879dc2160)
    assert act_y1 == Fp(0x3825fb730c259375175ff31b94dc36dcf031b13f3116bda725f1c98717739f1f)
    (act_x2, act_y2) = map_to_curve_simple_swu(tmp2)
    assert act_x2 == Fp(0x2c6e5aa1a88cd76c8a9d436438d2993244bf7704e4f322a86d0890bd6cee28ab)
    assert act_y2 == Fp(0x0b20c46efea44d15e4828808c86a72789d54328635ba4274d8e9b48d9654f65b)

    q0 = Pallas(act_x1, act_y1, Fp(1))
    q1 = Pallas(act_x2, act_y2, Fp(1))
    r = q0 + q1
    print("r is: ", r)  # WORKS!!!!!!!!!!!!!!!

    # TODO: now need to implement ISO mult; maybe check on curve? <<-- good!

    assert is_on_iso_pallas_curve(r)



    # TODO: NEED RANDOM POINT
    # hash to curve https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/


    # a = [Fq(randint(0, 100)) for x in range(100)]
    # G = [Pallas.base() * Fq(randint(0, 100)) for x in range(100)]
    #
    # xx = sum(map(lambda a, g: g * a, a, G), start=Pallas.neutral())
    #
    # zz = map_to_curve_simple_swu(Fp(1))
    # print(zz)
    #
    # # from past curves src/pallas.rs
    # exp_x = Fp(0x010cba5957e876534af5e967c026a1856d64b071068280837913b9a5a3561505) #* R_INV
    # exp_y = Fp(0x062fc61f9cd3118e7d6e65a065ebf46a547514d6b08078e976fa6d515dcc9c81) #* R_INV
    # exp_z = Fp(0x3f86cb8c311250c3101c4e523e7793605ccff5623de1753a7c75bc9a29a73688) #* R_INV
    #
    # asdf = (exp_x / exp_z)
    # print("asdf ", asdf)
    # print("zz[0] ", zz[0])
    # assert asdf == zz[0]

    # base = Pallas.base()
    # print(base)
    # base5 = Fq(2) * base
    # print(G)
    # print(xx)

x = 0x37913b9a5a356150556d64b071068280834af5e967c026a18010cba5957e8765
z = 0x7c75bc9a29a736885ccff5623de1753a101c4e523e7793603f86cb8c311250c3



# if res = x * r, then we need to mult by r_inv
# /// R = 2^256 mod p
# const R: Fp = Fp([
#     0x34786d38fffffffd,
#     0x992c350be41914ad,
#     0xffffffffffffffff,
#     0x3fffffffffffffff,
# ]);
# R = 0x3fffffffffffffffffffffffffffffff992c350be41914ad34786d38fffffffd
# R_INV = hex(pow(R, modulus - 2, modulus))  # Fp.modlus
# R_INV = 0x21f1c4ff1e2278d570cb2996efc89a65ac9fba6a4077fc57cf3f8e8753a769a9

############# TODO: PASTA CURVES MONTGOMERY REPRESENTATION IS GIVING ME MAJOR HEADACHES
#############  DO I NEED TO TRACK MULTIPLES OF R ACROSS INTERMEDIATE CALCULATIONS? (NO??)  (YES PROBABLY!!!!!!!!!!)
#### TODO: or better YET, USE NORMAL MATH THEN MULT BY R BEFORE COMPARING TO GOLAND



## Dec 18th rethink....
##  1. Need to create a random point on the curve  TODO <=== confirm this!
##     Presumably this needs a hash to curve
##     Then we must implement matching code to Pallas
##
## implementing bls12381 pro has test vectors, con needs different curve arith