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
    h = blake2b(digest_size=CHUNKLEN, person=b'\x00' * 16)
    h.update(b'\x00' * R_IN_BYTES)
    h.update(message)
    h.update(bytes([0, CHUNKLEN * 2, 0]))
    h.update(domain_prefix)
    h.update(b'-')
    h.update(curve_id)
    h.update(b'_XMD:BLAKE2b_SSWU_RO_')
    h.update(bytes([22 + len(curve_id) + len(domain_prefix)]))
    b_0 = h.hexdigest()
    return b_0  # TODO: CORRECT, BUT WE NEED BYTEARRAY HERE






def map_to_curve_simple_swu(u):
    z_u2 = (Z * u**2) * R
    # 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
    tv1 = Fp(1) / (Z**2 * u**4 + Z * u**2)
    # 2.  x1 = (-B / A) * (1 + tv1)
    x1 = ((Fp(0)-B)/A) * (Fp(1) + tv1)
    num_x1 = (Z**2 * u**4 + Z * u**2 + Fp(1)) * B * R  ##### TODO: THIS IS OFF BY JUST A BIT!?!?!?!?!?!
    # 3.  If tv1 == 0, set x1 = B / (Z * A)
    if tv1 == Fp(0): x1 = B / (Z * A)
    # 4. gx1 = x1^3 + A * x1 + B
    gx1 = x1**3 + A * x1 + B
    # 5.  x2 = Z * u^2 * x1
    x2 = Z * u**2 * x1
    # 6. gx2 = x2^3 + A * x2 + B
    gx2 = x2**3 + A * x2 + B
    # 7.  If is_square(gx1), set x = x1 and y = sqrt(gx1)
    if gx1.is_square():
        x = x1
        y = gx1.sqrt()
    # 8.  Else set x = x2 and y = sqrt(gx2)
    else:
        x = x2
        y = gx2.sqrt()
    # 9.  If sgn0(u) != sgn0(y), set y = -y
    if u.sgn0 != y.sgn0: y = Fp(0) - y
    # 10. return (x, y)
    return x, y






if __name__ == '__main__':

    b0 = hash_to_field(b'pallas', b'z.cash:test', b'Trans rights now!')

    # TODO: NEED RANDOM POINT
    # hash to curve https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/


    a = [Fq(randint(0, 100)) for x in range(100)]
    G = [Pallas.base() * Fq(randint(0, 100)) for x in range(100)]

    xx = sum(map(lambda a, g: g * a, a, G), start=Pallas.neutral())

    zz = map_to_curve_simple_swu(Fp(1))
    print(zz)

    # from past curves src/pallas.rs
    exp_x = Fp(0x010cba5957e876534af5e967c026a1856d64b071068280837913b9a5a3561505) #* R_INV
    exp_y = Fp(0x062fc61f9cd3118e7d6e65a065ebf46a547514d6b08078e976fa6d515dcc9c81) #* R_INV
    exp_z = Fp(0x3f86cb8c311250c3101c4e523e7793605ccff5623de1753a7c75bc9a29a73688) #* R_INV

    asdf = (exp_x / exp_z)
    print("asdf ", asdf)
    print("zz[0] ", zz[0])
    assert asdf == zz[0]

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